#********************** Licensing information ****************************
# Copyright (c) 2018 Guillaume Bareigts
# MIT License, see LICENSE
#*************************************************************************
##[ This module helps to solve the 1-D and radial Poisson-Boltzmann equation.

Examples
========

Simple solve
------------
  
Solves a PB equation for a planar geometry, a surface density of 1, using 
a fixed derivative as a left boundary condition.
Then outputs the resulting potential profile to a file.

.. code-block:: nim

    import polypbrenpkg / pbsolv

    let
      ds=0.1
      zv = [1.0, -1.0]
      conc = [0.5, 0.5]
    var eq = initPBEquation(smin = 0.0, smax = 10.0, sbin = 0.1,
                            kind=PlanePBE, zv = zv, concv=conc, 
                            bvL=BCder, bcR=BCDer, left = 1.0, right = 0.0, 
                            charge = 1.0, lambdab = 1.0)
    # initialize the potential profile to 1.
    eq.initialize(1.0)
    # solves the equation.
    let number_iterations = eq.solve(niter = 5000, tol = 0.00001)
    var fout = open("phi.dat",fmWrite)
    for i,x in eq.mesh:
      fout.writeLine(x, ' ',eq.phi[i])
    fout.close

With a condition
----------------

Solves a PB equation for radial geometry with a condition that
should converge to 0 after enough iterations.
Here it is a simplified pH condition.

.. code-block:: nim

    import polypbrenpkg / pbsolv
  
    let
      pKa = 7.7
      siteDens = 5.5  # density of sites
      dqSite = -1.0    # charge difference of a site
      pHref = 9.0     # imposed pH 

    proc pHcond(e: PBEquation): float =
      let
        alpha = abs(sigma(e))/siteDens   # sigma() computes the charge density (defined in this module)
        pH = pKa + (ln(alpha/(1.0-alpha)) + dqSite*e.phi[0])/ln(10.0)
      result = pH - pHref

    var eq = initPBEquation(smin = 8.0, smax = 30.0, sbin = 0.1,
                            kind=RadialPBE, zv = [1.0,-1.0], concv=[6.022e-3,6.022e-3], 
                            lambdab = 0.7105)
    # initialize the potential profile to -1.
    eq.initialize(-1.0)
    # solves the equation.
    let number_iterations = eq.solveByCondition(condition = pHcond, minVal = -10.0, maxVal = 0.0,
                                                  nitermax = 50, niterin = 5000, tolIn = 1.0e-4)
]##

import strutils, math
import macros

import pbsolv_common
export pbsolv_common

when defined(pbsolvSimdTesting):
  import cpuinfo_testing as cpuinfo
else:
  import cpuinfo

import pbsolv_sse2
import pbsolv_avx2
when defined(pbsolvNoAVX512):
  import pbsolv_avx2 as pbsolv_avx512
else:
  import pbsolv_avx512

type
  AlignedArray*[T] = object
    len*: int
    data: ptr UncheckedArray[T]
    align: int
    storage: seq[T]

# proc `=copy`*[T](a: var AlignedArray[T], b: AlignedArray[T]) =
#   if a.data != nil:
#     alignedDealloc(a.data, a.align)
#   a.data = cast[ptr UncheckedArray[T]](alignedAlloc0(b.len, b.align))
#   a.len = b.len
#   a.align = b.align
#   for i in 0..<b.len:
#     a.data[i] = b.data[i]

# proc `=destroy`*[T](a: var AlignedArray[T]) =
#   if a.data != nil:
#     alignedDealloc(a.data, a.align)

proc initAlignedArray*[T](len: int, align: int): AlignedArray[T] =
  result.len = len
  result.align = align
  result.storage.setLen(len + align - 1)
  let offset = align - ((cast[int](addr result.storage[1])) and (align - 1)) # we want 2nd elem to be aligned because loop starts with it
  result.data = cast[ptr UncheckedArray[T]](addr result.storage[offset div sizeof(T)])
  #doassert (cast[int](addr result.data[0]) and (align-1)) == 0
  doassert (cast[int](addr result.data[1]) and (align-1)) == 0


proc checkBounds(len, idx: int) =
  if idx < 0 or idx >= len:
    raise newException(when (NimMajor, NimMinor) <= (1,4): IndexError else: IndexDefect, 
      "Index $1 not in 0..$2" % [$idx, $(len-1)])

template `[]`*[T](a: AlignedArray[T], idx: int): T = 
  when compileOption("boundchecks"):
    checkBounds(a.len, idx)
  a.data[idx]

template `[]=`*[T](a: AlignedArray[T], idx: int, val: T) = 
  when compileOption("boundchecks"): 
    checkBounds(a.len, idx)
  a.data[idx] = val

func toSeq*[T](a: AlignedArray[T]): seq[T] =
  result.setLen(a.len)
  for i in 0..<a.len:
    result[i] = a.data[i]

iterator items*[T](a: AlignedArray[T]): T =
  for i in 0..<a.len:
    yield a[i]

iterator mitems*[T](a: AlignedArray[T]): var T =
  for i in 0..<a.len:
    yield a[i]

iterator pairs*[T](a: AlignedArray[T]): (int, T) =
  for i in 0..<a.len:
    yield (i, a[i])

iterator mpairs*[T](a: AlignedArray[T]): (int, var T) =
  for i in 0..<a.len:
    yield (i, a[i])


func high*(a: AlignedArray): int {.inline.} = a.len-1

template toOpenArray*[T](a: AlignedArray[T]): openArray[T] = toOpenArray(a.data, 0, a.len-1)

type
  PBEquation* = object
    pars: PBEquationPars
    grid: int
    h, h2: float
    phi*: AlignedArray[float]
    phi2: AlignedArray[float]

# proc phi*(pbe : PBEquation): seq()

proc defMesh*(pbe: var PBEquation, smin, smax,sbin: float, reDefSbin = true, gridStep = 2) =
  ## Defines (or redefines) the space grid for equation ``pbe``.
  ## ``smin`` is the lower bound, ``smax`` the higher, ``sbin``
  ## the space interval.
  ## If ``reDefSbin`` is set to true, then ``sbin`` is redefined
  ## to match the [``sbin``, ``smax``] interval and ``gridStep``.
  ## Dimensions of potential tables : N = int((smin-smax)/pas)
  ## raises ``ValueError`` if the grid is too small (number of points < 3)
  let vecSize = 64 div sizeof(float)  # max simd vector size
  var
    nv = (int((smax-smin)/sbin) ) div (vecSize*gridStep)
    n2 = nv*vecSize + 2
    n = gridStep * (n2 - 2) + 2
  if n2 < 3: raise newException(ValueError, "Grid too small (#points < 3*<vector size>): min=" & $smin & " max=" & $smax & "bin=" & $sbin)
  if reDefSbin:
    pbe.h = (smax-smin)/(float(n-1))
  else:
    pbe.h = sbin
  pbe.h2 = gridStep.float * pbe.h
  pbe.pars.s0 = smin
  pbe.phi = initAlignedArray[float](n, 64) #newSeq[float](n)
  pbe.phi2 = initAlignedArray[float](n2, 64)
  pbe.pars.N = n
  pbe.grid = gridStep

proc `charge=`*(e: var PBEquation, charge: float) =
  ## sets the surface charge density of the left (inner) surface
  e.pars.f0 = charge*(4*PI*e.pars.lb)

proc `rho=`*(e: var PBEquation, rho: float) =
  ## sets the density of an uniform background charge
  e.pars.rho = 4.0*PI*e.pars.lb*rho

proc `lambdab=`*(e: var PBEquation, lambdab: float) =
  ## changes the bjerrum length
  if e.pars.lb != 0.0:
    e.pars.rho = e.pars.rho/e.pars.lb*lambdab
    e.pars.f0 = e.pars.f0/e.pars.lb*lambdab
  e.pars.lb = lambdab

proc defBoundCond*(e: var PBEquation, bcL,bcR: BoundaryCondition, left = 0.0, right = 0.0) =
  ## Defines the boundary conditions : left (or inner radius) (type is ``bcL``, 
  ## value is ``left``),
  ## and right (or outer radius) (type is ``bcL``, value is ``left``) conditions.
  e.pars.bcL = bcL
  case bcL
  of BCval: e.pars.phiL = left
  of BCder: e.pars.phiLd = left
  e.pars.bcR = bcR
  case bcR
  of BCval: e.pars.phiR = right
  of BCder: e.pars.phiRd = right

proc initPBEquation*(smin, smax, sbin: float, zv, concv: openArray[float], kind = PlanePBE,
                    bcL=BCval, bcR=BCder, left=0.0, right = 0.0, lambdab = 0.715,
                    charge = 0.0 ,rho= 0.0, gridStep = 2): PBEquation =
  ## Returns a new PBEquation object.
  ## Parameters:
  ##   - ``smin``: minimum of radius/x coordinate
  ##   - ``smax``: maximum of radius/x coordinate
  ##   - ``sbin``: requested separation between points
  ##   - ``zv``: charges of ion species
  ##   - ``conv``: bulk concentration of ion species
  ##   - ``kind``: defines a planar geometry or a spherical geometry
  ##   - ``bcL``: type of boundary condition on the "left" (inner radius)
  ##   - ``bcR``: type of boundary condition on the "right" (inner radius)
  ##   - ``left``: value of the left boundary condition
  ##   - ``right``: value of the right boundary condition
  ##   - ``charge``: surface charge density of the left (inner) surface
  ##   - ``rho``: charge density of the background charge

  defMesh(result, smin,smax,sbin, gridStep = gridStep)
  defBoundCond(result, bcL, bcR, left, right)

  result.pars.kind = kind
  result.pars.f0 = charge*(4*PI*lambdab)
  result.pars.lb = lambdab
  result.pars.rho = 4.0*PI*lambdab*rho
  assert(zv.len >= 1 and concv.len >= 1)
  result.pars.nions = min(min(zv.len, concv.len),NionsMax)
  var truc=0.0
  for i in 0..<result.pars.nions:
    result.pars.z[i] = zv[i]
    result.pars.fac[i] = 4.0*PI*lambdab*zv[i]*concv[i]
    truc += zv[i]*result.pars.fac[i]
  if(result.pars.nions == 2):
    if zv[0] == -zv[1] : result.pars.saltKind = skIsoval
    else: result.pars.saltKind = skGeneric

proc initialize*(e: var PBEquation, f: proc (x: float): float {.closure.}) =
  ## Initialize the potential with a function.
  var x = e.pars.s0
  for phi in e.phi.mitems:
    phi = f(x)
    x+=e.h

proc initialize*(e: var PBEquation, x: float) =
  ## Initialize the potential with an uniform value.
  for phi in e.phi.mitems:
    phi = x

proc mesh*(e: var PBEquation): seq[float] =
  ## gives the grid used by ``e`` as a sequence of x coordinates (or radial coordinates).
  result = newSeq[float](e.phi.len)
  for i in 0..e.phi.high: result[i] = e.pars.s0 + float(i)*e.h

proc chargeTerm(e: PBEquation,phi: float): float {.noSideEffect.} =
  result = e.pars.rho
  for i in 0..<e.pars.nions:
    result += e.pars.fac[i]*exp(-e.pars.z[i]*phi)

proc ztotAdim(e: PBEquation): float {.noSideEffect.} =
  let
    h = e.h
  assert(e.phi.len >= 3)
  var sum = 0.0
  var s= e.pars.s0
  case e.pars.kind
  of RadialPBE:
    sum = 0.5*s*s*chargeTerm(e,e.phi[0])
    for i in 1..<e.phi.high:
      s+=h
      sum+=s*s*chargeTerm(e,e.phi[i])
    s+=h
    sum += 0.5*s*s*chargeTerm(e,e.phi[e.phi.high])
    result = e.pars.f0 + (h*sum)/(e.pars.s0*e.pars.s0)
  of PlanePBE:
    sum = 0.5*chargeTerm(e,e.phi[0])
    for i in 1..<e.phi.high:
      sum+=chargeTerm(e,e.phi[i])
    sum += chargeTerm(e,e.phi[e.phi.high])
    result = e.pars.f0 + h*sum
    
proc sigma*(e: PBEquation): float {.inline.} = (-ztotAdim(e)+e.pars.f0)/(4*PI*e.pars.lb)
  ## Computes the inner (ion+background) surface charge density of the system

proc sigmaCum*(e: PBEquation): seq[float] {.noSideEffect.} =
  ## Computes the cumulative (ion+background) surface charge density of the system
  result = newSeq[float](e.phi.len)
  let
    h = e.h
  assert(e.phi.len >= 3)
  var sum = 0.0
  var s= e.pars.s0
  result[0]=0.0
  case e.pars.kind
  of RadialPBE:
    var prev = s*s*chargeTerm(e,e.phi[0])
    for i in 1..e.phi.high:
      s+=h
      let next=s*s*chargeTerm(e,e.phi[i])
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= (0.5*h)/(e.pars.s0*e.pars.s0)
  of PlanePBE:
    var prev = chargeTerm(e,e.phi[0])
    for i in 1..e.phi.high:
      s+=h
      let next=chargeTerm(e,e.phi[i])
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= (0.5*h)
  for r in result.mitems: r /= (4*PI*e.pars.lb)
  for r in result.mitems: r += e.sigma

proc sigmaCumSansBack*(e: PBEquation): seq[float] {.noSideEffect.} =
  ## Computes the cumulative (ions only) surface charge density of the system
  result = newSeq[float](e.phi.len)
  let
    h = e.h
  assert(e.phi.len >= 3)
  var sum = 0.0
  var s= e.pars.s0
  result[0]=0.0
  case e.pars.kind
  of RadialPBE:
    var prev = s*s*(chargeTerm(e,e.phi[0])-e.pars.rho)
    for i in 1..e.phi.high:
      s+=h
      let next=s*s*(chargeTerm(e,e.phi[i])-e.pars.rho)
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= (0.5*h)/(e.pars.s0*e.pars.s0)
  of PlanePBE:
    var prev = (chargeTerm(e,e.phi[0])-e.pars.rho)
    for i in 1.. e.phi.high:
      s+=h
      let next=(chargeTerm(e,e.phi[i])-e.pars.rho)
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= -(0.5*h)
  for r in result.mitems: r /= (4*PI*e.pars.lb)
  for r in result.mitems: r += e.sigma


proc laplacian(e: PBEquation, i: int): float =
  assert(e.pars.N>=3)
  assert(i>0 and i<e.phi.high)
  let
    h = e.h
    s = i.float*h
    shi = 1.0/(s*h)
    hi2 = 1.0/(h*h)
    phiprev = e.phi[i-1]
    phi = e.phi[i]
    phinext = e.phi[i+1]
  case e.pars.kind
  of PlanePBE: result = (phinext - 2.0*phi + phiprev)*hi2
  of RadialPBE: result = (phinext - 2.0*phi + phiprev)*hi2 + (phinext - phiprev)*shi


func ngs_iter(RadialPlac: static[EqModelkind], saltKind: static[SaltKind], eq: var PBEquationPars,
    phi: var openArray[float], h: float): float {.inline.} =
  {.noSideEffect.}:
    when not defined(amd64) or defined(pbsolvNoSimd):
      result = ngs_iter_novec(RadialPlac, saltKind, eq, phi, h) # may be 10% than fake vectorized one
    else:
      result = 
        if cpuinfo.hasAVX512f:
          pbsolv_avx512.ngs_iter_vec(RadialPlac, saltKind, eq, phi, h)
        elif cpuinfo.hasAVX2:
          pbsolv_avx2.ngs_iter_vec(RadialPlac, saltKind, eq, phi, h)
        elif cpuinfo.hasSSE2:
          pbsolv_sse2.ngs_iter_vec(RadialPlac, saltKind, eq, phi, h)
        else:
          ngs_iter_novec(RadialPlac, saltKind, eq, phi, h)

proc residual*(e : PBEquation): float =
  result = 0.0
  for i in 1..<e.phi.high:
    let res = (laplacian(e, i) + chargeTerm(e, e.phi[i]))
    result += res*res
  result = sqrt(result/(e.phi.high-2).float)

template branchDispatchEqkind(kind: EqModelKind, saltKind: SaltKind, body) =
  const
    eqKind {.inject.} = kind
    eqSaltKind {.inject.} = saltKind
  body

macro dispatchEqKind(eqKind: EqModelKind, eqSaltKind: SaltKind, body: untyped): untyped =
  result = nnkCaseStmt.newTree(eqKind)
  for kind in EqModelKind:
    let nkind = newIdentNode($kind)
    let branch = nnkCaseStmt.newTree(eqSaltKind)
    for saltKind in SaltKind:
      let nskind = newIdentNode($saltKind)
      branch.add nnkOfBranch.newTree(nskind, getAst(branchDispatchEqkind(nkind, nskind, body)))
    result.add  nnkOfBranch.newTree(nkind, branch)
  when defined(pbsolvDebug):
    echo repr(result)
  
proc solve*(eq: var PBEquation, niter = 10000, tol = 0.00001): int {.noSideEffect.} =
  ## Solves the Equation for the boundary conditions previously set, in
  ## a maximum of ``niter`` iterations, for a tolerance of ``tol`` on the residue.
  ## Returns the number of iterations effectively computed.
  let tol2 = tol*tol
  var resi = tol2+1.0
  if eq.phi.len < 3: return 0
  if eq.pars.bcL == BCval: eq.phi[0] = eq.pars.phiL
  if eq.pars.bcR == BCval: eq.phi[eq.phi.high] = eq.pars.phiR
  var count = niter
  dispatchEqKind(eq.pars.kind, eq.pars.saltKind):
    eq.phi2[0] = eq.phi[0]
    for i in 1..<eq.phi2.len-1:
      eq.phi[i] = eq.phi[i*eq.grid]
    eq.phi2[eq.phi2.len-1] = eq.phi[eq.phi.len-1]
    for step in 1..niter:
      resi = ngs_iter(eqKind, eqSaltKind, eq.pars, eq.phi2.toOpenArray, eq.h2)
      if resi<=tol2:
        count = step
        break
    for i in 1..<eq.phi2.len-1:
      let delt = (eq.phi2[i+1] - eq.phi2[i])
      for j in 0..<eq.grid:
        eq.phi[i*eq.grid+j] = eq.phi2[i] + delt * (j/eq.grid)
    eq.phi[eq.phi2.len-1] = eq.phi[eq.phi2.len-1]

    for step in 1..niter:
      resi = ngs_iter(eqKind, eqSaltKind, eq.pars, eq.phi.toOpenArray, eq.h)
      if resi<=tol2:
        count = step
        break
  result = count
  when defined(pbsolvDebug):
    echo "Iters = ", count," ,Residual = ", eq.residual


proc solveByNeutralization*(eq: var PBEquation, tolF = 0.001, tolIn = 0.00001, nitermax = 500, niterin = 10000, minVal = 0.0, maxVal = 10.0): int =
  ## Solves the Equation by matching the input ``eq.charge`` (placed on left or inner surface) with
  ## the calculated ionic + outer static charge iteratively, with
  ## a maximum of ``nitermax`` iterations, with ``niterin`` steps
  ## for each iteration. Returns the number of iterations *effectively*
  ## computed.
  assert(eq.phi.len >= 3)
  let tolc = abs(tolF*eq.pars.f0)
  result = 0
  var
    phi0min = minVal
    phi0max= maxVal
  if phi0min > phi0max: swap phi0min, phi0max
  var compzin: float
  eq.pars.bcL= BCval
  eq.pars.phiL = phi0max
  discard solve(eq, niter = niterin, tol = tolIn)
  var compzoldin = ztotAdim(eq)
  var phi0 = 0.5*(phi0max + phi0min)
  for iter in 0..<nitermax:
    eq.pars.phiL = phi0
    let nitin = solve(eq, niter = niterin, tol = tolIn)
    if nitin >= niterin:
      stderr.writeLine("solver::Note: Reached max number of inner iteration ($1)" % $nitin)
    compzin = eq.ztotAdim()
    when defined(pbsolvDebug):
      let deriv = (eq.phi[0]-eq.phi[1])/eq.h
      echo "solveByNeutralization :: Ztot = ",compzin, " (",
          phi0, ',',deriv,',',eq.pars.f0,',', phi0min,',',phi0max

    inc(result)
    # Dichoto
    if compzin*compzoldin > 0.0: phi0max=phi0 #; compzoldin=compzin
    else: phi0min = phi0
    phi0 = 0.5*(phi0max + phi0min)
    if abs(compzin) <= tolc: break

proc solveByCondition*(eq: var PBEquation, condition: proc (e: PBEquation): float {.closure.} , tolF = 0.001, tolIn = 0.00001,
                       nitermax = 500, niterin = 10000, minVal = 0.0, maxVal = 10.0): int =
  ## Solves the Equation by nullifying iteratively the conditional function  ``condition``
  ## with a maximum of ``nitermax`` iterations, with ``niterin`` steps
  ## for each iteration. Returns the number of iterations actually
  ## computed.
  assert(eq.phi.len >= 3)
  let tolc = abs(tolF)  # Input must be fine-tuned ?
  result = 0
  var
    phi0min = minVal
    phi0max= maxVal
  if phi0min > phi0max: swap phi0min, phi0max
  var compzin: float
  eq.pars.bcL= BCval
  eq.pars.phiL = phi0max
  discard solve(eq, niter = niterin, tol = tolIn)
  var compzoldin = condition(eq)
  var phi0 = 0.5*(phi0max + phi0min)
  for iter in 0..<nitermax:
    eq.pars.phiL = phi0
    let nitin = solve(eq, niter = niterin, tol = tolIn)
    if nitin >= niterin:
      stderr.writeLine("solver::Note: Reached max number of inner iteration ($1)" % $nitin)
    compzin = condition(eq)
    when defined(pbsolvDebug):
      let deriv = (eq.phi[0]-eq.phi[1])/eq.h
      echo "solveByCondition :: Condition = ",compzin, " ChargeD = ", sigma(eq), " Phi0=", phi0, " Phi0m=", phi0min," Phi0M=", phi0max
          #" (",phi0, ',',deriv,',',eq.f0,')'
    inc(result)
    # Dichoto
    if compzin*compzoldin > 0.0: phi0max=phi0; compzoldin=compzin
    else: phi0min = phi0
    phi0 = 0.5*(phi0max + phi0min)
    if abs(compzin) <= tolc: break
    

    
proc residuals*(e : PBEquation): seq[float] =
  result = newSeq[float](e.phi.len)
  for i in 1..<e.phi.high:
    result[i] = laplacian(e, i) + chargeTerm(e, e.phi[i])



when isMainModule:
  let ds=0.1
  let zv = [1.0, -1.0]
  let conc = [0.5, 0.5]
  var eq = initPBEquation(smin = 0.0, smax = 10.0, sbin = 0.1, zv = zv, concv=conc, kind=PlanePBE ,left = 1.0, charge = 1.0, lambdab = 1.0)
  #eq.initialize(proc (x:float):float = exp(-x))
  eq.initialize(1.0)
  echo "Nb ofs its = ",eq.solve(niter = 5000, tol = 0.00001)
  let res = eq.residuals
  var fin = open("phi.dat",fmWrite)
  for i,x  in eq.mesh:
    fin.writeLine(x, ' ',eq.phi[i], ' ', res[i])
  fin.close

  echo ""
  var eq2 = initPBEquation(smin = 1.0, smax = 2.0, sbin = 0.005, zv = [-1.0], concv=[1.0], left = 1.0, charge = 1.0, lambdab = 1.0)
  eq2.initialize(1.0)
  echo "Nb ofs its = ", eq2.solveByNeutralization(nitermax = 500, niterin = 10000, maxVal = 10.0)
  var fin2 = open("phi2.dat",fmWrite)
  for i,phi in eq2.phi:
    fin2.writeLine(i.float*ds, ' ',phi)
  fin2.close
