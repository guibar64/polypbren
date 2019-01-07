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

import sequtils, strutils, math
import macros

const NionsMax = 8
type
  EqModelKind* = enum
    PlanePBE, RadialPBE
  SaltKind = enum
    skGeneric, skIsoVal
  BoundaryCondition* = enum
    BCval,BCder
  PBEquation* = object of RootObj
    kind*: EqModelKind
    saltKind: SaltKind
    bcL,bcR: BoundaryCondition
    N, nions: int
    s0, h, f0: float
    phiL, phiLd, phiR,phiRd: float
    lb, rho: float
    fac,z : array[Nionsmax, float]
    phi*: seq[float]

proc defMesh*(pbe: var PBEquation, smin, smax,sbin: float, reDefSbin = true) =
  ## Defines (or redefines) the space grid for equation ``pbe``.
  ## ``smin`` is the lower bound, ``smax`` the higher, ``sbin``
  ## the space interval.
  ## If ``reDefSbin`` is set to true, then ``sbin`` is redefined
  ## to match the [``sbin``, ``smax``] interval.
  ## Dimensions of potential tables : N = int((smin-smax)/pas)
  ## raises ``ValueError`` if the grid is too small (number of points < 3)
  var N = int((smax-smin)/sbin) + 1
  if N < 3: raise newException(ValueError, "Grid too small (#points < 3): min=" & $smin & " max=" & $smax & "bin=" & $sbin)
  if reDefSbin:
    pbe.h = (smax-smin)/(float(N-1))
  else:
    pbe.h = sbin
  pbe.s0 = smin
  pbe.phi = newSeq[float](N)
  pbe.N = N

proc `charge=`*(e: var PBEquation, charge: float) =
  ## sets the surface charge density of the left (inner) surface
  e.f0 = charge*(4*PI*e.lb)

proc `rho=`*(e: var PBEquation, rho: float) =
  ## sets the density of an uniform background charge
  e.rho = 4.0*PI*e.lb*rho

proc `lambdab=`*(e: var PBEquation, lambdab: float) =
  ## changes the bjerrum length
  if e.lb != 0.0:
    e.rho = e.rho/e.lb*lambdab
    e.f0 = e.f0/e.lb*lambdab
  e.lb = lambdab

proc defBoundCond*(e: var PBEquation, bcL,bcR: BoundaryCondition, left = 0.0, right = 0.0) =
  ## Defines the boundary conditions : left (or inner radius) (type is ``bcL``, 
  ## value is ``left``),
  ## and right (or outer radius) (type is ``bcL``, value is ``left``) conditions.
  e.bcL = bcL
  case bcL
  of BCval: e.phiL = left
  of BCder: e.phiLd = left
  e.bcR = bcR
  case bcR
  of BCval: e.phiR = right
  of BCder: e.phiRd = right

proc initPBEquation*(smin, smax, sbin: float, zv, concv: openArray[float], kind = PlanePBE,
                    bcL=BCval, bcR=BCder, left=0.0, right = 0.0, lambdab = 0.715,
                    charge = 0.0 ,rho= 0.0): PBEquation =
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

  defMesh(result, smin,smax,sbin)
  defBoundCond(result, bcL, bcR, left, right)

  result.kind = kind
  result.f0 = charge*(4*PI*lambdab)
  result.lb = lambdab
  result.rho = 4.0*PI*lambdab*rho
  assert(zv.len >= 1 and concv.len >= 1)
  result.nions = min(min(zv.len, concv.len),NionsMax)
  var truc=0.0
  for i in 0..<result.nions:
    result.z[i] = zv[i]
    result.fac[i] = 4.0*PI*lambdab*zv[i]*concv[i]
    truc += zv[i]*result.fac[i]
  if(result.nions == 2):
    if zv[0] == -zv[1] : result.saltKind = skIsoval
    else: result.saltKind = skGeneric

proc initialize*(e: var PBEquation, f: proc (x: float): float {.closure.}) =
  ## Initialize the potential with a function.
  var x = e.s0
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
  for i in 0..e.phi.high: result[i] = e.s0 + float(i)*e.h
      
proc chargeTerm(e: PBEquation,phi: float): float {.noSideEffect.} =
  result = e.rho
  for i in 0..<e.nions:
    result += e.fac[i]*exp(-e.z[i]*phi)

proc chargeAndDerivTerms(e: PBEquation,phi: float): tuple[ch: float,der: float] {.noSideEffect.} =
  result.ch = e.rho
  for i in 0..<e.nions:
    let
      z = e.z[i]
      fl = e.fac[i]*exp(-z*phi)
    result.ch += fl
    result.der -= z*fl

proc chargeAndDerivTermsIsoVal(e: PBEquation,phi: float): tuple[ch: float,der: float] {.noSideEffect.} =
  let z=e.z[0]
  let dump=e.fac[0]
  let dump2=exp(-z*phi)
  let dump3=1.0/dump2
  result.ch = e.rho + dump*(dump2-dump3)
  result.der = -z*dump*(dump2+dump3)

proc ztotAdim(e: PBEquation): float {.noSideEffect.} =
  let
    h = e.h
  assert(e.phi.len >= 3)
  var sum = 0.0
  var s= e.s0
  case e.kind
  of RadialPBE:
    sum = 0.5*s*s*chargeTerm(e,e.phi[0])
    for i in 1..<e.phi.high:
      s+=h
      sum+=s*s*chargeTerm(e,e.phi[i])
    s+=h
    sum += 0.5*s*s*chargeTerm(e,e.phi[e.phi.high])
    result = e.f0 + (h*sum)/(e.s0*e.s0)
  of PlanePBE:
    sum = 0.5*chargeTerm(e,e.phi[0])
    for i in 1..<e.phi.high:
      sum+=chargeTerm(e,e.phi[i])
    sum += chargeTerm(e,e.phi[e.phi.high])
    result = e.f0 + h*sum
    
proc sigma*(e: PBEquation): float {.inline.} = (-ztotAdim(e)+e.f0)/(4*PI*e.lb)
  ## Computes the inner (ion+background) surface charge density of the system

proc sigmaCum*(e: PBEquation): seq[float] {.noSideEffect.} =
  ## Computes the cumulative (ion+background) surface charge density of the system
  result = newSeq[float](e.phi.len)
  let
    h = e.h
  assert(e.phi.len >= 3)
  var sum = 0.0
  var s= e.s0
  result[0]=0.0
  case e.kind
  of RadialPBE:
    var prev = s*s*chargeTerm(e,e.phi[0])
    for i in 1.. e.phi.high:
      s+=h
      let next=s*s*chargeTerm(e,e.phi[i])
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= (0.5*h)/(e.s0*e.s0)
  of PlanePBE:
    var prev = chargeTerm(e,e.phi[0])
    for i in 1.. e.phi.high:
      s+=h
      let next=chargeTerm(e,e.phi[i])
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= (0.5*h)
  for r in result.mitems: r /= (4*PI*e.lb)
  for r in result.mitems: r += e.sigma

proc sigmaCumSansBack*(e: PBEquation): seq[float] {.noSideEffect.} =
  ## Computes the cumulative (ions only) surface charge density of the system
  result = newSeq[float](e.phi.len)
  let
    h = e.h
  assert(e.phi.len >= 3)
  var sum = 0.0
  var s= e.s0
  result[0]=0.0
  case e.kind
  of RadialPBE:
    var prev = s*s*(chargeTerm(e,e.phi[0])-e.rho)
    for i in 1.. e.phi.high:
      s+=h
      let next=s*s*(chargeTerm(e,e.phi[i])-e.rho)
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= (0.5*h)/(e.s0*e.s0)
  of PlanePBE:
    var prev = (chargeTerm(e,e.phi[0])-e.rho)
    for i in 1.. e.phi.high:
      s+=h
      let next=(chargeTerm(e,e.phi[i])-e.rho)
      sum+=prev+next
      result[i] = sum
      prev=next
    for r in result.mitems: r *= -(0.5*h)
  for r in result.mitems: r /= (4*PI*e.lb)
  for r in result.mitems: r += e.sigma


proc laplacian(e: PBEquation, i: int): float =
  assert(e.N>=3)
  assert(i>0 and i<e.phi.high)
  let
    h = e.h
    s = i.float*h
    shi = 1.0/(s*h)
    hi2 = 1.0/(h*h)
    phiprev = e.phi[i-1]
    phi = e.phi[i]
    phinext = e.phi[i+1]
  case e.kind
  of PlanePBE: result = (phinext - 2.0*phi + phiprev)*hi2
  of RadialPBE: result = (phinext - 2.0*phi + phiprev)*hi2 + (phinext - phiprev)*shi

proc ngs_iter(RadialPlac: static[EqModelkind], saltKind: static[SaltKind], eq: var PBEquation): float {.noSideEffect, inline.} =
# radial PBNewton Gauss-Seidel iteration
  let
    hh=eq.h
    hh2=hh*hh
    hi2=1.0/(eq.h*eq.h)
    moins2hi2= -2.0*hi2
    ss0=eq.s0
  result = 0.0
  var sh=hh*ss0
  for i in 1.. eq.phi.high-1:
    sh += hh2
    let (cht, numt) =
      when saltKind == skIsoVal: chargeAndDerivTermsIsoVal(eq, eq.phi[i])
      else: chargeAndDerivTerms(eq, eq.phi[i])
    let dlduo = moins2hi2 + numt
    var dent: float
    when RadialPlac == RadialPBE:
      let shi=1.0/sh
      dent = (eq.phi[i+1] - 2.0*eq.phi[i] + eq.phi[i-1])*hi2 +
           (eq.phi[i+1] - eq.phi[i-1])*shi + cht
    else:
      dent = (eq.phi[i+1] - 2.0*eq.phi[i] + eq.phi[i-1])*hi2 + cht
    eq.phi[i] -= dent/dlduo
    result += dent*dent
  result /= float(eq.phi.high-1)
  # BCs
  let nlast=eq.phi.high
  if eq.bcR == BCder: eq.phi[nlast] = eq.phi[nlast-1] + eq.phiRd*eq.h
  if eq.bcL == BCder: eq.phi[0] = eq.phi[1] - eq.phiLd*eq.h

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
  if eq.N < 3: return 0
  if eq.bcL == BCval: eq.phi[0] = eq.phiL
  if eq.bcR == BCval: eq.phi[eq.phi.high] = eq.phiR
  var count = niter
  dispatchEqKind(eq.kind, eq.saltKind):
    for step in 1..niter:
      resi = ngs_iter(eqKind, eqSaltKind, eq)
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
  assert(eq.N >= 3)
  let tolc = abs(tolF*eq.f0)
  result = 0
  var
    phi0min = minVal
    phi0max= maxVal
  if phi0min > phi0max: swap phi0min, phi0max
  var compzin: float
  eq.bcL= BCval
  eq.phiL = phi0max
  discard solve(eq, niter = niterin, tol = tolIn)
  var compzoldin = ztotAdim(eq)
  var phi0 = 0.5*(phi0max + phi0min)
  for iter in 0..<nitermax:
    eq.phiL = phi0
    let nitin = solve(eq, niter = niterin, tol = tolIn)
    if nitin >= niterin:
      stderr.writeLine("solver::Note: Reached max number of inner iteration ($1)" % $nitin)
    compzin = eq.ztotAdim()
    when defined(pbsolvDebug):
      let deriv = (eq.phi[0]-eq.phi[1])/eq.h
      echo "solveByNeutralization :: Ztot = ",compzin, " (",
          phi0, ',',deriv,',',eq.f0,',', phi0min,',',phi0max

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
  assert(eq.N >= 3)
  let tolc = abs(tolF)  # Input must be fine-tuned ?
  result = 0
  var
    phi0min = minVal
    phi0max= maxVal
  if phi0min > phi0max: swap phi0min, phi0max
  var compzin: float
  eq.bcL= BCval
  eq.phiL = phi0max
  discard solve(eq, niter = niterin, tol = tolIn)
  var compzoldin = condition(eq)
  var phi0 = 0.5*(phi0max + phi0min)
  for iter in 0..<nitermax:
    eq.phiL = phi0
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
  result = newSeq[float](e.N)
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
