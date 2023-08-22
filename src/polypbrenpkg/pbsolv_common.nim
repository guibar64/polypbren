import math

const NionsMax* = 8
type
  EqModelKind* = enum
    PlanePBE, RadialPBE
  SaltKind* = enum
    skGeneric, skIsoVal
  BoundaryCondition* = enum
    BCval,BCder
  PBEquationPars* = object
    kind*: EqModelKind
    saltKind*: SaltKind
    bcL*,bcR*: BoundaryCondition
    N*, nions*: int
    s0*, f0*: float
    phiL*, phiLd*, phiR*,phiRd*: float
    lb*, rho*: float
    fac*,z*: array[Nionsmax, float]

{.checks: not defined(release).}


proc chargeAndDerivTerms(e: PBEquationPars,phi: float): tuple[ch: float,der: float] {.inline.} =
  result.ch = e.rho
  for i in 0..<e.nions:
    let
      z = e.z[i]
      fl = e.fac[i]*exp(-z*phi)
    result.ch += fl
    result.der -= z*fl

proc chargeAndDerivTermsIsoVal(e: PBEquationPars,phi: float): tuple[ch: float,der: float] {.inline.} =
  let z=e.z[0]
  let dump=e.fac[0]
  let dump2=exp(-z*phi)
  let dump3=1.0/dump2
  result.ch = e.rho + dump*(dump2-dump3)
  result.der = -z*dump*(dump2+dump3)

func ngs_iter_novec*(RadialPlac: static[EqModelkind], saltKind: static[SaltKind], eq: var PBEquationPars,
  phi: var openArray[float], h: float): float {.inline.}  =
  # radial PBNewton Gauss-Seidel iteration
  let
    hh=h
    hh2=hh*hh
    hi2=1.0/(h*h)
    moins2hi2= -2.0*hi2
    ss0=eq.s0
  result = 0.0
  var sh=hh*ss0
  for i in 1..<phi.len-1:
    sh += hh2
    let (cht, numt) =
      when saltKind == skIsoVal: chargeAndDerivTermsIsoVal(eq, phi[i])
      else: chargeAndDerivTerms(eq, phi[i])
    let dlduo = moins2hi2 + numt
    var dent: float
    when RadialPlac == RadialPBE:
      let shi=1.0/sh
      dent = (phi[i+1] - 2.0*phi[i] + phi[i-1])*hi2 +
           (phi[i+1] - phi[i-1])*shi + cht
    else:
      dent = (phi[i+1] - 2.0*phi[i] + phi[i-1])*hi2 + cht
    phi[i] -= dent/dlduo
    result += dent*dent
  result /= float(phi.len-2)
  # BCs
  let nlast=phi.len-1
  if eq.bcR == BCder: phi[nlast] = phi[nlast-1] #+ eq.phiRd*h
  if eq.bcL == BCder: phi[0] = phi[1] - eq.phiLd*h


template def_ngs_iter_vec*()  {.dirty.} =
  func chargeAndDerivTermsVec(e: PBEquationPars, phi: VecF64): tuple[ch: VecF64,der: VecF64] {.inline.} =
    result.ch = vec(e.rho)
    result.der = vec(0.0)
    for i in 0..<e.nions:
      let
        z = vec(e.z[i])
        fl = vec(e.fac[i])*exp(vec(-1.0)*(z*phi))
      result.ch = result.ch + fl
      result.der = result.der - z*fl

  func chargeAndDerivTermsIsoValVec(e: PBEquationPars, phi: VecF64): tuple[ch: VecF64,der: VecF64] {.inline.} =
    let z=vec(e.z[0])
    let dump=vec(e.fac[0])
    let dump2=exp(vec(-1.0)*z*phi)
    let dump3=vec(1.0)/dump2
    result.ch = vec(e.rho) + dump*(dump2-dump3)
    result.der = vec(-1.0)*z*dump*(dump2+dump3)

  func ngs_iter_vec*(RadialPlac: static[EqModelkind], saltKind: static[SaltKind], eq: var PBEquationPars,
    phi: var openArray[float], h: float): float =
    # radial PBNewton Gauss-Seidel iteration
    let
      hh=h
      hh2=vec(hh*hh)
      hi2=vec(1.0/(h*h))
      moins2hi2= vec(-2.0)*hi2
      ss0=eq.s0
    var vresult = vec(0.0)
    let sh0=vec(hh*ss0)
    const vecSize = simdSize div sizeof(float)
    for i in countup(1, phi.len-2, vecSize):
      var shAr: array[vecSize, float]
      for k in 0..<vecSize: shAr[k] = (i+k).float
      let sh = sh0 +  vec(shAr)*hh2
      let phiI = load(addr phi[i]) # ith element is aligned
      let (cht, numt) =
        when saltKind == skIsoVal: chargeAndDerivTermsIsoValVec(eq, phiI)
        else: chargeAndDerivTermsVec(eq, phiI)
      let dlduo = moins2hi2 + numt
      let phiIm1 = loadu(addr phi[i-1]) # (i-1)th element is not aligned
      let phiIp1 = loadu(addr phi[i+1])
      var dent: VecF64
      when RadialPlac == RadialPBE:
        let shi=vec(1.0)/sh
        dent = (phiIp1 - vec(2.0)*phiI + phiIm1)*hi2 +
            (phiIp1 - phiIm1)*shi + cht
      else:
        dent = (phiIp1 - vec(2.0)*phiI + phiIm1)*hi2 + cht
      store addr phi[i], phiI - dent/dlduo # ith element is aligned
      vresult = vresult + dent*dent
    let vresult2 = toArray(vresult)
    result = 0.0
    for i in 0..<vresult2.len:
      result += vresult2[i] 
    result /= float(phi.len-2)
    # BCs
    let nlast=phi.len-1
    if eq.bcR == BCder: phi[nlast] = phi[nlast-1] #+ eq.phiRd*h
    if eq.bcL == BCder: phi[0] = phi[1] - eq.phiLd*h
