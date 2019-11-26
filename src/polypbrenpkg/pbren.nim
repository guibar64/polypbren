#********************** Licensing information ****************************
# Copyright (c) 2018 Guillaume Bareigts
# MIT License, see LICENSE
#*************************************************************************
## This module computes renormalized charges of polydisperse charges hard
## spheres and the dispersion effective debye screening length.
## The polydisperse cell model, jellium and the virtual cell matching
## are available.

import math
import threadpool_simple
import pbsolv, params, distrib

proc calcBackDensity*(phiB:float, zv, concv: openArray[float]): float =
  assert(concv.len >= zv.len)
  result = 0.0
  for i,z in zv:
    result -= z*concv[i]*exp(-z*phiB)

proc calcPhiBfromBackDensity*(rhoB: float, zv, concv: openArray[float], p0=0.0, eps=1.0e-6): float =
  assert(concv.len >= zv.len)
  let epsa=abs(eps*rhoB)
  var p = p0
  for it in 1..20:
    var
      ch=rhoB
      der = 0.0
    for i,z in zv:
      let fl = z*concv[i]*exp(-z*p)
      ch += fl
      der += z*fl
    p += ch/der
    if abs(ch)<epsa: break
  result = p

    
proc calcIonDensity*(phiB,lb: float, zv, concv: openArray[float]): float =
  assert(concv.len >= zv.len)
  result = 0.0
  for i,z in zv:
    result += concv[i]*exp(-z*phiB)

proc calcKappaEff*(phiB,lb: float, zv, concv: openArray[float]): float =
  assert(concv.len >= zv.len)
  var res = 0.0
  for i,z in zv:
    res += z*z*concv[i]*exp(-z*phiB)
  result = sqrt(res*4*PI*lb)

proc calcZeffJellium*(rv, pv: openArray[float], pD , rp, keff: float): float =
  #let minR = FactCalcZeffJellium/keff
  var
    sNum=0.0
    sDem=0.0
  for i,r in rv:
    let dp = pv[i] - pD
    if abs(dp) <= max_potential_zeff_jellium:
      sNum += dp
      sDem += exp(-keff*(r-rp))/r
  result = (1.0+keff*rp)/(4.0*Pi*rp*rp*lambdaB) * sNum/sDem
  

proc calcZeffCellModel*(phiR, lb, R, a: float, zv, concv: openArray[float] ): float =
  let
    lb = lambdaB
    keff = calcKappaEff(phiR, lb, zv, concv)
    keff2 = keff*keff
    gamma0 = calcBackDensity(phiR, zv, concv)*(4*Pi*lb)/keff2
    dump=keff*(R-a)
    dump2 = gamma0*((keff2*a*R - 1.0)*sinh(dump) + dump*cosh(dump))
  result = dump2/(4*PI*a*a*lb*keff)

proc solvation*(eq: var PBEquation, fam: FamComp, finalSig: var float): int =
  ## Solves the PB equation (represented in ``eq``), for a family ``fam``.
  ## The final charge is returned in ``finalSig`` if there is titration.
  ## Returns the number of steps effected.
  
  let rp=fam.r
  let chin=fam.ch
  
  proc deltaPH(e: PBEquation): float =
    let sig = sigma(e)
    let signSig = sig/abs(sig)
    let psiS = e.phi[0] + 4.0*Pi*lambdaB*sternL/(1.0+sternL/rp)*sig
    let alpha = min(abs(sig)/(siteDens)*pow((rp+sternL)/rp,2), 0.999999)
    let pH = pKa + (ln(alpha/(1.0-alpha)) + signSig*psiS)/ln(10.0)
    result = pH - pHInput

  eq.defBoundCond(bcL = BCDer, bcR= BCDer, left= -4*Pi*lambdaB*chin, right=0.0)
  result = eq.solve(tol = tolInLoop, niter = maxInIters)
  let
    sig1=sigma(eq)
    phi01=eq.phi[0]
  eq.defBoundCond(bcL = BCVal, bcR= BCDer, left= phi01, right=0.0)
  case solveCondition
  of PHCond:
    result = eq.solveByCondition(deltaPH, tolF = tolPhi0*abs(pHInput), tolIn = tolInLoop, nitermax = maxIters, niterin = maxInIters, maxVal = 2.0*eq.phi[0])
    finalSig = eq.sigma()
  of ChargeCond:
    let fac = factPhi0ChargeCond*(chin-sig1)/chin*phi01
    result = eq.solveByNeutralization(tolF = tolPhi0, tolIn = tolInLoop, nitermax = maxIters, niterin = maxInIters, minVal = phi01-fac, maxVal=phi01+fac)


type PbrenRes* = object
  volfrac*: float      ## final volume fraction
  kappaEff*: float     ## effective inverse screening length
  rhoEdge*: float       ## ion concentration at the edge
  sigmaEffv*: seq[float] ## effective charge densities 
  finalSig*: seq[float]  ## bare charge densities 
  finalPhiv*: seq[seq[float]] ## holds electrical potential profiles foreach family
  finalMesh*: seq[seq[float]] ## holds spatial grid foreach family

proc doCalculations*(distribution: Distrib, maxThreads = 1): PbrenRes =
  ## Calculates potential profiles, renormalized parameters for a given ``distribution``.
  let tp = newThreadPool(maxThreads)
  var (volmoy, nptot) = calcVolmNptot(distribution)

  let dlen = 1.0/calcKappaEff(0.0, lambdaB, ionChv, ionDensv)
  
  result.sigmaEffv = newSeq[float](distribution.len)
  result.finalPhiv = newSeq[seq[float]](distribution.len)
  result.finalMesh = newSeq[seq[float]](distribution.len)
  result.finalSig = newSeq[float](distribution.len)
  case model
  of Jellium,Jellbid:
    if verb >= 1:
      case model
      of Jellium:
        echo "Running jellium model..."    
      of Jellbid:
        echo "Running virtual-cell-matching jellium model..."
      else: discard
    var volfrac: float
    for step in 1..maxSteps:
      result.kappaEff = calcKappaEff(phiD, lambdaB, ionChv, ionDensv)
      var rhoBack = calcBackDensity(phiD, ionChv, ionDensv)

      type TaskRetType = tuple[virtVol, zmoytemp, finalSig, sigmaEffv: float, finalMesh, finalPhiv: seq[float]]   
      proc task(fam: FamComp, ionDensv, ionChv: seq[float], dlen, kappaEff, rhoBack: float): TaskRetType = 
        let
          s0 = fam.r
          sR = s0+slength
          Zr = 4.0*PI*s0*s0*fam.ch
        proc initer(x: float): float =
          result = lambdaB*Zr/(1.0 + s0/dlen)*exp(-(x-s0)/dlen)/x
        var eq = initPBEquation(s0, sR, sbin, zv = ionChv, concv=ionDensv, kind = eqKind,
                                lambdab = lambdaB, charge = fam.ch,rho = rhoBack, gridStep = gridStep)
        eq.initialize(initer)
        let ires = eq.residual
        var itereff = solvation(eq, fam, result.finalSig)
        let
          chin = 4*PI*s0*s0*sigma(eq)
          sigren = calcZeffJellium(eq.mesh(), eq.phi, phiD , fam.r, kappaEff)
          zren = sigren*4*PI*s0*s0
        if verb >= 2:
          echo "Nb of iterations: ", itereff, " R = ", fam.r," Z = ",Zr, " ChIn = ", chin
          echo "Sigma = ",fam.ch," SigmaEff = ",sigren, " Zeff = ", zren
          echo "Residue : Initial = ", ires, " Final = ", eq.residual
        result.zmoytemp += fam.np.float*zren
        result.sigmaEffv = sigren
        result.finalPhiv = eq.phi
        result.finalMesh = eq.mesh()

        case model
        of Jellium: discard
        of Jellbid:
          var cumCharge = eq.sigmaCumSansBack()
          var rtruc = result.finalMesh[eq.phi.high-1]
          var itruc=0
          var prev= cumCharge[0]
          for i in 1..eq.phi.high:
            var val = cumCharge[i]
            if val*prev <= 0.0:
              itruc=i
              break
            prev = val # is it really necessary ?
          if itruc > 0:
            let p = cumCharge[itruc]
            let pm1 = cumCharge[itruc-1]
            let h = result.finalMesh[itruc] - result.finalMesh[itruc-1]
            rtruc = s0 + h*(itruc.float-p/(p-pm1))
          result.virtVol = 4.0/3.0*Pi*fam.np.float*pow(rtruc,3)
        else: discard

      var res = newSeq[FlowVar[TaskRetType]](distribution.len)  
      for f, fam in distribution:
        res[f] = tp.spawn task(fam, ionDensv, ionChv, dlen, result.kappaEff, rhoBack)
      tp.sync()

      var virtVol = 0.0
      var zmoytemp = 0.0
      for f in 0..<distribution.len:
        let resf = ^res[f]
        virtVol += resf.virtVol
        zmoytemp += resf.zmoytemp
        result.sigmaEffv[f] = resf.sigmaEffv
        result.finalPhiv[f] = resf.finalPhiv
        result.finalMesh[f] = resf.finalMesh
        result.finalSig[f] = resf.finalSig

      var phiDnew=0.0
      case model
      of Jellium:
        zmoytemp /= nptot.float
        volfrac = rhoBack*volmoy/zmoytemp
        let rhoBackNew=zmoytemp/volmoy*volfracIn
        phiDNew = calcPhiBfromBackDensity(rhoBackNew, ionChv, ionDensv)
      of Jellbid:
        virtVol /= nptot.float
        volfrac = volmoy/virtVol
        phiDNew = phiD
      else: discard

      if verb >= 2:
        echo "Step ",step,". Volume Fraction = ", volfrac, ", New PhiD= ", phiDNew
      if(abs(volfrac-volfracIn)<tolVolFrac*volfracIn): break
      phiD=phiDNew
      rhoBack = calcBackDensity(phiD, ionChv, ionDensv)
  
    result.volfrac = volfrac
    result.rhoEdge = calcIonDensity(phiD,lambdaB, ionChv, ionDensv)

  of Cellmod:
    if verb >= 1: echo "Running cell model..."
    result.kappaEff = calcKappaEff(phiD, lambdaB, ionChv, ionDensv)
    type TaskRetType = tuple[npR3Cell, finalSig, sigmaEffv: float, finalMesh, finalPhiv: seq[float]]   
    proc task(fam: FamComp, kappaEff: float, ionChv, ionDensv: seq[float]): TaskRetType {.gcSafe.} = 
      var chin,sigren, zren, sigren_ana: float
      let rp = fam.r
      let Zr = 4.0*PI*rp*rp*fam.ch
      var compz,compzold = if phiD<0.0: -1.0 else: 1.0
      var slimmin=rp
      var slimmax=slimmin+2.0*slength
      var slim = 0.5*(slimmax+slimmin)
      var rcell = slim
      var effStep = maxSteps
      for step in 1..maxSteps:
        rcell = slim
        proc initerCell(x: float): float =        
          let
            fplus=(kappaEff*rcell+1)/(2.0*kappaEff)*exp(-kappaEff*rcell)
            fmoins=(kappaEff*rcell-1)/(2.0*kappaEff)*exp(+kappaEff*rcell)        
            pref=lambdaB*Zr/(fplus*(1-kappaEff*rp)*exp(kappaEff*rp)+fmoins*(1+kappaEff*rp)*exp(-kappaEff*rp))
          result = pref*(-1.0+(fplus*exp(+kappaEff*x)+fmoins*exp(-kappaEff*x))/x)

        var eq = initPBEquation(rp,slim, sbin, zv = ionChv, concv=ionDensv, kind = eqKind,
                              lambdab = lambdaB, charge = fam.ch,rho = 0.0, gridStep = gridStep)
        eq.initialize(initerCell)
        var itereff = solvation(eq, fam, result.finalSig)

        chin = 4*PI*rp*rp*sigma(eq)
        sigren = calcZeffCellModel(eq.phi[eq.phi.high],lambdaB,rcell,fam.r, ionChv, ionDensv)
        zren = sigren*4*PI*rp*rp
        if verb >= 2:
          echo "Nb of iterations: ", itereff, " R = ", fam.r," Z = ",Zr, " ChIn = ", chin, " Residue = ", eq.residual
          echo "Sigma = ",fam.ch," SigmaEff = ",sigren, " Zeff = ", zren
          echo " Rcell = ",rcell, " Phi(Rcell)= ", eq.phi[eq.phi.high]
        result.finalPhiv = eq.phi
        result.finalMesh = eq.mesh()
        compz =  eq.phi[eq.phi.high] - phiD
        if abs(compz) <= abs(tolCellm*phiD): effStep=step ; break
        if maxSteps == 1: # in case of a standard monodisperse cell model
          phiD = eq.phi[eq.phi.high]
          break 
        #if abs(slimmax-slimmin) <= Sbin: break     # In theory: you gotta change Sbin
        if compz*compzold > 0.0: slimmin = slim  #; compzold = compz
        else: slimmax = slim 
        slim = 0.5*(slimmax+slimmin)

      sigren_ana = calcZeffCellModel(phiD,lambdaB,rcell,fam.r, ionChv, ionDensv)
      if verb >= 1:
        echo "Steps = ", effStep, " Rcell = ", rcell, " PhiR = ", result.finalPhiv[^1]," (",phiD,")", " SigmaEff2= ", sigren_ana
      result.sigmaEffv = sigren_ana
      # sigmaEffvAlt[f] = sigren
      result.npR3Cell = fam.np.float*pow(slim,3)

    var res = newSeq[FlowVar[TaskRetType]](distribution.len)
    for f,fam in distribution:
      res[f] = tp.spawn task(fam, result.kappaEff, ionChv, ionDensv)
    tp.sync()

    var vmoycell = 0.0
    for f in 0..<distribution.len:
      let resf = ^res[f]
      vmoycell += resf.npR3Cell
      result.finalPhiv[f] = resf.finalPhiv
      result.finalMesh[f] = resf.finalMesh
      result.finalSig[f] = resf.finalSig
      result.sigmaEffv[f] = resf.sigmaEffv
    vmoycell *= (4*PI)/(3.0*nptot.float)
    result.volfrac = volmoy/vmoycell
    if verb >= 1:
      echo "Volume Fraction = ", result.volfrac, " Effective Debye = ", 1.0/result.kappaEff
    result.rhoEdge = calcIonDensity(phiD,lambdaB, ionChv, ionDensv)
  else:
    quit "Invalid model. Please check your input.", 6

