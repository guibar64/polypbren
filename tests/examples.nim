import std/[strutils, math, os, sequtils]
import polypbrenpkg / [pbren, params, distrib]

proc appr(x,y,eps: float): bool =
  abs(x-y) < eps*abs(y)

proc testExample*(dir: string, nth = 2) =
  loadParams(dir / "polypbren.cfg")

  let dist = readDistrib(dir / "distrib.in")

  let pbres = doCalculations(dist, nth)

  let peffExpected = readFile(dir / "expected_results" / "peffs.dat").splitWhitespace().mapIt(parseFloat(it))
  let (volmoy, _) = calcVolmNptot(dist)
  let press = kBpnm3 * temperature * (pbres.volfrac/(volmoy) + pbres.rhoEdge - ionDensv.sum)


  doAssert appr(pbres.volfrac, peffExpected[0], 2.0e-2)
  doAssert appr(pbres.kappaEff, peffExpected[1], 2.0e-2)
  doAssert appr(press, peffExpected[2], 2.0e-2)
