import polypbrenpkg/pbsolv

import math

proc appr(x,y,eps: float): bool =
  abs(x-y) < eps*abs(y)

let zv = [1.0, -1.0]
let conc = [0.5, 0.5]
var eq = initPBEquation(smin = 0.0, smax = 10.0, sbin = 0.02, zv = zv, concv=conc, kind=PlanePBE ,left = 1.0, lambdab = 1.0)
eq.initialize(1.0)
echo eq.solve(niter = 10000, tol = 0.00001)
doAssert appr(eq.sigma, 0.293997, 1.0e-3)

var eq2 =  initPBEquation(smin = 0.0, smax = 10.0, sbin = 0.01, zv = zv, concv=conc, kind=PlanePBE, lambdab = 1.0,
  bcL = BCDer, bcR= BCDer, left= -4*Pi*1.0*0.3, right=0.0)
eq2.initialize(1.0)
echo eq2.solve(niter = 20000, tol = 0.00001)
doAssert appr(eq2.sigma, 0.3, 5.0e-2)
doAssert appr(eq2.phi[0],1.0188318, 2.0e-2)
