import polypbrenpkg/pbsolv

import math

proc appr(x,y,eps: float): bool =
  abs(x-y) < eps*abs(y)

let zv = [1.0, -1.0]
let conc = [0.5, 0.5]
var eq = initPBEquation(smin = 0.0, smax = 10.0, sbin = 0.02, zv = zv, concv=conc, kind=PlanePBE ,left = 1.0, lambdab = 1.0)
eq.initialize(1.0)

eq.charge = 0.4
echo eq.solveByNeutralization(niterin = 10000, tolIn = 0.00001)
doAssert appr(eq.sigma, 0.4, 5.0e-4) 

eq.initialize(1.0)
echo eq.solveByCondition(proc(e: PBEquation):float = sigma(e) - 0.4, niterin = 10000, tolIn = 0.00001)
doAssert appr(eq.sigma, 0.4, 5.0e-4) 
