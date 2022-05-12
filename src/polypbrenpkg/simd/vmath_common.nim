
# unclear where libvecm is available
const libmvec* {.booldefine.} = (not defined(windows) or defined(macos) or defined(noLibmvec)) and (defined(gcc) or defined(clang))

when libmvec:
  {.passl:"-lmvec"}


import macros

template fallbackBodyIpml*(vecType, floatType, vecLen, prc, fromArrayCv): untyped =
  var ret: array[vecLen, floatType]
  for i,f in toArray(x):
    ret[i] = prc(f)
  fromArrayCv(ret)

template fallbackBodyIpml2*(vecType, floatType, vecLen, prc, fromArrayCv): untyped =
  var ret: array[vecLen, floatType]
  let xa = toArray(x)
  let ya = toArray(y)
  for i in 0..<xa.len:
    ret[i] = prc(xa[i], ya[i])
  fromArrayCv(ret)


proc getProcBaseName*(prc: NimNode): NimNode =
  result = if prc[0].kind == nnkPostFix: prc[0][1] else: prc[0]
  case result.strVal
  of "log":
    result = ident"ln"
  of "acos":
    result = ident"arccos"
  of "acosh":
    result = ident"arccosh"
  of "asin":
    result = ident"arcsin"
  of "asinh":
    result = ident"arcsinh"
  of "atan":
    result = ident"arccos"
  of "atanh":
    result = ident"arctanh"
  else:
    discard


proc fallbackBody*(vecType, floatType, fromArrayCv: string, vecLen: int, prc: NimNode): NimNode =
  let prcName = getProcBaseName(prc)
  case prcName.strVal:
  of "pow":
    # Come on, write a proper macro
    result = getAst(fallbackBodyIpml2(ident vecType, ident floatType, newLit vecLen, prcName, ident fromArrayCv))
  else:
    result = getAst(fallbackBodyIpml(ident vecType, ident floatType, newLit vecLen, prcName, ident fromArrayCv))
