## Test cpu features
# see https://en.wikipedia.org/wiki/CPUID
import std/bitops
when defined(amd64) or defined(i386):

  {.push stackTrace:off.}
  proc cpuid(eax: int32, ecx: int32 = 0): tuple[eax,ebx,ecx,edx: int32] =
    # a in eax, and b in edx
    var reax, rebx, recx, redx: int32
    asm """
      cpuid
      : "=a"(`reax`),"=b"(`rebx`),"=c"(`recx`),"=d"(`redx`)
      : "a"(`eax`),"c"(`ecx`)
    """
    (reax, rebx, recx, redx)
  {.pop.}

  let
    featInfo = cpuid(eax = 1)
    featExt = cpuid(eax = 7, ecx = 0)

  let 
    hasSSE2* = featInfo.edx.testBit(26)
    hasAVX2* = featExt.ebx.testBit(5)
    hasAVX512f* = featExt.ebx.testBit(16)

else:
  let 
    hasSSE2* = false
    hasAVX2* = false
    hasAVX512f* = false


