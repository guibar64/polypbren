import ./avx512

## See `avx512` module for sqrt, ceil floor, round intrinsics

{.localPassc: simdCFlag.}

# fallback
import ./vmath_common
import math
import macros

macro importVF32Fallback(prc: untyped) =
  result = prc
  result.body = fallbackBody("M512", "float32", "m512", 16, prc)

macro importVF64Fallback(prc: untyped) =
  result = prc
  result.body = fallbackBody("M512d", "float64", "m512d", 8, prc)


when libmvec:
  {.pragma: importVF32, importc: "_ZGVeN16v_$1f".}
  {.pragma: importVF64, importc: "_ZGVeN8v_$1".}
elif defined(icc):
  {.pragma: importVF32, importc: "_mm512_$1_ps", header: "<emmintrin.h>"}
  {.pragma: importVF64, importc: "_mm512_$1_pd", header: "<emmintrin.h>"}
else:
  template importVF32(prc: untyped) = importVF32Fallback(prc)
  template importVF64(prc: untyped) = importVF64Fallback(prc)

func exp*(x: M512): M512 {.importVF32.}
func pow*(x, y: M512): M512 {.importVF32.}
func log*(x: M512): M512 {.importVF32.}
func sin*(x: M512): M512 {.importVF32.}
func cos*(x: M512): M512 {.importVF32.}

# The following does not seem to be in libmvec
when defined(icc):
  func tan*(x: M512): M512 {.importVF32.}
  func acos*(x: M512): M512 {.importVF32.}
  func asin*(x: M512): M512 {.importVF32.}
  func atan*(x: M512): M512 {.importVF32.}
else:
  func tan*(x: M512): M512 {.inline.} =
    sin(x)/cos(x)
  func acos*(x: M512): M512 {.importVF32Fallback.}
  func asin*(x: M512): M512 {.importVF32Fallback.}
  func atan*(x: M512): M512 {.importVF32Fallback.}

func exp*(x: M512d): M512d {.importVF64.}
func pow*(x, y: M512d): M512d {.importVF64.}
func log*(x: M512d): M512d {.importVF64.}
func sin*(x: M512d): M512d {.importVF64.}
func cos*(x: M512d): M512d {.importVF64.}

when defined(icc):
  func tan*(x: M512d): M512d {.importVF64.}
  func acos*(x: M512d): M512d {.importVF64.}
  func asin*(x: M512d): M512d {.importVF64.}
  func atan*(x: M512d): M512d {.importVF64.}
else:
  func tan*(x: M512d): M512d {.inline.} =
    sin(x)/cos(x)
  func acos*(x: M512d): M512d {.importVF64Fallback.}
  func asin*(x: M512d): M512d {.importVF64Fallback.}
  func atan*(x: M512d): M512d {.importVF64Fallback.}

