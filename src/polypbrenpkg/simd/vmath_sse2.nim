import ./sse2

## See `sse41` module for sqrt, ceil floor, round intrinsics

{.localPassc: simdCFlag.}

# fallback
import ./vmath_common
import math
import macros

macro importVF32Fallback(prc: untyped) =
  result = prc
  result.body = fallbackBody("M128", "float32", "m128", 4, prc)

macro importVF64Fallback(prc: untyped) =
  result = prc
  result.body = fallbackBody("M128d", "float64", "m128d", 2, prc)


when libmvec:
  {.pragma: importVF32, importc: "_ZGVbN4v_$1f".}
  {.pragma: importVF64, importc: "_ZGVbN2v_$1".}
elif defined(icc):
  {.pragma: importVF32, importc: "_mm_$1_ps", header: "<xmmintrin.h>"}
  {.pragma: importVF64, importc: "_mm_$1_pd", header: "<xmmintrin.h>"}
else:
  template importVF32(prc: untyped) = importVF32Fallback(prc)
  template importVF64(prc: untyped) = importVF64Fallback(prc)

func exp*(x: M128): M128 {.importVF32.}
func pow*(x, y: M128): M128 {.importVF32.}
func log*(x: M128): M128 {.importVF32.}
func sin*(x: M128): M128 {.importVF32.}
func cos*(x: M128): M128 {.importVF32.}

# The following does not seem to be in libmvec
when defined(icc):
  func tan*(x: M128): M128 {.importVF32.}
  func acos*(x: M128): M128 {.importVF32.}
  func asin*(x: M128): M128 {.importVF32.}
  func atan*(x: M128): M128 {.importVF32.}
else:
  func tan*(x: M128): M128 {.inline.} =
    sin(x)/cos(x)
  func acos*(x: M128): M128 {.importVF32Fallback.}
  func asin*(x: M128): M128 {.importVF32Fallback.}
  func atan*(x: M128): M128 {.importVF32Fallback.}

func exp*(x: M128d): M128d {.importVF64.}
func pow*(x, y: M128d): M128d {.importVF64.}
func log*(x: M128d): M128d {.importVF64.}
func sin*(x: M128d): M128d {.importVF64.}
func cos*(x: M128d): M128d {.importVF64.}

when defined(icc):
  func tan*(x: M128d): M128d {.importVF64.}
  func acos*(x: M128d): M128d {.importVF64.}
  func asin*(x: M128d): M128d {.importVF64.}
  func atan*(x: M128d): M128d {.importVF64.}
else:
  func tan*(x: M128d): M128d {.inline.} =
    sin(x)/cos(x)
  func acos*(x: M128d): M128d {.importVF64Fallback.}
  func asin*(x: M128d): M128d {.importVF64Fallback.}
  func atan*(x: M128d): M128d {.importVF64Fallback.}

