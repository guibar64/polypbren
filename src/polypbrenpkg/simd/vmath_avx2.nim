from os import parentDir
import ./avx2


## See `avx2` module for sqrt, ceil floor, round intrinsics

{.localPassc: simdCFlag.}

# fallback
import ./vmath_common
import math
import macros

macro importVF32Fallback(prc: untyped) =
  result = prc
  result.body = fallbackBody("M256", "float32", "m256", 8, prc)

macro importVF64Fallback(prc: untyped) =
  result = prc
  result.body = fallbackBody("M256d", "float64", "m256d", 4, prc)


when libmvec:
  {.pragma: importVF32, importc: "_ZGVdN8v_$1f".}
  {.pragma: importVF64, importc: "_ZGVdN4v_$1".}
elif defined(icc):
  {.pragma: importVF32, importc: "_mm256_$1_ps", header: "<immintrin.h>"}
  {.pragma: importVF64, importc: "_mm256_$1_pd", header: "<immintrin.h>"}
else:
  template importVF32(prc: untyped) = importVF32Fallback(prc)
  template importVF64(prc: untyped) = importVF64Fallback(prc)

func exp*(x: M256): M256 {.importVF32.}
func pow*(x, y: M256): M256 {.importVF32.}
func log*(x: M256): M256 {.importVF32.}
func sin*(x: M256): M256 {.importVF32.}
func cos*(x: M256): M256 {.importVF32.}

# The following does not seem to be in libmvec
when defined(icc):
  func tan*(x: M256): M256 {.importVF32.}
  func acos*(x: M256): M256 {.importVF32.}
  func asin*(x: M256): M256 {.importVF32.}
  func atan*(x: M256): M256 {.importVF32.}
else:
  func tan*(x: M256): M256 {.inline.} =
    sin(x)/cos(x)
  func acos*(x: M256): M256 {.importVF32Fallback.}
  func asin*(x: M256): M256 {.importVF32Fallback.}
  func atan*(x: M256): M256 {.importVF32Fallback.}

func exp*(x: M256d): M256d {.importVF64.}
func pow*(x, y: M256d): M256d {.importVF64.}
func log*(x: M256d): M256d {.importVF64.}
func sin*(x: M256d): M256d {.importVF64.}
func cos*(x: M256d): M256d {.importVF64.}

when defined(icc):
  func tan*(x: M256d): M256d {.importVF64.}
  func acos*(x: M256d): M256d {.importVF64.}
  func asin*(x: M256d): M256d {.importVF64.}
  func atan*(x: M256d): M256d {.importVF64.}
else:
  func tan*(x: M256d): M256d {.inline.} =
    sin(x)/cos(x)
  func acos*(x: M256d): M256d {.importVF64Fallback.}
  func asin*(x: M256d): M256d {.importVF64Fallback.}
  func atan*(x: M256d): M256d {.importVF64Fallback.}

