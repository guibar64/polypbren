import sse2

export sse2

const
  simdCFlag* = when defined(gcc) or defined(clang) or defined(icc): "-msse4.1" else: ""

{.localPassc: simdCFlag.}


func ceil*(a: M128): M128 {.importc: "_mm_ceil_ps", header: "<smmintrin.h>".}
func floor*(a: M128): M128 {.importc: "_mm_floor_ps", header: "<smmintrin.h>".}
func round*(a: M128, b: int32): M128 {.importc: "_mm_round_ps", header: "<smmintrin.h>".}
  ## ``b`` must be a 4-bit constant

func ceil*(a: M128d): M128d {.importc: "_mm_ceil_pd", header: "<smmintrin.h>".}
func floor*(a: M128d): M128d {.importc: "_mm_floor_pd", header: "<smmintrin.h>".}
func round*(a: M128d, b: int32): M128d {.importc: "_mm_round_pd", header: "<smmintrin.h>".}

when isMainModule:
  block:
    let v1 = m128 [1.2f32,2.6f32, -1.3f32, -4.8f32]
    echo ceil(v1)
    echo floor(v1)
    echo round(v1, 0b1100)
  block:
    let v1 = m128d [1.2f64,2.6]
    echo ceil(v1)
    echo floor(v1)
