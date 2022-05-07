
const
  simdCFlag* = when defined(gcc) or defined(clang) or defined(icc):  "-msse2" else: ""

{.localPassc: simdCFlag.}

type
  M128* {.importc:"__m128", header: "<xmmintrin.h>".} = object

func m128*(): M128 {.importc: "_mm_setzeroes_ps", header: "<xmmintrin.h>".}
func m128*(a: float32): M128 {.importc: "_mm_set1_ps", header: "<xmmintrin.h>".}
func m128*(a,b,c,d: float32): M128 {.importc: "_mm_setr_ps", header: "<xmmintrin.h>".}

func load*(a: ptr float32): M128 {.importc: "_mm_load_ps", header: "<xmmintrin.h>".}
func loadu*(a: ptr float32): M128 {.importc: "_mm_loadu_ps", header: "<xmmintrin.h>".}
func store*(a: ptr float32, b: M128) {.importc: "_mm_store_ps", header: "<xmmintrin.h>".}
func storeu*(a: ptr float32, b: M128) {.importc: "_mm_storeu_ps", header: "<xmmintrin.h>".}

func store1*(a: var float32, b: M128) {.importc: "_mm_store1_ps", header: "<xmmintrin.h>".}

func m128*(a: array[4, float32]): M128 {.inline.} =
  loadu(unsafeaddr a[0])

func toArray*(a: M128): array[4, float32] {.inline.} =
  storeu(addr result[0], a)

func `$`*(a: M128): string =
  $toArray(a)

func `+`*(a, b: M128): M128 {.importc: "_mm_add_ps", header: "<xmmintrin.h>".}
func `-`*(a, b: M128): M128 {.importc: "_mm_sub_ps", header: "<xmmintrin.h>".}
func `*`*(a, b: M128): M128 {.importc: "_mm_mul_ps", header: "<xmmintrin.h>".}
func `/`*(a, b: M128): M128 {.importc: "_mm_div_ps", header: "<xmmintrin.h>".}

func min*(a: M128): M128 {.importc: "_mm_min_ps", header: "<xmmintrin.h>".}
func max*(a: M128): M128 {.importc: "_mm_max_ps", header: "<xmmintrin.h>".}
func sqrt*(a: M128): M128 {.importc: "_mm_sqrt_ps", header: "<xmmintrin.h>".}


type
  M128d* {.importc:"__m128d", header: "<emmintrin.h>".} = object

func m128d*(): M128d {.importc: "_mm_setzero_pd", header: "<emmintrin.h>".}
func m128d*(a: float64): M128d {.importc: "_mm_set1_pd", header: "<emmintrin.h>".}
func m128d*(a,b: float64): M128d {.importc: "_mm_setr_pd", header: "<emmintrin.h>".}

func load*(a: ptr float64): M128d {.importc: "_mm_load_pd", header: "<emmintrin.h>".}
func loadu*(a: ptr float64): M128d {.importc: "_mm_loadu_pd", header: "<emmintrin.h>".}
func store*(a: ptr float64, b: M128d) {.importc: "_mm_store_pd", header: "<emmintrin.h>".}
func storeu*(a: ptr float64, b: M128d) {.importc: "_mm_storeu_pd", header: "<emmintrin.h>".}

# func store1*(a: var float64, b: M128d) {.importc: "_mm_store1_pd", header: "<emmintrin.h>".}

func m128d*(a: array[2, float64]): M128d {.inline.} =
  loadu(unsafeaddr a[0])

func toArray*(a: M128d): array[2, float64] {.inline.} =
  storeu(addr result[0], a)

func `$`*(a: M128d): string =
  $toArray(a)

func `+`*(a, b: M128d): M128d {.importc: "_mm_add_pd", header: "<emmintrin.h>".}
func `-`*(a, b: M128d): M128d {.importc: "_mm_sub_pd", header: "<emmintrin.h>".}
func `*`*(a, b: M128d): M128d {.importc: "_mm_mul_pd", header: "<emmintrin.h>".}
func `/`*(a, b: M128d): M128d {.importc: "_mm_div_pd", header: "<emmintrin.h>".}

func min*(a,b: M128d): M128d {.importc: "_mm_min_pd", header: "<emmintrin.h>".}
func max*(a,b: M128d): M128d {.importc: "_mm_max_pd", header: "<emmintrin.h>".}
func sqrt*(a: M128d): M128d {.importc: "_mm_sqrt_pd", header: "<emmintrin.h>".}

  
type
  M128i* {.importc:"__m128i", header: "<emmintrin.h>".} = object

const
  simdSize* = 16
type
  VecF32* = M128
  VecF64* = M128d
template vec*(a: array[simdSize div 4, float32]): VecF32 = m128(a)
template vec*(a: array[simdSize div 8, float64]): VecF64 = m128d(a)
template vec*(a: float32): VecF32 = m128(a)
template vec*(a: float64): VecF64 = m128d(a)

when isMainModule:
  block:
    let v1 = m128([1.0f32,2,3,4])
    let v2 = m128(-1f32,-2, -6, -8)
    let v3 = v1 + v2
    echo toArray(v3)
    echo sqrt(v1)
    echo v2.toArray[0]
  block:
    let v1 = m128d([1.0f64,2])
    let v2 = m128d(-1f64,-4)
    let v3 = v1 + v2
    echo toArray(v3)
    echo sqrt(v1)
    echo m128d(1.0f64).toArray[1]
