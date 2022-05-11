
type
  M512* {.importc:"__m512", header: "<immintrin.h>".} = object

const
  simdCFlag* = when defined(gcc) or defined(clang) or defined(icc): "-mavx512f" else: ""

{.localPassc:simdCFlag.}


func m512*(): M512 {.importc: "_mm512_setzero_ps", header: "<immintrin.h>".}
func m512*(a: float32): M512 {.importc: "_mm512_set1_ps", header: "<immintrin.h>".}
func m512*(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15: float32): M512 {.importc: "_mm512_setr_ps", header: "<immintrin.h>".}

func load*(a: ptr float32): M512 {.importc: "_mm512_load_ps", header: "<immintrin.h>".}
func loadu*(a: ptr float32): M512 {.importc: "_mm512_loadu_ps", header: "<immintrin.h>".}
func store*(a: ptr float32, b: M512) {.importc: "_mm512_store_ps", header: "<immintrin.h>".}
func storeu*(a: ptr float32, b: M512) {.importc: "_mm512_storeu_ps", header: "<immintrin.h>".}

func store1*(a: var float32, b: M512) {.importc: "_mm512_store1_ps", header: "<immintrin.h>".}



func m512*(a: array[16, float32]): M512 {.inline, noinit.} =
  loadu(unsafeaddr a[0])

func toArray*(a: M512): array[16, float32] {.inline, noinit.} =
  storeu(addr result[0], a)

func `$`*(a: M512): string =
  let temp = a
  result = $toArray(temp)

func `+`*(a, b: M512): M512 {.importc: "_mm512_add_ps", header: "<immintrin.h>".}
func `-`*(a, b: M512): M512 {.importc: "_mm512_sub_ps", header: "<immintrin.h>".}
func `*`*(a, b: M512): M512 {.importc: "_mm512_mul_ps", header: "<immintrin.h>".}
func `/`*(a, b: M512): M512 {.importc: "_mm512_div_ps", header: "<immintrin.h>".}

func min*(a: M512): M512 {.importc: "_mm512_min_ps", header: "<immintrin.h>".}
func max*(a: M512): M512 {.importc: "_mm512_max_ps", header: "<immintrin.h>".}
func sqrt*(a: M512): M512 {.importc: "_mm512_sqrt_ps", header: "<immintrin.h>".}

func ceil*(a: M512): M512 {.importc: "_mm512__ceil_ps", header: "<immintrin.h>".}
func floor*(a: M512): M512 {.importc: "_mm512__floor_ps", header: "<immintrin.h>".}
func round*(a: M512, b: int32): M512 {.importc: "_mm512__round_ps", header: "<immintrin.h>".}
  ## ``b`` must be a 8-bit constant



type
  M512d* {.importc:"__m512d", header: "<immintrin.h>".} = object

func m512d*(): M512d {.importc: "_mm512_setzero_pd", header: "<immintrin.h>".}
func m512d*(a: float64): M512d {.importc: "_mm512_set1_pd", header: "<immintrin.h>".}
func m512d*(a0,a1,a2,a3,a4,a5,a6,a7: float64): M512d {.importc: "_mm512_setr_pd", header: "<immintrin.h>".}

func load*(a: ptr float64): M512d {.importc: "_mm512_load_pd", header: "<immintrin.h>".}
func loadu*(a: ptr float64): M512d {.importc: "_mm512_loadu_pd", header: "<immintrin.h>".}
func store*(a: ptr float64, b: M512d) {.importc: "_mm512_store_pd", header: "<immintrin.h>".}
func storeu*(a: ptr float64, b: M512d) {.importc: "_mm512_storeu_pd", header: "<immintrin.h>".}

# func store1*(a: var float64, b: M512d) {.importc: "_mm512_store1_pd", header: "<immintrin.h>".}

func m512d*(a: array[8, float64]): M512d {.inline.} =
  loadu(unsafeaddr a[0])

func toArray*(a: M512d): array[8, float64] {.inline.} =
  storeu(addr result[0], a)


func `$`*(a: M512d): string =
  let temp = a
  result = $toArray(temp)


func `+`*(a, b: M512d): M512d {.importc: "_mm512_add_pd", header: "<immintrin.h>".}
func `-`*(a, b: M512d): M512d {.importc: "_mm512_sub_pd", header: "<immintrin.h>".}
func `*`*(a, b: M512d): M512d {.importc: "_mm512_mul_pd", header: "<immintrin.h>".}
func `/`*(a, b: M512d): M512d {.importc: "_mm512_div_pd", header: "<immintrin.h>".}

func min*(a,b: M512d): M512d {.importc: "_mm512_min_pd", header: "<immintrin.h>".}
func max*(a,b: M512d): M512d {.importc: "_mm512_max_pd", header: "<immintrin.h>".}
func sqrt*(a: M512d): M512d {.importc: "_mm512_sqrt_pd", header: "<immintrin.h>".}

func ceil*(a: M512d): M512d {.importc: "_mm512__ceil_pd", header: "<immintrin.h>".}
func floor*(a: M512d): M512d {.importc: "_mm512__floor_pd", header: "<immintrin.h>".}
func round*(a: M512d, b: int32): M512d {.importc: "_mm512__round_pd", header: "<immintrin.h>".}

const
  simdSize* = 64
type
  VecF32* = M512
  VecF64* = M512d
template vec*(a: array[simdSize div 4, float32]): VecF32 = m512(a)
template vec*(a: array[simdSize div 8, float64]): VecF64 = m512d(a)
template vec*(a: float32): VecF32 = m512(a)
template vec*(a: float64): VecF64 = m512d(a)

when isMainModule:
  block:
    let v1 = m512([1.0f32,2,3,4,5,6,7,8,1.0f32,2,3,4,5,6,7,8])
    let v2 = m512([-1f32,-2, -6, -8, 1,1,1,1,-1f32,-2, -6, -8, 1,1,1,1])
    let v3 = v1 + v2
    echo toArray(v3)
    # why echo sqrt(v1) crash on win ?
    echo sqrt(v1)
    echo v2.toArray[0]
  block:
    let v1 = m512d([1.0f64,2,3,4,1.0f64,2,3,4])
    let v2 = m512d([-1f64,-4, -6, 1,-1f64,-4, -6, 1])
    let v3 = v1 + v2
    echo toArray(v3)
    let v4 = sqrt(v1)
    echo v4
    echo m512d(1.0f64).toArray[3]
