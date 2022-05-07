
const
  simdCFlag* = when defined(gcc) or defined(clang) or defined(icc): "-mavx2" else: ""

{.localPassc: simdCFlag.}


type
  M256* {.importc:"__m256", header: "<immintrin.h>".} = object


func m256*(): M256 {.importc: "_mm256_setzero_ps", header: "<immintrin.h>".}
func m256*(a: float32): M256 {.importc: "_mm256_set1_ps", header: "<immintrin.h>".}
func m256*(a0,a1,a2,a3,a4,a5,a6,a7: float32): M256 {.importc: "_mm256_setr_ps", header: "<immintrin.h>".}

func load*(a: ptr float32): M256 {.importc: "_mm256_load_ps", header: "<immintrin.h>".}
func loadu*(a: ptr float32): M256 {.importc: "_mm256_loadu_ps", header: "<immintrin.h>".}
func store*(a: ptr float32, b: M256) {.importc: "_mm256_store_ps", header: "<immintrin.h>".}
func storeu*(a: ptr float32, b: M256) {.importc: "_mm256_storeu_ps", header: "<immintrin.h>".}

func store1*(a: var float32, b: M256) {.importc: "_mm256_store1_ps", header: "<immintrin.h>".}



func m256*(a: array[8, float32]): M256 {.inline, noinit.} =
  loadu(unsafeaddr a[0])

func toArray*(a: M256): array[8, float32] {.inline, noinit.} =
  storeu(addr result[0], a)

func `$`*(a: M256): string =
  let temp = a
  result = $toArray(temp)

func `+`*(a, b: M256): M256 {.importc: "_mm256_add_ps", header: "<immintrin.h>".}
func `-`*(a, b: M256): M256 {.importc: "_mm256_sub_ps", header: "<immintrin.h>".}
func `*`*(a, b: M256): M256 {.importc: "_mm256_mul_ps", header: "<immintrin.h>".}
func `/`*(a, b: M256): M256 {.importc: "_mm256_div_ps", header: "<immintrin.h>".}

func min*(a: M256): M256 {.importc: "_mm256_min_ps", header: "<immintrin.h>".}
func max*(a: M256): M256 {.importc: "_mm256_max_ps", header: "<immintrin.h>".}
func sqrt*(a: M256): M256 {.importc: "_mm256_sqrt_ps", header: "<immintrin.h>".}

func ceil*(a: M256): M256 {.importc: "_mm256__ceil_ps", header: "<immintrin.h>".}
func floor*(a: M256): M256 {.importc: "_mm256__floor_ps", header: "<immintrin.h>".}
func round*(a: M256, b: int32): M256 {.importc: "_mm256__round_ps", header: "<immintrin.h>".}
  ## ``b`` must be a 8-bit constant



type
  M256d* {.importc:"__m256d", header: "<immintrin.h>".} = object

func m256d*(): M256d {.importc: "_mm256_setzero_pd", header: "<immintrin.h>".}
func m256d*(a: float64): M256d {.importc: "_mm256_set1_pd", header: "<immintrin.h>".}
func m256d*(a0,a1,a2,a3: float64): M256d {.importc: "_mm256_setr_pd", header: "<immintrin.h>".}

func load*(a: ptr float64): M256d {.importc: "_mm256_load_pd", header: "<immintrin.h>".}
func loadu*(a: ptr float64): M256d {.importc: "_mm256_loadu_pd", header: "<immintrin.h>".}
func store*(a: ptr float64, b: M256d) {.importc: "_mm256_store_pd", header: "<immintrin.h>".}
func storeu*(a: ptr float64, b: M256d) {.importc: "_mm256_storeu_pd", header: "<immintrin.h>".}

# func store1*(a: var float64, b: M256d) {.importc: "_mm256_store1_pd", header: "<immintrin.h>".}

func m256d*(a: array[4, float64]): M256d {.inline.} =
  loadu(unsafeaddr a[0])

func toArray*(a: M256d): array[4, float64] {.inline.} =
  storeu(addr result[0], a)


func `$`*(a: M256d): string =
  let temp = a
  result = $toArray(temp)


func `+`*(a, b: M256d): M256d {.importc: "_mm256_add_pd", header: "<immintrin.h>".}
func `-`*(a, b: M256d): M256d {.importc: "_mm256_sub_pd", header: "<immintrin.h>".}
func `*`*(a, b: M256d): M256d {.importc: "_mm256_mul_pd", header: "<immintrin.h>".}
func `/`*(a, b: M256d): M256d {.importc: "_mm256_div_pd", header: "<immintrin.h>".}

func min*(a,b: M256d): M256d {.importc: "_mm256_min_pd", header: "<immintrin.h>".}
func max*(a,b: M256d): M256d {.importc: "_mm256_max_pd", header: "<immintrin.h>".}
func sqrt*(a: M256d): M256d {.importc: "_mm256_sqrt_pd", header: "<immintrin.h>".}

func ceil*(a: M256d): M256d {.importc: "_mm256__ceil_pd", header: "<immintrin.h>".}
func floor*(a: M256d): M256d {.importc: "_mm256__floor_pd", header: "<immintrin.h>".}
func round*(a: M256d, b: int32): M256d {.importc: "_mm256__round_pd", header: "<immintrin.h>".}

const
  simdSize* = 32
type
  VecF32* = M256
  VecF64* = M256d
template vec*(a: array[simdSize div 4, float32]): VecF32 = m256(a)
template vec*(a: array[simdSize div 8, float64]): VecF64 = m256d(a)
template vec*(a: float32): VecF32 = m256(a)
template vec*(a: float64): VecF64 = m256d(a)

when isMainModule:
  block:
    let v1 = m256([1.0f32,2,3,4,5,6,7,8])
    let v2 = m256(-1f32,-2, -6, -8, 1,1,1,1)
    let v3 = v1 + v2
    echo toArray(v3)
    # why echo sqrt(v1) crash ?
    echo sqrt(v1)
    echo v2.toArray[0]
  block:
    let v1 = m256d([1.0f64,2,3,4])
    let v2 = m256d(-1f64,-4, -6, 1)
    let v3 = v1 + v2
    echo toArray(v3)
    let v4 = sqrt(v1)
    echo v4
    echo m256d(1.0f64).toArray[3]
