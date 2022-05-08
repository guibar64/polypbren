import std/os

let simdEnv = getEnv("PBSOLV_SIMD", "NONE")
let 
  hasSSE2* = simdEnv == "SSE2"
  hasAVX2* = simdEnv == "AVX2"
  hasAVX512f* = simdEnv == "AVX512"



