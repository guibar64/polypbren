import pbsolv_common
import simd/[avx512, vmath_avx512]
{.localPassc: simdCFlag.}

{.checks: not defined(release).}

def_ngs_iter_vec()
