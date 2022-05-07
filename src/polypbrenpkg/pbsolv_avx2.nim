import pbsolv_common
import simd/[avx2, vmath_avx2]
{.localPassc: simdCFlag.}

{.checks: not defined(release).}

def_ngs_iter_vec()
