import pbsolv_common
import simd/[sse2, vmath_sse2]
{.localPassc: simdCFlag.}

{.checks: not defined(release).}

def_ngs_iter_vec()
