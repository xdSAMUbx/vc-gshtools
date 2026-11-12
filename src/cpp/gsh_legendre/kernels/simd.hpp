// simd.hpp: the lane type V used by the column kernels, one definition per instruction set.
// Each translation unit defines exactly one of GSH_ISA_SCALAR, GSH_ISA_AVX2 or GSH_ISA_AVX512 before including
// this file; the definitions then live in namespace gsh::isa_scalar, gsh::isa_avx2 or gsh::isa_avx512, so the same
// kernel template (column_kernel.inl) compiles three times without symbol clashes.
// Every lane primitive below operates on WIDTH doubles at once: set1, load, store, add/sub/mul/div, fused
// multiply-add (fmadd, fnmadd, fmsub), abs, comparisons producing a mask M, mask logic (mand, mor, any) and blend.
#pragma once

// GSH_INLINE forces inlining of the lane primitives, so that the accumulator arrays of the kernels stay in
// registers. MSVC has no always_inline attribute, hence __forceinline there.
#if defined(_MSC_VER) && !defined(__clang__)
#define GSH_INLINE __forceinline
#else
#define GSH_INLINE inline __attribute__((always_inline))
#endif

// AVX-512 (F, DQ, VL): V wraps __m512d, WIDTH = 8. Compare instructions return an opmask register (__mmask8),
// so blend uses mask_blend. Requires the CPU and the OS to support ZMM state (checked in shs.cpp).
#if defined(GSH_ISA_AVX512)
#include <immintrin.h>
#define GSH_ISA_NS isa_avx512
namespace gsh::GSH_ISA_NS
{
  struct V
  {
    __m512d v;
  };
  using M = __mmask8;
  inline constexpr int WIDTH = 8;
  GSH_INLINE V set1(double a) { return {_mm512_set1_pd(a)}; }
  GSH_INLINE V load(const double *p) { return {_mm512_loadu_pd(p)}; }
  GSH_INLINE void store(double *p, V a) { _mm512_storeu_pd(p, a.v); }
  GSH_INLINE V add(V a, V b) { return {_mm512_add_pd(a.v, b.v)}; }
  GSH_INLINE V sub(V a, V b) { return {_mm512_sub_pd(a.v, b.v)}; }
  GSH_INLINE V mul(V a, V b) { return {_mm512_mul_pd(a.v, b.v)}; }
  GSH_INLINE V div(V a, V b) { return {_mm512_div_pd(a.v, b.v)}; }
  GSH_INLINE V fmadd(V a, V b, V c) { return {_mm512_fmadd_pd(a.v, b.v, c.v)}; }   // a*b + c
  GSH_INLINE V fnmadd(V a, V b, V c) { return {_mm512_fnmadd_pd(a.v, b.v, c.v)}; } // c - a*b
  GSH_INLINE V fmsub(V a, V b, V c) { return {_mm512_fmsub_pd(a.v, b.v, c.v)}; }   // a*b - c
  GSH_INLINE V abs(V a) { return {_mm512_abs_pd(a.v)}; }
  GSH_INLINE M gt(V a, V b) { return _mm512_cmp_pd_mask(a.v, b.v, _CMP_GT_OQ); }
  GSH_INLINE M lt(V a, V b) { return _mm512_cmp_pd_mask(a.v, b.v, _CMP_LT_OQ); }
  GSH_INLINE M ne(V a, V b) { return _mm512_cmp_pd_mask(a.v, b.v, _CMP_NEQ_UQ); }
  GSH_INLINE M mand(M a, M b) { return static_cast<M>(a & b); }
  GSH_INLINE M mor(M a, M b) { return static_cast<M>(a | b); }
  GSH_INLINE bool any(M a) { return a != 0; }
  GSH_INLINE V blend(M m, V f, V t) { return {_mm512_mask_blend_pd(m, f.v, t.v)}; } // m ? t : f
}
// AVX2 + FMA: V wraps __m256d, WIDTH = 4. Masks are full-width double vectors (all bits set or clear per lane),
// so any() uses movemask and blend uses blendv. abs clears the sign bit by and-not with -0.0.
#elif defined(GSH_ISA_AVX2)
#include <immintrin.h>
#define GSH_ISA_NS isa_avx2
namespace gsh::GSH_ISA_NS
{
  struct V
  {
    __m256d v;
  };
  using M = __m256d;
  inline constexpr int WIDTH = 4;
  GSH_INLINE V set1(double a) { return {_mm256_set1_pd(a)}; }
  GSH_INLINE V load(const double *p) { return {_mm256_loadu_pd(p)}; }
  GSH_INLINE void store(double *p, V a) { _mm256_storeu_pd(p, a.v); }
  GSH_INLINE V add(V a, V b) { return {_mm256_add_pd(a.v, b.v)}; }
  GSH_INLINE V sub(V a, V b) { return {_mm256_sub_pd(a.v, b.v)}; }
  GSH_INLINE V mul(V a, V b) { return {_mm256_mul_pd(a.v, b.v)}; }
  GSH_INLINE V div(V a, V b) { return {_mm256_div_pd(a.v, b.v)}; }
  GSH_INLINE V fmadd(V a, V b, V c) { return {_mm256_fmadd_pd(a.v, b.v, c.v)}; }
  GSH_INLINE V fnmadd(V a, V b, V c) { return {_mm256_fnmadd_pd(a.v, b.v, c.v)}; }
  GSH_INLINE V fmsub(V a, V b, V c) { return {_mm256_fmsub_pd(a.v, b.v, c.v)}; }
  GSH_INLINE V abs(V a) { return {_mm256_andnot_pd(_mm256_set1_pd(-0.0), a.v)}; }
  GSH_INLINE M gt(V a, V b) { return _mm256_cmp_pd(a.v, b.v, _CMP_GT_OQ); }
  GSH_INLINE M lt(V a, V b) { return _mm256_cmp_pd(a.v, b.v, _CMP_LT_OQ); }
  GSH_INLINE M ne(V a, V b) { return _mm256_cmp_pd(a.v, b.v, _CMP_NEQ_UQ); }
  GSH_INLINE M mand(M a, M b) { return _mm256_and_pd(a, b); }
  GSH_INLINE M mor(M a, M b) { return _mm256_or_pd(a, b); }
  GSH_INLINE bool any(M a) { return _mm256_movemask_pd(a) != 0; }
  GSH_INLINE V blend(M m, V f, V t) { return {_mm256_blendv_pd(f.v, t.v, m)}; }
}
// Portable scalar fallback: V wraps one double, WIDTH = 1, masks are plain bool. It has the same semantics as the
// SIMD sets, so the same kernel template runs on any CPU. Tests compare every ISA against reference_point().
#elif defined(GSH_ISA_SCALAR)
#include <cmath>
#define GSH_ISA_NS isa_scalar
namespace gsh::GSH_ISA_NS
{
  struct V
  {
    double v;
  };
  using M = bool;
  inline constexpr int WIDTH = 1;
  GSH_INLINE V set1(double a) { return {a}; }
  GSH_INLINE V load(const double *p) { return {*p}; }
  GSH_INLINE void store(double *p, V a) { *p = a.v; }
  GSH_INLINE V add(V a, V b) { return {a.v + b.v}; }
  GSH_INLINE V sub(V a, V b) { return {a.v - b.v}; }
  GSH_INLINE V mul(V a, V b) { return {a.v * b.v}; }
  GSH_INLINE V div(V a, V b) { return {a.v / b.v}; }
  GSH_INLINE V fmadd(V a, V b, V c) { return {a.v * b.v + c.v}; }
  GSH_INLINE V fnmadd(V a, V b, V c) { return {c.v - a.v * b.v}; }
  GSH_INLINE V fmsub(V a, V b, V c) { return {a.v * b.v - c.v}; }
  GSH_INLINE V abs(V a) { return {std::fabs(a.v)}; }
  GSH_INLINE M gt(V a, V b) { return a.v > b.v; }
  GSH_INLINE M lt(V a, V b) { return a.v < b.v; }
  GSH_INLINE M ne(V a, V b) { return a.v != b.v; }
  GSH_INLINE M mand(M a, M b) { return a && b; }
  GSH_INLINE M mor(M a, M b) { return a || b; }
  GSH_INLINE bool any(M a) { return a; }
  GSH_INLINE V blend(M m, V f, V t) { return m ? t : f; }
}
#else
#error "Define GSH_ISA_SCALAR, GSH_ISA_AVX2 or GSH_ISA_AVX512"
#endif
