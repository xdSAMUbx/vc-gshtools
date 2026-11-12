// kernels_avx512.cpp: AVX-512 build of the column kernels. Defines the lane type (simd.hpp) and instantiates
// column_kernel.inl for this ISA; it is compiled with the AVX-512 instruction-set flags of the CMake target.
#define GSH_ISA_AVX512
#include "simd.hpp"
#include "column_kernel.inl"

const gsh::kernels::Table &gsh::kernels::table_avx512()
{
  static const Table t = GSH_ISA_NS::make_table("avx512");
  return t;
}
