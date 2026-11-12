// kernels_avx2.cpp: AVX2 + FMA build of the column kernels. Defines the lane type (simd.hpp) and instantiates
// column_kernel.inl for this ISA; it is compiled with the AVX2 instruction-set flags of the CMake target.
#define GSH_ISA_AVX2
#include "simd.hpp"
#include "column_kernel.inl"

const gsh::kernels::Table &gsh::kernels::table_avx2()
{
  static const Table t = GSH_ISA_NS::make_table("avx2");
  return t;
}
