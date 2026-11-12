// kernels_scalar.cpp: portable scalar build of the column kernels. Used when no SIMD ISA is available or requested.
#define GSH_ISA_SCALAR
#include "simd.hpp"
#include "column_kernel.inl"

const gsh::kernels::Table &gsh::kernels::table_scalar()
{
  static const Table t = GSH_ISA_NS::make_table("scalar");
  return t;
}
