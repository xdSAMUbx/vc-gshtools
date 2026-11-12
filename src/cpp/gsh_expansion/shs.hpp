// shs.hpp: spherical-harmonic synthesis (SHS) at scattered points.
// Returns the RAW SUMS of Bucha and Janák (2019) (A_n, A^φ_n, A^λ_n, ...,
// without the GM/r^k factor). functionals.hpp turns these sums into the geopotential functionals.
//
// Units and conventions shared by everything in this header:
//   - Latitudes are SPHERICAL (geocentric) and in radians; longitudes in radians; radii in metres.
//   - Coefficients are fully normalised (4*pi), dimensionless; GM in m^3 s^-2; R in metres.
//   - The raw sums are dimensionless. Field order is kernels::Field (kernels.hpp).
#pragma once
#include <memory>
#include <string>
#include <vector>
#include "core.hpp"
#include "legendre.hpp"
#include "model.hpp"
#include "synth_types.hpp"

namespace gsh
{
  // Recurrence tables for one maximum degree, shared between calls that use the same nmax (saves rebuilding them).
  //   ab: interleaved (a, b) per packed triangle index, ab[2i] = a, ab[2i+1] = b (legendre::Coefs);
  //   f:  factors for dP/dφ; d: sectorial factors d[m], nmax + 1 entries.
  struct Recurrence
  {
    int nmax = -1;
    AlignedVec<double> ab, f, d;
  };
  // Builds the recurrence tables for nmax (through legendre::make_coefs; throws std::invalid_argument for nmax < 0).
  std::shared_ptr<const Recurrence> make_recurrence(int nmax);

  // A model prepared for the kernels: coefficients truncated to [nmin, nmax] and repacked for nmax.
  struct Prepared
  {
    int nmin = 0, nmax = -1;
    double GM = 0.0, R = 0.0;             // copied from the model (m^3 s^-2 and m)
    std::shared_ptr<const Recurrence> rec; // recurrence tables, possibly shared between Prepared objects
    AlignedVec<double> cs;                // [C, S] per packed (n, m); zero for n < nmin
  };
  // Prepares model for synthesis. nmin < 0 is not allowed; nmax < 0 means model.nmax. Degrees below nmin are zeroed.
  // Pass an existing rec (for the same nmax) to share recurrence tables; nullptr builds them.
  // Throws std::invalid_argument if nmax > model.nmax, nmin is outside [0, nmax], or rec has another degree.
  Prepared prepare(const Model &model, int nmin = 0, int nmax = -1, std::shared_ptr<const Recurrence> rec = nullptr);

  // Raw sums per point, one array per field, in kernels::Field order. grad 0 -> 1 field, 1 -> 4, 2 -> 10.
  struct RawFields
  {
    int grad = 0;
    std::vector<AlignedVec<double>> f;
  };

  // Raw sums at every point, in input order. grad = 0, 1 or 2 (derivative order of the Legendre functions).
  // Uses the CPU kernels. Results are dimensionless.
  // With grad >= 1 the values are not reliable when cos φ is very small (the derivative recurrence divides by
  // cos φ); compute_functionals() handles those points with the oracle (reference_point).
  // Throws std::invalid_argument if grad is outside [0, 2]. Parallel over tiles of points with OpenMP.
  RawFields synth_points(const Prepared &prep, const Points &pts, int grad, const SynthOptions &opts = {});

  // Name of the instruction set that synth_points() would use for isa: "avx512", "avx2" or "scalar".
  // Isa::Auto picks the best one the CPU supports. Throws std::runtime_error if a forced ISA is not supported.
  std::string active_isa_name(Isa isa = Isa::Auto);

  // Sequential oracle at one point: full triangle computed with MSFC and the literal formulas of the specification,
  // with no singularity at the poles. Slow (O(nmax^2) per point); used for validation and for near-polar points.
  struct ReferencePoint
  {
    double raw[10];                     // raw sums, kernels::Field order (dimensionless)
    // Local tensor in the LNOF frame (x north, y west, z up), in s^-2, already multiplied by GM/r^3 (Bucha and Janák 2019).
    double vxx, vxy, vxz, vyy, vyz, vzz;
  };
  // Evaluates the oracle at spherical latitude lat (rad), longitude lon (rad) and radius r (m).
  // coefs must have coefs.nmax >= nmax, and nmax must not exceed model.nmax. Degrees below nmin contribute nothing.
  ReferencePoint reference_point(const Model &model, const legendre::Coefs &coefs, double lat, double lon, double r, int nmin = 0, int nmax = -1);
}
