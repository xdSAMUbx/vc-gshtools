// legendre.hpp: fully normalised (4*pi) associated Legendre functions P̄_nm(t) and their latitude derivatives.
//
// Conventions:
//   - Arguments: t = sin(φ) = cos(θ) and u = cos(φ) = sin(θ), where φ is the GEOCENTRIC (spherical) latitude
//     and θ = 90° - φ the geocentric colatitude. Geodetic latitudes must be converted to geocentric before
//     calling these routines (see Ellipsoid::geodetic_of_spherical in ellipsoid.hpp).
//   - Normalisation: P̄_nm = sqrt((2 - δ_0m)(2n+1)(n-m)!/(n+m)!) P_nm (fully normalised, 4*pi). This is NOT
//     the Schmidt semi-normalisation. The Condon-Shortley phase is NOT included (geodetic convention).
//   - Storage: a triangle packed by columns m. Element (n, m), 0 <= m <= n <= nmax, is at
//     idx_tri(n, m, nmax) = col_offset(m, nmax) + (n - m). A buffer for degree nmax needs tri_size(nmax)
//     doubles. All values are dimensionless.
//   - Derivatives are taken with respect to the spherical latitude φ, in radians (per radian).
//
// Algorithms (the one used by compute() is chosen by degree):
//   - SFC  (Standard Forward Column), Holmes and Featherstone (2002): column-by-column three-term recurrence in plain
//     doubles. Valid up to nmax ≈ 1800; above that the sectorial term underflows (measured in this project, benchmark/alf).
//   - MSFC (Modified SFC), after the X-numbers of Fukushima (2012): each value is stored as mantissa * R^-k with R = 2^500, which removes
//     underflow. It gives the same values as SFC wherever SFC is valid. Used for ultra-high degree.
//   - ERA  (Extended-Range Arithmetic), Smith, Olver and Lozier (1981), as used by Wittwer et al. (2008): every number carries its own integer
//     exponent. It has no scaling factors and works for all latitudes and degrees, but it is slow. It is the
//     precision reference used by the tests, not a production path.
//
// Recurrence formulas (Holmes and Featherstone 2002):
//   P̄_00 = 1; P̄_mm = sqrt((2m+1)/(2m)) u P̄_{m-1,m-1} (m >= 2), P̄_11 = sqrt(3) u;
//   P̄_{m+1,m} = sqrt(2m+3) t P̄_mm; P̄_nm = a_nm t P̄_{n-1,m} - b_nm P̄_{n-2,m} (n >= m+2).
#pragma once
#include <span>
#include "core.hpp"

namespace gsh::legendre
{
  // Recurrence coefficients for one maximum degree, stored in the same packed triangle layout as P.
  //   a(n,m) = sqrt((2n+1)(2n-1)/((n-m)(n+m)))                          for n >= m+1 (a(m,m) stays 0)
  //   b(n,m) = sqrt((2n+1)(n+m-1)(n-m-1)/((n-m)(n+m)(2n-3)))            for n >= m+2 (0 otherwise)
  //   f(n,m) = sqrt((n-m)(n+m)(2n+1)/(2n-1))                            used by the derivative (dP/dφ) column kernels
  //   d[0] = 1, d[1] = sqrt(3), d[m] = sqrt((2m+1)/(2m)) for m >= 2     sectorial scale factors
  struct Coefs
  {
    int nmax = -1;
    AlignedVec<double> a, b, f; // packed triangles, see make_coefs()
    AlignedVec<double> d;       // one entry per order m, size nmax + 1
  };

  // Builds the coefficients for degree nmax. Throws std::invalid_argument if nmax < 0.
  // Cost O(nmax^2) time and memory; the loop is parallel over m (OpenMP, dynamic chunks of 16 columns).
  Coefs make_coefs(int nmax);

  // Standard forward column (SFC), Holmes and Featherstone (2002). Fills P column by column: sectorial P(m,m), then the first
  // subdiagonal P(m+1,m), then the upward recurrence in n. Accurate to about 1e-13 up to nmax = 1800; beyond that
  // the sectorial term u^m underflows in double precision (Holmes and Featherstone (2002) measures a sum-rule error of 3e-4 at
  // nmax = 2000 and 0.46 at nmax = 2700, both for u ≈ 0.34). Use msfc() or compute() for larger degrees.
  //   t = sin φ, u = cos φ (geocentric), c = make_coefs(nmax), P has at least tri_size(c.nmax) elements.
  //   Throws std::invalid_argument if P is too small. At a pole (u = 0) all m >= 1 entries are exactly 0, no NaN.
  void sfc(double t, double u, const Coefs &c, std::span<double> P);

  // Modified standard forward column (MSFC) (the scaling idea is the one of Fukushima's X-numbers, Fukushima 2012).
  // Same recurrence as sfc(), but each value is kept
  // as a mantissa x with an implicit integer exponent k: value = x * R^-k, R = 2^500 (R_BIG / R_INV in core.hpp).
  // The mantissa is rescaled up when it becomes small and the last two column values are rescaled down when a
  // term grows large. Values with k >= 3 are returned as 0 (magnitude about 2^-1000 ≈ 1e-301). Same arguments and exceptions as sfc().
  void msfc(double t, double u, const Coefs &c, std::span<double> P);

  // Extended-range arithmetic (ERA) evaluation of P, Smith, Olver and Lozier (1981), ACM Trans. Math. Softw. 7, 93-105
  // (Lozier and Smith 1981, Algorithm 567, ACM TOMS 7, 141-146). Each number is an (exponent, mantissa) pair renormalised after
  // every operation. No scale factors, valid at every latitude and degree, but slow: a precision reference only.
  // GrafLab selects the same recurrences with its option fnALFs (1 sfc, 2 msfc, 3 era; Bucha and Janák 2013).
  // Same arguments and exceptions as sfc().
  void era(double t, double u, const Coefs &c, std::span<double> P);

  // Dispatcher: sfc() for c.nmax <= 1800, msfc() otherwise. The 1800 limit is the SFC validity bound measured in
  // this project (benchmark/alf).
  void compute(double t, double u, const Coefs &c, std::span<double> P);

  // Latitude derivative dP/dφ of a packed triangle P (same layout, same nmax), written to dP.
  // Uses the non-singular identities of Bucha and Janák (2019) (equivalent to those of Bosch 2000 and Fukushima 2012, II):
  //   m = 0:  dP(n,0)/dφ = sqrt(n(n+1)/2) P(n,1)
  //   m >= 1: dP(n,m)/dφ = 1/2 [ -sqrt((n+m)(n-m+1)) sqrt(1+δ_m1) P(n,m-1) + sqrt((n-m)(n+m+1)) P(n,m+1) ]
  //   with P(n, n+1) = 0. There is no division by cos φ, so the result is finite at the poles.
  // Applying dphi twice gives d²P/dφ². The ODE form P'' = tanφ P' - (n(n+1) - m²/cos²φ) P is singular at the
  // poles and is deliberately not used. Parallel over m (OpenMP) when nmax > 256.
  // Throws std::invalid_argument if either buffer holds fewer than tri_size(nmax) elements.
  void dphi(int nmax, std::span<const double> P, std::span<double> dP);
}
