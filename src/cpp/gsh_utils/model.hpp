// model.hpp: spherical-harmonic coefficient model (ICGEM .gfc format) and its reader.
//
// Conventions:
//   - Coefficients are fully normalized (the reader rejects any other "norm" value).
//   - C and S are stored in a triangular layout: index idx_tri(n, m, nmax) (core.hpp).
//   - GM in m^3/s^2, R (reference radius) in metres.
#pragma once
#include <string>
#include "core.hpp"

namespace gsh
{
  // Harmonic coefficient set up to degree nmax. Coefficients not present in the file are zero.
  struct Model
  {
    double GM = 0.0;    // gravitational constant [m^3/s^2]
    double R = 0.0;     // reference radius [m]
    int nmax = -1;      // maximum degree stored in C and S (may be below file_nmax)
    int file_nmax = -1; // max_degree declared in the file header
    std::string name, norm, tide_system; // header fields: modelname, norm, tide_system (may be empty)
    AlignedVec<double> C, S; // packed with idx_tri(n, m, nmax); C is fully normalized C̄, S is S̄

    // Accessors without bounds checking. Valid for 0 <= m <= n <= nmax.
    double c(int n, int m) const { return C[idx_tri(n, m, nmax)]; }
    double s(int n, int m) const { return S[idx_tri(n, m, nmax)]; }

    // Model with all coefficients zero, degree nmax, and the given GM and R.
    static Model zeros(double GM, double R, int nmax);
    // Copy with the maximum degree reduced to new_nmax (coefficients of degree <= new_nmax are kept).
    // Throws std::invalid_argument if new_nmax > nmax.
    Model truncated(int new_nmax) const;
  };

  // Reads an ICGEM .gfc file in parallel. The header (up to the line starting with "end_of_head")
  // provides earth_gravity_constant, radius, max_degree, and optionally modelname, norm and
  // tide_system. Data lines start with "gfc" followed by n, m, C, S (and two extra fields that are
  // ignored). Fortran exponents such as 1.0d0 are accepted.
  //   nmax < 0: keep all degrees up to max_degree; otherwise read only degrees <= nmax.
  // Coefficients absent from the file stay zero.
  // Throws std::runtime_error for a missing file, an incomplete or invalid header, an unsupported
  // norm, or malformed coefficient lines; throws std::invalid_argument if nmax > max_degree.
  Model load_gfc(const std::string &path, int nmax = -1);

  // Same reader, but C and S of the returned model hold the formal standard deviations sigma C and sigma S
  // (the 5th and 6th columns of the "gfc" lines) instead of the coefficients. GM, R and the header are as in
  // load_gfc. Throws std::runtime_error if a coefficient line has no sigma columns.
  Model load_gfc_sigma(const std::string &path, int nmax = -1);
}
