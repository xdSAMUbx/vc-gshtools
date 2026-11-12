// coordinates.hpp: conversions between geodetic coordinates on a reference ellipsoid (latitude phi,
// longitude lambda, ellipsoidal height h), geocentric Cartesian coordinates (ECEF X, Y, Z), and
// geocentric spherical coordinates (geocentric latitude phi', longitude lambda, radius r).
//
// Conventions:
//   - Angles are in DEGREES at this API (inputs and outputs), like the rest of the gsh public API.
//     Internal trigonometry uses radians. The one exception is gsh::Spherical (core.hpp), which
//     stores radians; this module does not return it.
//   - Lengths are metres. The ECEF origin is the Earth's centre, Z points to the north pole, X lies
//     in the equatorial plane at lambda = 0 (IERS Conventions 2010, IERS Technical Note 36).
//   - Longitude is returned in (-180, 180] by ecef_to_geodetic and is returned unchanged by the
//     forward functions.
//   - The default reference ellipsoid is GRS80 (Ellipsoid::GRS80(), ellipsoid.hpp).
//
// Formulas (geodetic to ECEF, with N the prime-vertical radius of curvature):
//   N = a / sqrt(1 - e^2 sin^2 phi)
//   X = (N + h) cos phi cos lambda
//   Y = (N + h) cos phi sin lambda
//   Z = (N (1 - e^2) + h) sin phi
// Source: Bucha and Janák (2013).
//
// Geocentric latitude: tan(phi') = (1 - e^2) tan(phi), derived from the formulas above.
//
// Vector overloads (std::vector) use OpenMP (`parallel for`, static schedule) and call the scalar
// function element by element, so they return bit-identical results to the scalar functions.
#pragma once
#include <vector>
#include "core.hpp"
#include "ellipsoid.hpp"

namespace gsh::coords
{
  // Geocentric Cartesian (ECEF) position in metres.
  struct Ecef
  {
    double x, y, z;
  };

  // Geodetic position: latitude and longitude in degrees, ellipsoidal height in metres.
  struct Geodetic
  {
    double lat_deg, lon_deg, h;
  };

  // Geocentric spherical position: geocentric latitude and longitude in degrees, radius in metres.
  struct SphericalDeg
  {
    double lat_deg, lon_deg, r;
  };

  // ---- Scalar functions -----------------------------------------------------------------------

  // Geodetic (lat_deg, lon_deg, h) to ECEF. Closed form; no iteration; no error is possible for
  // finite inputs. Inputs: latitude and longitude in degrees, h in metres.
  Ecef geodetic_to_ecef(double lat_deg, double lon_deg, double h, const Ellipsoid &ell = Ellipsoid::GRS80());

  // ECEF to geodetic. Returns latitude in [-90, 90] degrees, longitude in (-180, 180] degrees and the
  // ellipsoidal height in metres.
  //
  // Algorithm (iterative, converges in a few steps to machine precision for any point, including
  // the poles and points at large heights):
  //   1. p = sqrt(X^2 + Y^2), lambda = atan2(Y, X).
  //   2. Initial latitude from a closed-form approximation of Bowring's type (Bowring 1976, Eq. 4,
  //      Survey Review 23(181), 323-327):
  //        theta = atan2(Z a, p b),
  //        phi0  = atan2(Z + (e'^2 b) sin^3 theta, p - e^2 a cos^3 theta),  e'^2 = e^2 / (1 - e^2).
  //   3. Fixed-point iteration, which is exact at the solution because
  //        Z + e^2 N sin(phi) = (N + h) sin(phi)  and  p = (N + h) cos(phi):
  //        phi_{k+1} = atan2(Z + e^2 N(phi_k) sin(phi_k), p),  N(phi) = a / sqrt(1 - e^2 sin^2 phi).
  //      Stops when |phi_{k+1} - phi_k| < 1e-15 rad (at most 20 steps).
  //   4. h = p cos(phi) + Z sin(phi) - a sqrt(1 - e^2 sin^2 phi).
  // At the poles p = 0 and atan2 gives +-90 degrees exactly; the formulas stay valid there.
  // The rounding error of h is limited by the representation of X, Y, Z (about 1e-9 m at |r| ~ 6.4e6 m).
  Geodetic ecef_to_geodetic(const Ecef &p, const Ellipsoid &ell = Ellipsoid::GRS80());

  // Geodetic latitude (degrees) to geocentric latitude (degrees): phi' = atan2((1 - e^2) sin phi, cos phi).
  // The atan2 form is exact at +-90 degrees (no division by cos phi).
  double geodetic_to_geocentric_latitude(double lat_deg, const Ellipsoid &ell = Ellipsoid::GRS80());
  // Geocentric latitude (degrees) to geodetic latitude (degrees): phi = atan2(sin phi', (1 - e^2) cos phi').
  // Implemented with Ellipsoid::geodetic_of_spherical (ellipsoid.hpp), the inverse of the same relation.
  double geocentric_to_geodetic_latitude(double lat_deg, const Ellipsoid &ell = Ellipsoid::GRS80());

  // Geodetic (lat_deg, lon_deg, h) to geocentric spherical (phi', lambda, r). The result is computed
  // by gsh::geodetic_to_spherical (core.hpp), so the formula is not duplicated. The longitude is
  // returned as given. r is the geocentric radius in metres.
  SphericalDeg geodetic_to_spherical(double lat_deg, double lon_deg, double h, const Ellipsoid &ell = Ellipsoid::GRS80());

  // ---- Vector functions (OpenMP) --------------------------------------------------------------
  // Element-wise versions of the functions above. Functions with several input vectors throw
  // std::invalid_argument if the sizes differ. The output has the size of the input vector.

  std::vector<Ecef> geodetic_to_ecef(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg,
                                     const std::vector<double> &h, const Ellipsoid &ell = Ellipsoid::GRS80());
  std::vector<Geodetic> ecef_to_geodetic(const std::vector<Ecef> &p, const Ellipsoid &ell = Ellipsoid::GRS80());
  std::vector<double> geodetic_to_geocentric_latitude(const std::vector<double> &lat_deg,
                                                      const Ellipsoid &ell = Ellipsoid::GRS80());
  std::vector<double> geocentric_to_geodetic_latitude(const std::vector<double> &lat_deg,
                                                      const Ellipsoid &ell = Ellipsoid::GRS80());
  std::vector<SphericalDeg> geodetic_to_spherical(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg,
                                                  const std::vector<double> &h, const Ellipsoid &ell = Ellipsoid::GRS80());
}
