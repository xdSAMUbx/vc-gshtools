// coordinates.cpp: geodetic, ECEF and geocentric spherical conversions. See coordinates.hpp for the
// conventions, formulas and sources.
#include "coordinates.hpp"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace gsh::coords
{
  namespace
  {
    // Upper bound on the fixed-point iterations of ecef_to_geodetic. The iteration converges quickly
    // from the Bowring start, so the bound is only a safeguard against non-convergence.
    constexpr int MAX_ITER = 20;
    // Convergence threshold on the latitude update, in radians (about 6e-9 m at the Earth surface).
    constexpr double PHI_TOL = 1e-15;
  }

  // Standard closed form given in coordinates.hpp (section "Formulas").
  Ecef geodetic_to_ecef(double lat_deg, double lon_deg, double h, const Ellipsoid &ell)
  {
    const double phi = deg2rad(lat_deg), lam = deg2rad(lon_deg);
    const double s = std::sin(phi), c = std::cos(phi);
    const double N = ell.a / std::sqrt(1.0 - ell.e2 * s * s); // prime-vertical radius of curvature [m]
    return {(N + h) * c * std::cos(lam),
            (N + h) * c * std::sin(lam),
            (N * (1.0 - ell.e2) + h) * s};
  }

  // Iterative ECEF to geodetic conversion (algorithm in coordinates.hpp).
  Geodetic ecef_to_geodetic(const Ecef &q, const Ellipsoid &ell)
  {
    const double a = ell.a, e2 = ell.e2;
    // Semi-minor axis b = a sqrt(1 - e^2) and the second eccentricity squared e'^2 = e^2 / (1 - e^2).
    const double b = a * std::sqrt(1.0 - e2);
    const double ep2 = e2 / (1.0 - e2);

    const double p = std::hypot(q.x, q.y);  // distance from the Z axis [m]
    const double lam = std::atan2(q.y, q.x); // longitude [rad], in (-pi, pi]

    // Step 2: closed-form starting latitude (Bowring-type; see coordinates.hpp). theta is the
    // parametric (reduced) latitude of the point.
    const double theta = std::atan2(q.z * a, p * b);
    const double st = std::sin(theta), ct = std::cos(theta);
    double phi = std::atan2(q.z + ep2 * b * st * st * st, p - e2 * a * ct * ct * ct);

    // Step 3: fixed-point iteration phi = atan2(Z + e^2 N sin(phi), p).
    for (int it = 0; it < MAX_ITER; ++it)
    {
      const double s = std::sin(phi);
      const double N = a / std::sqrt(1.0 - e2 * s * s);
      const double next = std::atan2(q.z + e2 * N * s, p);
      const double change = std::fabs(next - phi);
      phi = next;
      if (change < PHI_TOL)
        break;
    }

    // Step 4: height as the signed distance along the ellipsoidal normal.
    const double s = std::sin(phi), c = std::cos(phi);
    const double h = p * c + q.z * s - a * std::sqrt(1.0 - e2 * s * s);
    return {rad2deg(phi), rad2deg(lam), h};
  }

  // phi' = atan2((1 - e^2) sin phi, cos phi): same as atan(tan(phi) (1 - e^2)) but valid at +-90 deg.
  double geodetic_to_geocentric_latitude(double lat_deg, const Ellipsoid &ell)
  {
    const double phi = deg2rad(lat_deg);
    return rad2deg(std::atan2((1.0 - ell.e2) * std::sin(phi), std::cos(phi)));
  }

  // Reuses the inverse relation already implemented in the Ellipsoid class.
  double geocentric_to_geodetic_latitude(double lat_deg, const Ellipsoid &ell)
  {
    return rad2deg(ell.geodetic_of_spherical(deg2rad(lat_deg)));
  }

  // Reuses gsh::geodetic_to_spherical (core.cpp). Only the output units are converted here.
  SphericalDeg geodetic_to_spherical(double lat_deg, double lon_deg, double h, const Ellipsoid &ell)
  {
    const Spherical s = gsh::geodetic_to_spherical(lat_deg, lon_deg, h, ell.a, ell.e2);
    return {rad2deg(s.lat), lon_deg, s.r};
  }

  // ---- Vector functions -----------------------------------------------------------------------
  // Each loop body calls the scalar function for one index. The index type is signed (std::ptrdiff_t),
  // which OpenMP 2.0 (MSVC) requires. Size checks happen before the parallel region, so nothing
  // throws inside it.

  std::vector<Ecef> geodetic_to_ecef(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg,
                                     const std::vector<double> &h, const Ellipsoid &ell)
  {
    if (lat_deg.size() != lon_deg.size() || lat_deg.size() != h.size())
      throw std::invalid_argument("coords::geodetic_to_ecef: input vectors have different sizes");
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(lat_deg.size());
    std::vector<Ecef> out(lat_deg.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < n; ++i)
      out[i] = geodetic_to_ecef(lat_deg[i], lon_deg[i], h[i], ell);
    return out;
  }

  std::vector<Geodetic> ecef_to_geodetic(const std::vector<Ecef> &p, const Ellipsoid &ell)
  {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(p.size());
    std::vector<Geodetic> out(p.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < n; ++i)
      out[i] = ecef_to_geodetic(p[i], ell);
    return out;
  }

  std::vector<double> geodetic_to_geocentric_latitude(const std::vector<double> &lat_deg, const Ellipsoid &ell)
  {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(lat_deg.size());
    std::vector<double> out(lat_deg.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < n; ++i)
      out[i] = geodetic_to_geocentric_latitude(lat_deg[i], ell);
    return out;
  }

  std::vector<double> geocentric_to_geodetic_latitude(const std::vector<double> &lat_deg, const Ellipsoid &ell)
  {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(lat_deg.size());
    std::vector<double> out(lat_deg.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < n; ++i)
      out[i] = geocentric_to_geodetic_latitude(lat_deg[i], ell);
    return out;
  }

  std::vector<SphericalDeg> geodetic_to_spherical(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg,
                                                  const std::vector<double> &h, const Ellipsoid &ell)
  {
    if (lat_deg.size() != lon_deg.size() || lat_deg.size() != h.size())
      throw std::invalid_argument("coords::geodetic_to_spherical: input vectors have different sizes");
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(lat_deg.size());
    std::vector<SphericalDeg> out(lat_deg.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < n; ++i)
      out[i] = geodetic_to_spherical(lat_deg[i], lon_deg[i], h[i], ell);
    return out;
  }
}
