// core.cpp: implementation of the geodetic-to-spherical conversion declared in core.hpp.
#include "core.hpp"
#include <cmath>

namespace gsh
{
  // Geodetic (phi, lambda, h) to geocentric spherical (phi', lambda, r) on the ellipsoid (a, e2).
  //
  // With phi = geodetic latitude, lambda = longitude, N = prime-vertical radius of curvature:
  //   N = a / sqrt(1 - e2 sin^2 phi)
  //   p = (N + h) cos phi                  (distance from the Z axis)
  //   Z = [N (1 - e2) + h] sin phi         (height above the equatorial plane)
  //   phi' = atan2(Z, p)                   (geocentric latitude)
  //   r   = sqrt(p^2 + Z^2)                (geocentric radius)
  // Formulas: Bucha and Janák (2013), which
  // gives X = p cos(lambda), Y = p sin(lambda), Z as above and r = sqrt(X^2 + Y^2 + Z^2).
  Spherical geodetic_to_spherical(double lat_deg, double lon_deg, double h, double a, double e2)
  {
    const double phi = deg2rad(lat_deg);
    const double s = std::sin(phi), c = std::cos(phi);
    const double N = a / std::sqrt(1.0 - e2 * s * s);
    const double p = (N + h) * c;
    const double z = (N * (1.0 - e2) + h) * s;
    return {std::atan2(z, p), deg2rad(lon_deg), std::hypot(p, z)};
  }
}
