// ellipsoid.hpp: reference ellipsoid (GRS80 by default), its normal gravity field and the
// coefficients that remove the ellipsoidal (normal) part from a global geopotential model.
//
// Units: lengths in metres, GM in m^3/s^2, omega in rad/s. Angles in the normal-field functions
// are in radians unless the name says "_deg".
//
// Sources:
//   - GRS80 constants: Moritz, H. (2000), "Geodetic Reference System 1980", Journal of Geodesy 74,
//     128-133.
//   - Zonal coefficients and disturbing coefficients: Bucha and Janák (2019) (the "PDF n" labels in the comments
//     refer to the pages of that document).
#pragma once
#include "core.hpp"
#include "model.hpp"

namespace gsh
{
  // Reference ellipsoid of revolution with its gravitational constants. Default values are GRS80.
  struct Ellipsoid
  {
    double a = 6378137.0;           // semi-major axis [m]
    double e2 = 0.00669438002290;   // first eccentricity squared [-]
    double GM = 3.986005e14;        // gravitational constant of the ellipsoid [m^3/s^2]
    double J2 = 1.08263e-3;         // dynamical form factor J2 [-]
    double omega = 7.292115e-5;     // angular velocity [rad/s]

    // Returns the GRS80 ellipsoid (the default-constructed values).
    static Ellipsoid GRS80() { return {}; }

    // Fully normalized zonal coefficient C̄^Ell_{2k,0} of the ellipsoid, k = 0..10. Odd-degree
    // coefficients and all S coefficients are zero. Formula in
    // Bucha and Janák (2019).
    // Throws std::out_of_range if k is outside 0..10.
    double zonal(int k) const;
    // Geocentric radius of the ellipsoid at geocentric latitude lat [rad]. Returns metres.
    // r(phi') = a sqrt((1 - e2) / (1 - e2 cos^2 phi')).
    double radius(double lat) const;
    // Geodetic latitude [rad] of the ellipsoid point whose geocentric latitude is lat [rad].
    // Inverse of the relation tan(phi') = (1 - e2) tan(phi). Used by gsh::coords.
    double geodetic_of_spherical(double lat) const;

    struct Normal
    {
      double g_x, g_z, gamma; // north component, up component, magnitude [m/s^2]
    };
    // Normal gravity at geocentric latitude lat [rad] and geocentric radius r [m], evaluated by
    // the spherical-harmonic series of the ellipsoid (gamma_SH). Includes the centrifugal term.
    Normal normal_gravity(double lat, double r) const;

    // Somigliana closed-form normal gravity [m/s^2] at geodetic latitude lat_geodetic [rad]:
    //   gamma(phi) = (a gamma_e cos^2 phi + b gamma_p sin^2 phi) / sqrt(a^2 cos^2 phi + b^2 sin^2 phi)
    // with b = a sqrt(1 - e2). Intended for validation only (Moritz 2000). Defaults are the GRS80
    // equatorial and polar normal gravity values.
    double somigliana(double lat_geodetic, double gamma_e = 9.7803267715, double gamma_p = 9.8321863685) const;
  };

  // Disturbing (anomalous) coefficients: copy of `model` with the ellipsoidal zonal part removed,
  //   dC_{2k,0} = C_{2k,0} - C^Ell_{2k,0} (GM_Ell / GM) (a_Ell / R)^(2k),  k = 0..10 and 2k <= nmax
  // (degree 0 included). Only these zonal terms are modified; the S coefficients are unchanged.
  // Formulas: Bucha and Janák (2019).
  Model disturbing_model(const Model &model, const Ellipsoid &ell);
}
