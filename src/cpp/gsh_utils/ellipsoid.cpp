// ellipsoid.cpp: reference ellipsoid, zonal coefficients, normal gravity series and disturbing
// coefficients. See ellipsoid.hpp for units and sources.
#include "ellipsoid.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gsh
{
  // Zonal coefficient of degree 2k of the ellipsoid (fully normalized):
  //   k = 0: 1;   k = 1: C_20 = -J2 / sqrt(5);
  //   k >= 2: C_{2k,0} = (-1)^k 3 e^(2k) / [(2k+1)(2k+3) sqrt(4k+1)] * (1 - k - 5 sqrt(5) k C_20 / e^2)
  // Source: Bucha and Janák (2019); Moritz (2000).
  double Ellipsoid::zonal(int k) const
  {
    if (k < 0 || k > 10)
      throw std::out_of_range("Ellipsoid::zonal: k must be in 0..10");
    if (k == 0)
      return 1.0;
    const double c20 = -J2 / std::sqrt(5.0);
    if (k == 1)
      return c20;
    const double K = k;
    const double sign = (k % 2 == 0) ? 1.0 : -1.0;
    return sign * 3.0 * std::pow(e2, K) / ((2 * K + 1) * (2 * K + 3) * std::sqrt(4 * K + 1)) *
           (1.0 - K - 5.0 * std::sqrt(5.0) * K * c20 / e2);
  }

  // Geocentric radius of the reference ellipsoid: r = a sqrt((1 - e2) / (1 - e2 cos^2 lat)).
  // Formula: Bucha and Janák (2013).
  double Ellipsoid::radius(double lat) const
  {
    const double c = std::cos(lat);
    return a * std::sqrt((1.0 - e2) / (1.0 - e2 * c * c));
  }

  // Inverse of tan(phi') = (1 - e2) tan(phi): phi = atan2(sin phi', (1 - e2) cos phi').
  // atan2 keeps the result finite and correct at the poles (cos phi' = 0).
  double Ellipsoid::geodetic_of_spherical(double lat) const
  {
    return std::atan2(std::sin(lat), (1.0 - e2) * std::cos(lat));
  }

  // Normal gravity from the spherical-harmonic expansion of the normal potential (gamma_SH).
  //
  // Steps:
  //   1. Fully normalized associated Legendre functions P̄_{n,0} and P̄_{n,1} of t = sin(lat) and
  //      u = cos(lat), for n = 0..20. The P̄_{n,0} and P̄_{n,1} are built with the standard
  //      three-term recurrences; P̄_{n,1} is used in place of dP̄_{n,0}/dphi through the identity
  //      dP̄_{n,0}/dphi = sqrt(n(n+1)/2) P̄_{n,1}, which avoids the singularity at the poles.
  //   2. Radial and latitudinal sums over the even zonal coefficients C̄_{2k,0}, k = 0..10, with
  //      the factor (a/r)^(2k) accumulated in qn.
  //   3. Radial and latitudinal derivatives of the zonal potential sum (as coded in sr and sp):
  //        dVr = -GM/r^2 sum (n+1) C̄ P̄_{n,0}      dVp = GM/r sum C̄ sqrt(n(n+1)/2) P̄_{n,1}
  //   4. Centrifugal terms of the rotating frame: dVcr = omega^2 r cos^2(lat) (radial) and
  //      dVcp = -omega^2 r^2 sin(lat) cos(lat) (latitudinal). The north (g_x) and up (g_z)
  //      components combine these with dVp and dVr as written in the code; gamma is their norm.
  // Formulas: Bucha and Janák (2019) (normal gravity gamma_SH and zonal coefficients).
  Ellipsoid::Normal Ellipsoid::normal_gravity(double lat, double r) const
  {
    // At the exact poles cos(lat) = 0 and sin(lat) = +-1; handled explicitly to avoid rounding noise.
    const bool pole = std::fabs(lat) == PI / 2;
    const double t = pole ? (lat > 0 ? 1.0 : -1.0) : std::sin(lat);
    const double u = pole ? 0.0 : std::cos(lat);
    const double q = a / r;

    // Columns m = 0 and m = 1 up to degree 20; dP_n0/dphi = sqrt(n(n+1)/2) P_n1 (no singularity).
    double p0[21], p1[21];
    p0[0] = 1.0;
    p0[1] = std::sqrt(3.0) * t;
    p1[0] = 0.0;
    p1[1] = std::sqrt(3.0) * u;
    for (int n = 2; n <= 20; ++n)
    {
      const double N = n;
      p0[n] = std::sqrt(4 * N * N - 1) / N * t * p0[n - 1] - (N - 1) / N * std::sqrt((2 * N + 1) / (2 * N - 3)) * p0[n - 2];
      const double a1 = std::sqrt((2 * N + 1) * (2 * N - 1) / ((N - 1) * (N + 1)));
      const double b1 = n >= 3 ? std::sqrt((2 * N + 1) * N * (N - 2) / ((N - 1) * (N + 1) * (2 * N - 3))) : 0.0;
      p1[n] = a1 * t * p1[n - 1] - b1 * p1[n - 2];
    }

    double sr = 0, sp = 0, qn = 1.0;
    const double q2 = q * q;
    for (int k = 0; k <= 10; ++k)
    {
      const int n = 2 * k;
      const double N = n, c = zonal(k);
      sr += qn * (N + 1) * c * p0[n];
      sp += qn * c * std::sqrt(N * (N + 1) / 2) * p1[n];
      qn *= q2;
    }
    const double dVr = -GM / (r * r) * sr;
    const double dVp = GM / r * sp;
    const double dVcr = omega * omega * r * u * u;
    const double dVcp = -omega * omega * r * r * u * t;
    Normal out;
    out.g_x = (dVp + dVcp) / r;
    out.g_z = dVr + dVcr;
    out.gamma = std::hypot(out.g_x, out.g_z);
    return out;
  }

  // Somigliana closed-form normal gravity on the surface of the ellipsoid, at geodetic latitude phi:
  //   gamma(phi) = (a gamma_e cos^2 phi + b gamma_p sin^2 phi) / sqrt(a^2 cos^2 phi + b^2 sin^2 phi)
  // with b = a sqrt(1 - e2). Used only as an independent check of normal_gravity().
  // Source: Moritz (2000), Somigliana's closed form.
  double Ellipsoid::somigliana(double phi, double ge, double gp) const
  {
    const double s = std::sin(phi), c = std::cos(phi);
    const double b = a * std::sqrt(1.0 - e2);
    return (a * ge * c * c + b * gp * s * s) / std::sqrt(a * a * c * c + b * b * s * s);
  }

  // Zonal terms removed from the model, delta_k = C^Ell_{2k,0} (GM_Ell/GM)(a_Ell/R)^(2k), k = 0..10, 2k <= min(20, nmax).
  // This is the only place where delta_k is computed; disturbing_model subtracts it and compute_functionals uses it as
  // the model it adds back. `sn` carries (GM_Ell/GM)(a/R)^(2k) between steps.
  Model normal_zonal_model(const Model &model, const Ellipsoid &ell)
  {
    const int nz = std::min(20, model.nmax);
    Model z = Model::zeros(model.GM, model.R, nz);
    z.file_nmax = model.file_nmax;
    z.name = model.name;
    z.norm = model.norm;
    z.tide_system = model.tide_system;
    const double ratio = ell.GM / model.GM;
    const double scale = ell.a / model.R;
    double sn = ratio; // (GM_ell/GM) (a/R)^n
    for (int k = 0; k <= 10 && 2 * k <= nz; ++k)
    {
      z.C[idx_tri(2 * k, 0, nz)] = ell.zonal(k) * sn;
      sn *= scale * scale;
    }
    return z;
  }

  // Disturbing coefficients: the ellipsoidal zonal terms are scaled to the GM and R of the model,
  //   dC_{2k,0} = C_{2k,0} - C^Ell_{2k,0} (GM_Ell / GM) (a_Ell / R)^(2k).
  // The loop runs over k = 0..10 and 2k <= nmax (the subtracted terms come from normal_zonal_model).
  // Formula: Bucha and Janák (2019).
  Model disturbing_model(const Model &model, const Ellipsoid &ell)
  {
    Model d = model.truncated(model.nmax);
    const Model z = normal_zonal_model(model, ell);
    for (int k = 0; k <= 10 && 2 * k <= model.nmax; ++k)
      d.C[idx_tri(2 * k, 0, model.nmax)] -= z.C[idx_tri(2 * k, 0, z.nmax)];
    return d;
  }
}
