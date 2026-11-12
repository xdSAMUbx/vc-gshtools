// tides.cpp: conversions between the permanent-tide systems tide-free (NT), zero-tide (ZT) and
// mean-tide (MT). The operation order follows the original reference implementation (not part of this
// release), so that the values reproduce its results; tests/cpp/test_tides.cpp compares against reference
// tables generated with that implementation.
//
// Formulas and sources (see tides.hpp for the full list):
//   (1) dC̄20 = offset(to) - offset(from), offset(NT) = 0, offset(ZT) = k2 A0 H0.
//       Sources: the C̄20 shift proportional to k20 is in Makinen (2021) (Eqs. 28-30); the numerical values
//       of A0 and H0 are documented in tides.hpp.
//       MT has no C̄20 of its own (NO_C20_MT).
//   (2) D(phi) = A0 H0 (GM / (gamma(phi) R)) sqrt(5) P2(sin phi), the direct permanent term of the geoid
//       height. N(NT) - N(ZT) = -k D and N(MT) - N(ZT) = D. The per-unit-C̄20 factor
//       (GM / (gamma R)) sqrt(5) P2 is the degree-2 zonal geoid response. Sign and magnitude are
//       consistent with H_T(phi) of Makinen (2021) at the
//       equator and at the pole.
//   (3) h_T, H_T, g_T: polynomials in sin^2 phi and sin^4 phi (Makinen 2021).
//   (4) h(MT) = h(ZT) = h(NT) + h_T (IHRF convention; Makinen 2021).
#include "tides.hpp"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include "core.hpp"

namespace gsh::tides
{
  namespace
  {
    constexpr double PI = 3.141592653589793238462643383279502884;
    constexpr double DEG2RAD = PI / 180.0; // same value as numpy.radians

    // sin^2(phi) and sin^4(phi) for a latitude in degrees.
    inline void sin2_sin4(double lat_deg, double &s2, double &s4)
    {
      const double s = std::sin(lat_deg * DEG2RAD);
      s2 = s * s;
      s4 = s2 * s2;
    }

    // Legendre polynomial P2(sin phi) = (3 sin^2 phi - 1) / 2.
    inline double p2(double lat_deg)
    {
      const double s = std::sin(lat_deg * DEG2RAD);
      return (3.0 * s * s - 1.0) / 2.0;
    }

    // Normal gravity on the GRS80 surface (Somigliana), gamma(phi) [m/s^2]:
    //   gamma = gamma_e (1 + k sin^2 phi) / sqrt(1 - e2 sin^2 phi)
    // Constants: Moritz (2000), Geodetic Reference System 1980 (see tides.hpp).
    inline double gamma_normal(double lat_deg)
    {
      const double s = std::sin(lat_deg * DEG2RAD);
      const double s2 = s * s;
      return GAMMA_E * (1.0 + GRS80_K * s2) / std::sqrt(1.0 - GRS80_E2 * s2);
    }

    // Geoid response per unit of C̄20: (GM / (gamma R)) sqrt(5) P2(sin phi) [m].
    inline double geoid_per_c20(double lat_deg)
    {
      return (GM / (gamma_normal(lat_deg) * R)) * std::sqrt(5.0) * p2(lat_deg);
    }

    // D(phi): direct permanent term of the geoid height [m] (item (2) of the file header).
    inline double direct_geoid_term(double lat_deg)
    {
      return A0 * H0 * geoid_per_c20(lat_deg);
    }

    // Offset of C̄20 in system s relative to tide-free [dimensionless]. Throws for MeanTide.
    inline double c20_offset(TideSystem s, double k2)
    {
      switch (s)
      {
      case TideSystem::TideFree: return 0.0;
      case TideSystem::ZeroTide: return k2 * A0 * H0;
      case TideSystem::MeanTide: break;
      }
      throw std::invalid_argument(NO_C20_MT);
    }

    // Offset of the geoid height in system s relative to zero-tide [m].
    inline double n_offset(TideSystem s, double lat_deg, double k)
    {
      switch (s)
      {
      case TideSystem::ZeroTide: return 0.0;
      case TideSystem::TideFree: return -k * direct_geoid_term(lat_deg);
      case TideSystem::MeanTide: return direct_geoid_term(lat_deg);
      }
      return 0.0; // unreachable: the enumeration is closed
    }

    // Vertical crust displacement h_T(phi) [m]; coefficients in mm, converted with 1e-3.
    inline double crust_shift(double lat_deg)
    {
      double s2, s4;
      sin2_sin4(lat_deg, s2, s4);
      return (60.34 - 179.01 * s2 - 1.82 * s4) * 1e-3;
    }

    // Mean-minus-zero geoid separation H_T(phi) [m]; coefficients in mm, converted with 1e-3.
    inline double geoid_separation(double lat_deg)
    {
      double s2, s4;
      sin2_sin4(lat_deg, s2, s4);
      return (99.40 - 295.41 * s2 - 0.42 * s4) * 1e-3;
    }

    // Direct tidal gravity term g_T(phi) [mGal]; coefficients in microGal, converted with 1e-3.
    inline double gravity_tide(double lat_deg)
    {
      double s2, s4;
      sin2_sin4(lat_deg, s2, s4);
      return (-30.49 + 90.95 * s2 + 0.31 * s4) * 1e-3;
    }

    // Gravity offset of each system relative to tide-free, direct term only [mGal].
    inline double gravity_offset(TideSystem s, double lat_deg)
    {
      return s == TideSystem::MeanTide ? gravity_tide(lat_deg) : 0.0;
    }

    // Ellipsoidal height offset of each system relative to tide-free [m].
    inline double ellipsoidal_offset(TideSystem s, double lat_deg)
    {
      return s == TideSystem::TideFree ? 0.0 : crust_shift(lat_deg);
    }

    // Applies f(lat[i]) to every element, in parallel.
    template <class F>
    std::vector<double> map_lat(const std::vector<double> &lat, F f)
    {
      const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(lat.size());
      std::vector<double> out(lat.size());
#pragma omp parallel for schedule(static)
      for (std::ptrdiff_t i = 0; i < n; ++i)
        out[i] = f(lat[i]);
      return out;
    }

    // Applies f(a[i], lat[i]) to every element, in parallel. Throws std::invalid_argument, naming
    // `what`, if the sizes differ.
    template <class F>
    std::vector<double> map_pair(const std::vector<double> &a, const std::vector<double> &lat, const char *what, F f)
    {
      if (a.size() != lat.size())
        throw std::invalid_argument(std::string(what) + ": vectors have different sizes");
      const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(lat.size());
      std::vector<double> out(lat.size());
#pragma omp parallel for schedule(static)
      for (std::ptrdiff_t i = 0; i < n; ++i)
        out[i] = f(a[i], lat[i]);
      return out;
    }
  } // namespace

  // Conversion of C̄20 with an explicit source system. Copies the model and changes only the (2,0)
  // coefficient. Operation order: the shift is computed first, so mean-tide inputs throw before the
  // degree check.
  Model convert_c20_from(const Model &model, TideSystem src, TideSystem dst, double k2)
  {
    if (src == dst)
      return model;
    const double shift = c20_shift(src, dst, k2); // throws if either system is mean-tide
    if (model.nmax < 2)
      throw std::invalid_argument("The model has maximum degree < 2: there is no C̄20 to convert");
    Model out = model;
    out.C[idx_tri(2, 0, model.nmax)] = model.c(2, 0) + shift;
    out.tide_system = to_string(dst);
    return out;
  }

  // ---- Scalar functions ----------------------------------------------------------------------

  double c20_shift(TideSystem from, TideSystem to, double k2)
  {
    return c20_offset(to, k2) - c20_offset(from, k2);
  }

  double height_shift(double lat_deg, TideSystem from, TideSystem to, double k)
  {
    return n_offset(to, lat_deg, k) - n_offset(from, lat_deg, k);
  }

  double convert_geoid_height(double N, double lat_deg, TideSystem from, TideSystem to, double k)
  {
    return N + height_shift(lat_deg, from, to, k);
  }

  double crust_height_shift_m(double lat_deg) { return crust_shift(lat_deg); }

  double geoid_tide_separation_m(double lat_deg) { return geoid_separation(lat_deg); }

  double gravity_tide_mgal(double lat_deg) { return gravity_tide(lat_deg); }

  double gravity_tide_shift_mgal(double lat_deg, TideSystem from, TideSystem to)
  {
    return gravity_offset(to, lat_deg) - gravity_offset(from, lat_deg);
  }

  double convert_ellipsoidal_height(double h, double lat_deg, TideSystem from, TideSystem to)
  {
    return h + (ellipsoidal_offset(to, lat_deg) - ellipsoidal_offset(from, lat_deg));
  }

  double ihrf_separation_surface(double zeta_zero_tide, double lat_deg)
  {
    return zeta_zero_tide - crust_shift(lat_deg) + geoid_separation(lat_deg);
  }

  // Source system taken from model.tide_system; an empty or unknown name is rejected by parse_tide.
  Model convert_model_c20(const Model &model, TideSystem to, double k2)
  {
    if (model.tide_system.empty())
      throw std::invalid_argument("Model has no tide system (tide_system is empty): cannot convert C̄20");
    return convert_c20_from(model, parse_tide(model.tide_system), to, k2);
  }

  // ---- Vectorized functions ------------------------------------------------------------------

  std::vector<double> crust_height_shift_m(const std::vector<double> &lat_deg)
  {
    return map_lat(lat_deg, [](double x) { return crust_shift(x); });
  }

  std::vector<double> geoid_tide_separation_m(const std::vector<double> &lat_deg)
  {
    return map_lat(lat_deg, [](double x) { return geoid_separation(x); });
  }

  std::vector<double> gravity_tide_mgal(const std::vector<double> &lat_deg)
  {
    return map_lat(lat_deg, [](double x) { return gravity_tide(x); });
  }

  std::vector<double> height_shift(const std::vector<double> &lat_deg, TideSystem from, TideSystem to, double k)
  {
    return map_lat(lat_deg, [=](double x) { return height_shift(x, from, to, k); });
  }

  std::vector<double> convert_geoid_height(const std::vector<double> &N, const std::vector<double> &lat_deg,
                                           TideSystem from, TideSystem to, double k)
  {
    return map_pair(N, lat_deg, "convert_geoid_height",
                    [=](double n, double x) { return n + height_shift(x, from, to, k); });
  }

  std::vector<double> gravity_tide_shift_mgal(const std::vector<double> &lat_deg, TideSystem from, TideSystem to)
  {
    return map_lat(lat_deg, [=](double x) { return gravity_tide_shift_mgal(x, from, to); });
  }

  std::vector<double> convert_ellipsoidal_height(const std::vector<double> &h, const std::vector<double> &lat_deg,
                                                 TideSystem from, TideSystem to)
  {
    return map_pair(h, lat_deg, "convert_ellipsoidal_height",
                    [=](double hh, double x) { return convert_ellipsoidal_height(hh, x, from, to); });
  }

  std::vector<double> ihrf_separation_surface(const std::vector<double> &zeta_zero_tide,
                                              const std::vector<double> &lat_deg)
  {
    return map_pair(zeta_zero_tide, lat_deg, "ihrf_separation_surface",
                    [](double z, double x) { return ihrf_separation_surface(z, x); });
  }
}
