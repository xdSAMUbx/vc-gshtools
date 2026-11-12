// test_sanity_identities.cpp: analytic identities that the synthesis must satisfy.
//   1) Null ellipsoid: a model that contains only the ellipsoid's own even zonal coefficients (k = 1..4) has a
//      disturbing potential T equal to zero. disturbing_model() (ellipsoid.cpp) subtracts exactly those terms.
//   2) Laplace: outside the masses the gravitational tensor in the LNOF frame is trace-free,
//      V_xx + V_yy + V_zz = 0 (the centrifugal term is not part of V, so it does not enter).
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include "ellipsoid.hpp"
#include "functionals.hpp"
#include "model.hpp"
#include "synth_types.hpp"

using namespace gsh;
using F = Functional;
namespace fs = std::filesystem;

namespace
{
  constexpr double kEps = 2.220446049250313e-16; // double-precision unit roundoff times 2

  // WGS84 defining constants. J2 is the value that gives the WGS84 normalized C20 = -484.16685e-6 (J2 = -sqrt(5) C20).
  // The zonal series of Ellipsoid::zonal depends on e2 and J2 only, so this is a consistent parameter set.
  Ellipsoid wgs84()
  {
    Ellipsoid e;
    e.a = 6378137.0;
    const double f = 1.0 / 298.257223563;
    e.e2 = f * (2.0 - f);
    e.GM = 3.986004418e14;
    e.J2 = 1.08262982e-3;
    e.omega = 7.292115e-5;
    return e;
  }

  // Cartesian product of geodetic latitude, longitude and height (degrees, degrees, metres) on the given ellipsoid.
  Points geodetic_grid(const std::vector<double> &lat, const std::vector<double> &lon, const std::vector<double> &h, const Ellipsoid &ell)
  {
    std::vector<double> la, lo, hh;
    for (double p : lat)
      for (double l : lon)
        for (double z : h)
        {
          la.push_back(p);
          lo.push_back(l);
          hh.push_back(z);
        }
    return Points::from_geodetic(la, lo, hh, ell.a, ell.e2);
  }

  // Random fully normalized model with sigma_C = sigma_S = 1e-5 / n (C00 = 0), up to degree nmax.
  Model random_model(int nmax, unsigned seed)
  {
    Model m = Model::zeros(3.986004415e14, 6378136.3, nmax);
    std::mt19937 rng(seed);
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (int n = 1; n <= nmax; ++n)
      for (int k = 0; k <= n; ++k)
      {
        m.C[idx_tri(n, k, nmax)] = 1e-5 / n * gauss(rng);
        m.S[idx_tri(n, k, nmax)] = k == 0 ? 0.0 : 1e-5 / n * gauss(rng);
      }
    return m;
  }

  struct TraceStats
  {
    double maxTrace = 0;    // max |V_xx + V_yy + V_zz| over the points [s^-2]
    double maxDiagonal = 0; // max |V_xx|, |V_yy|, |V_zz| over the points [s^-2]
  };

  TraceStats laplace_trace(const Model &model, const Points &pts)
  {
    FunctionalRequest rq;
    rq.outputs = {F::V_xx, F::V_yy, F::V_zz};
    rq.nmin = 0;
    rq.nmax = model.nmax;
    const auto res = compute_functionals(model, pts, rq);
    TraceStats s;
    for (std::int64_t i = 0; i < pts.size(); ++i)
    {
      const double xx = res.values.at(F::V_xx)[i], yy = res.values.at(F::V_yy)[i], zz = res.values.at(F::V_zz)[i];
      s.maxTrace = std::max(s.maxTrace, std::abs(xx + yy + zz));
      s.maxDiagonal = std::max({s.maxDiagonal, std::abs(xx), std::abs(yy), std::abs(zz)});
    }
    return s;
  }
} // namespace

// ---------------------------------------------------------------------------------------------------------------
// 1) Null ellipsoid.
// Model: GM = GM_ell, R = a_ell, C00 = 1, C_{2k,0} = C^Ell_{2k,0} (k = 1..4), everything else zero, nmax = 8.
// disturbing_model() gives dC_{2k,0} = C - C^Ell (GM_ell/GM)(a/R)^(2k) = C - C^Ell * 1 * 1 = 0 exactly, including
// the degree-0 term (1 - 1 = 0), because the ratios are exactly 1 in double precision. So T should be exactly 0.
// Units: T in m^2 s^-2; the potential scale is GM/a ~ 6.25e7 m^2 s^-2.
// Rounding argument for the second bound: if the disturbing coefficients were not exactly zero but carried a
// relative error eps in each term of the potential sum, T would be of order eps * GM/a ~ 1.4e-8 m^2 s^-2, which is
// far above 1e-15. The absolute 1e-15 check is therefore only met when that cancellation is exact.
// ---------------------------------------------------------------------------------------------------------------
TEST(SanityNullEllipsoid, DisturbingPotentialVanishes)
{
  const int nmax = 8;
  const std::vector<std::pair<std::string, Ellipsoid>> ells = {{"GRS80", Ellipsoid::GRS80()}, {"WGS84", wgs84()}};
  for (const auto &[name, ell] : ells)
  {
    Model m = Model::zeros(ell.GM, ell.a, nmax);
    m.C[idx_tri(0, 0, nmax)] = 1.0;
    for (int k = 1; k <= 4; ++k)
      m.C[idx_tri(2 * k, 0, nmax)] = ell.zonal(k);

    // Latitudes avoid the exact poles. Geodetic points are converted with the ellipsoid under test.
    const Points p = geodetic_grid({0.0, 15.0, 45.0, 60.0, -70.0, 85.0, -89.0}, {0.0, 30.0, -120.0, 170.0, 10.0, -45.0, 90.0},
                                   {0.0, 5000.0}, ell);
    FunctionalRequest rq;
    rq.outputs = {F::T};
    rq.nmin = 0;
    rq.nmax = nmax;
    rq.ellipsoid = ell;
    const auto res = compute_functionals(m, p, rq);
    double maxT = 0;
    for (double t : res.values.at(F::T))
      maxT = std::max(maxT, std::abs(t));
    const double potential = ell.GM / ell.a;
    std::printf("[METRIC] sanity_null_%s_maxT %.3e m2/s2 (GM/a %.3e, ratio %.3e)\n", name.c_str(), maxT, potential, maxT / potential);

    EXPECT_LT(maxT, 1e-15) << name << ": requested absolute bound";
    // Achievable bound: a few units of roundoff relative to the potential scale GM/a (see comment above).
    EXPECT_LE(maxT, 8.0 * kEps * potential) << name << ": relative roundoff bound";
  }
}

// ---------------------------------------------------------------------------------------------------------------
// 2) Laplace trace of the gravitational tensor (LNOF), exterior space.
// The user tolerance is |trace| < 1e-15 s^-2. The tensor components are of order 1e-6 s^-2 for these models, so
// the tolerance corresponds to a relative level of about 1e-9 with respect to max |V_ii|.
// ---------------------------------------------------------------------------------------------------------------
TEST(SanityLaplaceTrace, RandomModelNmax60)
{
  const Model model = random_model(60, 2024);
  const Points pts = geodetic_grid({0.0, 30.0, -45.0, 70.0, -85.0, 89.0}, {0.0, 120.0, -75.0, 200.0}, {0.0, 5000.0},
                                   Ellipsoid::GRS80());
  const TraceStats s = laplace_trace(model, pts);
  const double rel = s.maxTrace / s.maxDiagonal;
  std::printf("[METRIC] sanity_laplace_random60_maxTrace %.3e s^-2 (max|V_ii| %.3e, relative %.3e)\n", s.maxTrace, s.maxDiagonal, rel);
  EXPECT_LT(s.maxTrace, 1e-15);
  EXPECT_LT(rel, 1e-13); // measured: 3.4e-15 at EGM2008 n = 2190
}

namespace
{
  // Path of an extracted model in the repository data cache, or empty if the file is not present.
  std::string cached_model(const std::string &file)
  {
    const std::string stem = file.substr(0, file.find_last_of('.'));
    const std::string p = std::string(GSH_REPO_DATA) + "/cache/modelos/" + stem + "/" + file;
    return fs::exists(p) ? p : std::string{};
  }

  void check_real_model(int nmax)
  {
    const std::string path = cached_model("EGM2008.gfc");
    if (path.empty())
      GTEST_SKIP() << "EGM2008.gfc is not in data/cache/modelos/EGM2008";
    const Model model = load_gfc(path, nmax);
    const Ellipsoid ell = Ellipsoid::GRS80();
    const Points pts = geodetic_grid({0.0, 38.2, -61.0, 80.0, -89.0}, {254.4, 0.0, 120.0}, {0.0, 5000.0}, ell);
    const TraceStats s = laplace_trace(model, pts);
    const double rel = s.maxTrace / s.maxDiagonal;
    std::printf("[METRIC] sanity_laplace_EGM2008_%d_maxTrace %.3e s^-2 (max|V_ii| %.3e, relative %.3e)\n", nmax, s.maxTrace,
                s.maxDiagonal, rel);
    EXPECT_LT(s.maxTrace, 1e-15);
    EXPECT_LT(rel, 1e-13); // measured: 3.4e-15 at EGM2008 n = 2190
  }
} // namespace

TEST(SanityLaplaceTrace, EGM2008Nmax360) { check_real_model(360); }

TEST(SanityLaplaceTrace, EGM2008Nmax2190) { check_real_model(2190); }
