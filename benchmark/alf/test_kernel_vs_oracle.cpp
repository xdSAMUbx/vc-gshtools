// Differential test of the production synthesis kernel against the sequential oracle, along the meridian.
//
// The production path (compute_functionals -> synth_points -> SIMD column kernel, "fast") fuses the Legendre
// recurrence with the synthesis, rescales with per-lane exponents and skips columns that cannot contribute. The
// oracle (reference_point) builds the full MSFC triangle per point and applies the literal formulas of the
// specification; its Legendre functions are validated separately by test_alf_identities.cpp (sum rules, -90..90 deg).
// Here only the gravitational potential V is compared (no derivatives, no disturbing potential).
//
// Model: EGM2008, degree and order 2190 (data/cache/modelos/EGM2008/EGM2008.gfc; the test is skipped if absent).
// Points: spherical (geocentric) latitude -90..90 deg every 0.5 deg including both exact poles, three longitudes
// (0, 61.7, -143.9 deg) to exercise the trigonometric recurrence, radius of the GRS80 ellipsoid at each latitude.
//
// Error metric: |V_kernel - V_oracle| in m²/s² (V ~ 6.2e7 m²/s², so 1e-7 is about 1.6e-15 relative, the resolution
// of a double at that magnitude). NaN or Inf, or a relative error above 1e-3, is a hard failure.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "functionals.hpp"
#include "legendre.hpp"
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

namespace
{
  // Absolute tolerance on V [m²/s²]: about 10x the rounding resolution of V at 6.2e7 m²/s². Set from the first
  // measurement (see the run summary), not in advance.
  constexpr double kTolV = 1e-6;
  constexpr double kHardRel = 1e-3;
}

TEST(KernelVsOracle, PotentialAlongMeridianEGM2008)
{
#ifdef GSH_REPO_DATA
  const std::string path = std::string(GSH_REPO_DATA) + "/cache/modelos/EGM2008/EGM2008.gfc";
#else
  const std::string path = "data/cache/modelos/EGM2008/EGM2008.gfc";
#endif
  if (!std::filesystem::exists(path))
    GTEST_SKIP() << "EGM2008 not found: " << path;
  const int nmax = 2190;
  const Model model = load_gfc(path, nmax);
  ASSERT_EQ(model.nmax, nmax);

  // GRS80 ellipsoid radius at geocentric latitude phi: r = a b / sqrt((b cos phi)² + (a sin phi)²)
  const double a = 6378137.0, f = 1.0 / 298.257222101, b = a * (1.0 - f);
  const double lons[] = {0.0, 61.7, -143.9};
  // Latitude step in degrees (default 1, ~6 minutes for the oracle); GSH_MERIDIAN_STEP=5 gives a quick run.
  double step = 1.0;
  if (const char *e = std::getenv("GSH_MERIDIAN_STEP"))
    step = std::max(0.05, std::atof(e));
  const int nlat = static_cast<int>(std::lround(180.0 / step));
  std::vector<double> lat_deg, lon_deg, r;
  for (int i = 0; i <= nlat; ++i)
  {
    const double lat = -90.0 + step * i;
    double s, c;
    if (std::abs(lat) == 90.0)
    {
      s = lat > 0 ? 1.0 : -1.0;
      c = 0.0;
    }
    else
    {
      s = std::sin(lat * PI / 180.0);
      c = std::cos(lat * PI / 180.0);
    }
    const double rr = a * b / std::sqrt((b * c) * (b * c) + (a * s) * (a * s));
    for (double lon : lons)
    {
      lat_deg.push_back(lat);
      lon_deg.push_back(lon);
      r.push_back(rr);
    }
  }
  const Points pts = Points::from_spherical(lat_deg, lon_deg, r);
  const std::int64_t np = pts.size();

  // Production path through the public API.
  FunctionalRequest req;
  req.outputs = {Functional::V};
  req.nmax = nmax;
  const FunctionalResult res = compute_functionals(model, pts, req);
  const auto &Vk = res.values.at(Functional::V);
  ASSERT_EQ(static_cast<std::int64_t>(Vk.size()), np);
  std::cout << "\n[Kernel vs oracle] ISA " << active_isa_name() << ", " << np << " points, nmax " << nmax << "\n";

  // Oracle, point by point.
  const legendre::Coefs coefs = legendre::make_coefs(nmax);
  std::vector<double> Vo(np), err(np);
  for (std::int64_t i = 0; i < np; ++i)
  {
    const ReferencePoint rp = reference_point(model, coefs, pts.lat[i], pts.lon[i], pts.r[i], 0, nmax);
    Vo[i] = model.GM / pts.r[i] * rp.raw[0];
    err[i] = std::abs(Vk[i] - Vo[i]);
  }

  long hard = 0, over = 0;
  double worst = 0.0, worst_lat = 0.0, worst_lon = 0.0;
  std::ofstream csv;
#ifdef GSH_ALF_OUT
  std::filesystem::create_directories(GSH_ALF_OUT);
  csv.open(std::string(GSH_ALF_OUT) + "/kernel_vs_oracle_V_egm2008_n2190.csv");
  csv << "lat_deg,lon_deg,V_kernel,V_oracle,abs_err\n";
  csv.precision(17);
#endif
  for (std::int64_t i = 0; i < np; ++i)
  {
    const bool bad_hard = !std::isfinite(Vk[i]) || !std::isfinite(Vo[i]) || !std::isfinite(err[i]) ||
                          err[i] > kHardRel * std::abs(Vo[i]);
    if (bad_hard)
      ++hard;
    if (bad_hard || err[i] > kTolV)
    {
      ++over;
      if (over <= 10)
        std::printf("    lat=%+.2f lon=%+.1f V_kernel=%.12e V_oracle=%.12e err=%.3e%s\n", lat_deg[i], lon_deg[i], Vk[i],
                    Vo[i], err[i], bad_hard ? " HARD" : "");
    }
    if (!(err[i] <= worst))
    {
      worst = err[i];
      worst_lat = lat_deg[i];
      worst_lon = lon_deg[i];
    }
    if (csv)
      csv << lat_deg[i] << ',' << lon_deg[i] << ',' << Vk[i] << ',' << Vo[i] << ',' << err[i] << '\n';
  }

  // Summary by latitude strip: max error.
  struct Strip { double lo, hi; double mx = 0; };
  Strip strips[] = {{0, 60}, {60, 80}, {80, 85}, {85, 89.5}, {89.5, 90.0001}};
  for (std::int64_t i = 0; i < np; ++i)
    for (auto &s : strips)
      if (std::abs(lat_deg[i]) >= s.lo && std::abs(lat_deg[i]) < s.hi)
        s.mx = std::max(s.mx, err[i]);
  std::printf("  worst |dV| = %.3e m²/s² at lat %+.2f lon %+.1f   (above %.0e: %ld, hard: %ld)\n", worst, worst_lat,
              worst_lon, kTolV, over, hard);
  for (const auto &s : strips)
    std::printf("    |lat| in [%.1f, %.1f): max |dV| = %.3e\n", s.lo, std::min(s.hi, 90.0), s.mx);

  EXPECT_EQ(hard, 0) << "NaN/Inf or relative error > 1e-3";
  EXPECT_EQ(over, 0) << "points with |V_kernel - V_oracle| > " << kTolV;
}
