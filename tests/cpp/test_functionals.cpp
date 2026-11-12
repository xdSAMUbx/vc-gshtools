// Acceptance tests from the checks of Bucha and Janák (2019), plus checks of the
// ellipsoid, the height handling and the Height grid reader.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>
#include "functionals.hpp"
#include "height_model.hpp"
#include "kernels.hpp"
using namespace gsh;
using F = Functional;

namespace
{
  Model random_model(int nmax, unsigned seed)
  {
    Model m = Model::zeros(3.986004415e14, 6378136.3, nmax);
    std::mt19937_64 g(seed);
    std::normal_distribution<double> nd;
    m.C[idx_tri(0, 0, nmax)] = 1.0;
    m.C[idx_tri(2, 0, nmax)] = -4.84e-4;
    for (int mm = 0; mm <= nmax; ++mm)
      for (int n = std::max(mm, 2); n <= nmax; ++n)
      {
        const double amp = 1e-6 / (n * n);
        m.C[idx_tri(n, mm, nmax)] += amp * nd(g);
        if (mm > 0)
          m.S[idx_tri(n, mm, nmax)] = amp * nd(g);
      }
    return m;
  }
  Points some_points(int seed = 1)
  {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> lat(-80, 80), lon(-180, 180), h(0, 10000);
    Points p;
    for (int i = 0; i < 25; ++i)
    {
      p.lat.push_back(deg2rad(lat(g)));
      p.lon.push_back(deg2rad(lon(g)));
      p.r.push_back(6378000.0 + h(g));
    }
    return p;
  }
}

TEST(Ellipsoid, Coefficients)
{
  const auto e = Ellipsoid::GRS80();
  EXPECT_NEAR(e.zonal(1), -4.841669e-4, 1e-9);
  EXPECT_NEAR(e.zonal(2), 7.903041e-7, 1e-11);
}

// Acceptance check 4 of Bucha and Janák (2019): normal gravity against the closed Somigliana formula (GRS80).
TEST(Ellipsoid, NormalGravityAgainstSomigliana)
{
  const auto e = Ellipsoid::GRS80();
  for (double phi_g : {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0, 90.0})
  {
    const double lg = deg2rad(phi_g);
    const auto sph = geodetic_to_spherical(phi_g, 0.0, 0.0, e.a, e.e2);
    const double lat = phi_g == 90.0 ? PI / 2 : sph.lat;
    const auto n = e.normal_gravity(lat, phi_g == 90.0 ? e.radius(PI / 2) : sph.r);
    EXPECT_NEAR(n.gamma, e.somigliana(lg), 5e-10) << phi_g;
  }
  EXPECT_NEAR(e.normal_gravity(0.0, e.a).gamma, 9.7803267715, 5e-10);
  EXPECT_NEAR(e.normal_gravity(PI / 2, e.radius(PI / 2)).gamma, 9.8321863685, 5e-10);
}

// Acceptance check 1 of Bucha and Janák (2019): Laplace's equation for V and T, and the raw-sum tensor against the §4.2 oracle.
TEST(Functionals, LaplaceAndOracleTensor)
{
  const int nmax = 150;
  const Model model = random_model(nmax, 2);
  const Points pts = some_points();
  FunctionalRequest rq;
  rq.outputs = {F::V_xx, F::V_xy, F::V_xz, F::V_yy, F::V_yz, F::V_zz, F::T_xx, F::T_yy, F::T_zz};
  const auto res = compute_functionals(model, pts, rq);
  const auto coefs = legendre::make_coefs(nmax);
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const double s = std::abs(res.values.at(F::V_zz)[i]);
    EXPECT_LT(std::abs(res.values.at(F::V_xx)[i] + res.values.at(F::V_yy)[i] + res.values.at(F::V_zz)[i]), 1e-11 * s);
    EXPECT_LT(std::abs(res.values.at(F::T_xx)[i] + res.values.at(F::T_yy)[i] + res.values.at(F::T_zz)[i]), 1e-11 * s);
    const auto r = reference_point(model, coefs, pts.lat[i], pts.lon[i], pts.r[i]);
    EXPECT_NEAR(res.values.at(F::V_xx)[i], r.vxx, 2e-11 * s);
    EXPECT_NEAR(res.values.at(F::V_xy)[i], r.vxy, 2e-11 * s);
    EXPECT_NEAR(res.values.at(F::V_xz)[i], r.vxz, 2e-11 * s);
    EXPECT_NEAR(res.values.at(F::V_yy)[i], r.vyy, 2e-11 * s);
    EXPECT_NEAR(res.values.at(F::V_yz)[i], r.vyz, 2e-11 * s);
    EXPECT_NEAR(res.values.at(F::V_zz)[i], r.vzz, 2e-11 * s);
  }
}

TEST(Functionals, PolarPointsUseNonSingularTensor)
{
  const int nmax = 80;
  const Model model = random_model(nmax, 2);
  Points p;
  for (double lat : {PI / 2, -PI / 2, PI / 2 - 1e-5})
  {
    p.lat.push_back(lat);
    p.lon.push_back(0.7);
    p.r.push_back(6370000.0);
  }
  FunctionalRequest rq;
  rq.outputs = {F::V_xx, F::V_xy, F::V_xz, F::V_yy, F::V_yz, F::V_zz};
  const auto res = compute_functionals(model, p, rq);
  EXPECT_EQ(res.polar_points, 3);
  FunctionalRequest bad;
  bad.outputs = {F::eta};
  EXPECT_THROW(compute_functionals(model, p, bad), std::invalid_argument);
  for (const auto &kv : res.values)
    for (double v : kv.second)
      EXPECT_TRUE(std::isfinite(v)) << functional_name(kv.first);
  for (std::int64_t i = 0; i < p.size(); ++i)
    EXPECT_LT(std::abs(res.values.at(F::V_xx)[i] + res.values.at(F::V_yy)[i] + res.values.at(F::V_zz)[i]),
              1e-10 * std::abs(res.values.at(F::V_zz)[i]));
}

// Acceptance check 5 of Bucha and Janák (2019): the ellipsoid against itself. All disturbing functionals must vanish.
TEST(Functionals, EllipsoidAgainstItself)
{
  const auto e = Ellipsoid::GRS80();
  Model m = Model::zeros(e.GM, e.a, 20);
  for (int k = 0; k <= 10; ++k)
    m.C[idx_tri(2 * k, 0, 20)] = e.zonal(k);
  const Points pts = some_points();
  FunctionalRequest rq;
  rq.outputs = {F::T, F::dg_sa, F::Dg_sa, F::T_rr, F::T_xx, F::T_xy, F::T_xz, F::T_yy, F::T_yz, F::T_zz, F::dg, F::xi, F::eta};
  const auto res = compute_functionals(m, pts, rq);
  for (const auto &kv : res.values)
    for (double v : kv.second)
      EXPECT_NEAR(v, 0.0, 1e-9) << functional_name(kv.first);
}

// Spec §9, items 6 and 7: a single coefficient, and the relation between functionals.
TEST(Functionals, SingleCoefficient)
{
  const auto e = Ellipsoid::GRS80();
  const double x = 1e-5;
  Ellipsoid zero = e;
  // With a very small GM_ell the normal field vanishes, so ΔC̄_20 = C̄_20.
  zero.GM = 1e-300;
  Model m = Model::zeros(3.986e14, 6378136.3, 2);
  m.C[idx_tri(2, 0, 2)] = x;
  const Points pts = some_points();
  FunctionalRequest rq;
  rq.ellipsoid = zero;
  rq.outputs = {F::T, F::dg_sa, F::Dg_sa};
  const auto res = compute_functionals(m, pts, rq);
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const double r = pts.r[i], q = m.R / r, t = std::sin(pts.lat[i]);
    const double P20 = std::sqrt(5.0) / 2 * (3 * t * t - 1);
    const double T = m.GM / r * q * q * x * P20;
    EXPECT_NEAR(res.values.at(F::T)[i], T, 1e-10 * std::abs(T) + 1e-14);
    EXPECT_NEAR(res.values.at(F::Dg_sa)[i], T / r, 1e-10 * std::abs(T / r) + 1e-18);
    EXPECT_NEAR(res.values.at(F::Dg_sa)[i], res.values.at(F::dg_sa)[i] - 2 * T / r, 1e-10 * std::abs(T / r) + 1e-18);
  }
}

// Acceptance check 8 of Bucha and Janák (2019): g_x, g_y, g_z against the finite-difference gradient of W.
TEST(Functionals, GravityVectorAgainstFiniteDifferences)
{
  const int nmax = 40;
  const Model model = random_model(nmax, 5);
  Points base = some_points(3);
  FunctionalRequest rq;
  rq.outputs = {F::g_x, F::g_y, F::g_z};
  const auto res = compute_functionals(model, base, rq);

  FunctionalRequest rw;
  rw.outputs = {F::W};
  const double dl = 1e-6, dr = 5.0;
  const auto shift = [&](double dlat, double dlon, double drad)
  {
    Points p = base;
    for (std::int64_t i = 0; i < p.size(); ++i)
    {
      p.lat[i] += dlat;
      p.lon[i] += dlon;
      p.r[i] += drad;
    }
    return compute_functionals(model, p, rw).values.at(F::W);
  };
  const auto wn = shift(dl, 0, 0), ws = shift(-dl, 0, 0), we = shift(0, dl, 0), ww = shift(0, -dl, 0), wu = shift(0, 0, dr), wd = shift(0, 0, -dr);
  for (std::int64_t i = 0; i < base.size(); ++i)
  {
    const double r = base.r[i], c = std::cos(base.lat[i]);
    EXPECT_NEAR(res.values.at(F::g_x)[i], (wn[i] - ws[i]) / (2 * dl * r), 3e-5);
    EXPECT_NEAR(res.values.at(F::g_y)[i], -(we[i] - ww[i]) / (2 * dl * r * c), 3e-5);
    EXPECT_NEAR(res.values.at(F::g_z)[i], (wu[i] - wd[i]) / (2 * dr), 3e-5);
  }
}

TEST(Functionals, AllRunAndHeightHandling)
{
  const int nmax = 60;
  const Model model = random_model(nmax, 9);
  const Points pts = some_points();

  // Synthetic height model: a constant H0 = C00.
  const double H0 = 1200.0;
  Model hm = Model::zeros(1.0, 1.0, 2);
  hm.C[0] = H0;
  FunctionalRequest rq;
  rq.outputs = all_functionals();
  rq.height = height_from_model(hm);
  const auto res = compute_functionals(model, pts, rq);
  EXPECT_EQ(res.values.size(), static_cast<std::size_t>(Functional::Count));
  EXPECT_EQ(res.points_without_height, 0);
  for (const auto &kv : res.values)
    for (double v : kv.second)
      EXPECT_TRUE(std::isfinite(v)) << functional_name(kv.first);

  // N from its direct definition (PDF 76) at each point.
  FunctionalRequest rt;
  rt.outputs = {F::T};
  const auto &e = rq.ellipsoid;
  Points foot = pts;
  for (std::int64_t i = 0; i < pts.size(); ++i)
    foot.r[i] = e.radius(pts.lat[i]);
  const auto Tf = compute_functionals(model, foot, rt).values.at(F::T);
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const double g = e.normal_gravity(pts.lat[i], foot.r[i]).gamma;
    EXPECT_NEAR(res.values.at(F::N)[i], (Tf[i] - 2 * PI * rq.G * rq.density * H0 * H0) / g, 1e-8);
  }
  rq.height = nullptr;
  EXPECT_THROW(compute_functionals(model, pts, rq), std::invalid_argument);
}

TEST(Height, GridParseInterpolateAndCoverage)
{
  namespace fs = std::filesystem;
  const auto path = fs::temp_directory_path() / "gsh_grid.grd";
  {
    std::ofstream o(path);
    o << "   35.0 36.0 250.0 252.0 0.5 1.0\n\n";
    o << "100 110 120\n200 210 220\n300 310 320\n";
  }
  const auto g = load_gravsoft_grid(path.string());
  EXPECT_EQ(g.nlat, 3);
  EXPECT_EQ(g.nlon, 3);
  EXPECT_NEAR(g.interpolate(36.0, 250.0), 100.0, 1e-5); // first row = north
  EXPECT_NEAR(g.interpolate(35.0, 252.0), 320.0, 1e-5);
  EXPECT_NEAR(g.interpolate(35.5, -109.0), 210.0, 1e-5); // negative longitude
  EXPECT_NEAR(g.interpolate(35.75, 250.5), 155.0, 1e-5);
  EXPECT_TRUE(std::isnan(g.interpolate(34.9, 251.0)));
  EXPECT_TRUE(std::isnan(g.interpolate(35.5, 253.0)));
  const auto m = g.block_mean(1.0);
  EXPECT_EQ(m.nlat, 2);
  fs::remove(path);
}

TEST(Height, ColoradoDemIfPresent)
{
  const std::string path = std::string(GSH_REPO_DATA) + "/cache/tablas/Colorado_dem/Colorado_dem.grd";
  if (!std::filesystem::exists(path))
    GTEST_SKIP();
  const auto g = load_gravsoft_grid(path);
  EXPECT_EQ(g.nlat, 6001);
  EXPECT_EQ(g.nlon, 9601);
  EXPECT_NEAR(g.dlat(), 5.0 / 6000.0, 1e-12);
  float mx = 0, mn = 1e9;
  std::size_t at = 0;
  for (std::size_t i = 0; i < g.h.size(); ++i)
  {
    if (g.h[i] > mx)
    {
      mx = g.h[i];
      at = i;
    }
    mn = std::min(mn, g.h[i]);
  }
  EXPECT_EQ(mx, 4385.0f);
  EXPECT_EQ(mn, 932.0f);
  EXPECT_EQ(at, 1059ull * 9601 + 4266); // Mount Elbert
}
