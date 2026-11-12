// Kernels against the sequential oracle (reference_point), and the local tensor against its non-singular form
// and against finite differences of the potential.
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include "kernels.hpp"
#include "shs.hpp"
using namespace gsh;
using namespace gsh::kernels;

namespace
{
  Model random_model(int nmax, unsigned seed)
  {
    Model m = Model::zeros(3.986004415e14, 6378136.3, nmax);
    std::mt19937_64 g(seed);
    std::normal_distribution<double> nd;
    m.C[idx_tri(0, 0, nmax)] = 1.0;
    for (int mm = 0; mm <= nmax; ++mm)
      for (int n = std::max(mm, 2); n <= nmax; ++n)
      {
        const double amp = 1e-5 / (n * n);
        m.C[idx_tri(n, mm, nmax)] = amp * nd(g);
        if (mm > 0)
          m.S[idx_tri(n, mm, nmax)] = amp * nd(g);
      }
    return m;
  }
  Points random_points(std::int64_t np, unsigned seed, double latmax = 85.0)
  {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> lat(-latmax, latmax), lon(-180, 180), h(-500, 20000);
    Points p;
    for (std::int64_t i = 0; i < np; ++i)
    {
      p.lat.push_back(deg2rad(lat(g)));
      p.lon.push_back(deg2rad(lon(g)));
      p.r.push_back(6371000.0 + h(g));
    }
    return p;
  }
  double max_rel(const AlignedVec<double> &a, const std::vector<double> &b)
  {
    double scale = 0, d = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
      scale = std::max(scale, std::abs(b[i]));
      d = std::max(d, std::abs(a[i] - b[i]));
    }
    return d / std::max(scale, 1e-300);
  }
}

class KernelVsOracle : public ::testing::TestWithParam<std::tuple<Isa, int>>
{
};

TEST_P(KernelVsOracle, AllFields)
{
  const auto [isa, nmax] = GetParam();
  // A forced ISA that the CPU lacks makes the synthesis throw (e.g. CI runners without AVX-512): skip, do not fail.
  try
  {
    (void)active_isa_name(isa);
  }
  catch (const std::runtime_error &e)
  {
    GTEST_SKIP() << e.what();
  }
  const Model model = random_model(nmax, 7);
  const auto prep = prepare(model);
  const auto coefs = legendre::make_coefs(nmax);
  // Point count that is not a multiple of any tile width, so the padding path is exercised.
  const Points pts = random_points(37, 11);
  SynthOptions o;
  o.isa = isa;
  o.polar_optimization = false;
  std::vector<std::vector<double>> ref(F_COUNT, std::vector<double>(pts.size()));
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const auto r = reference_point(model, coefs, pts.lat[i], pts.lon[i], pts.r[i]);
    for (int k = 0; k < F_COUNT; ++k)
      ref[k][i] = r.raw[k];
  }
  for (int grad = 0; grad <= 2; ++grad)
  {
    const auto raw = synth_points(prep, pts, grad, o);
    for (int k = 0; k < fields_for_grad(grad); ++k)
      EXPECT_LT(max_rel(raw.f[k], ref[k]), 3e-12) << "grad " << grad << " field " << k;
  }
}

INSTANTIATE_TEST_SUITE_P(Isas, KernelVsOracle,
                         ::testing::Combine(::testing::Values(Isa::Scalar, Isa::Avx2, Isa::Avx512),
                                            ::testing::Values(2, 60, 360)));

TEST(Synthesis, HighDegreeMatchesOracle)
{
  // Ultra-high degree with exponent rescaling and difficult colatitudes: forbidden zone and MSFC.
  const int nmax = 2200;
  const Model model = random_model(nmax, 3);
  const auto prep = prepare(model);
  const auto coefs = legendre::make_coefs(nmax);
  Points pts;
  for (double lat : {70.1, 38.3, -61.0, 80.0, 5.0, 89.9})
  {
    pts.lat.push_back(deg2rad(lat));
    pts.lon.push_back(deg2rad(lat * 1.7));
    pts.r.push_back(6371000.0 + 9000.0);
  }
  std::vector<std::vector<double>> ref(F_COUNT, std::vector<double>(pts.size()));
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const auto r = reference_point(model, coefs, pts.lat[i], pts.lon[i], pts.r[i]);
    for (int k = 0; k < F_COUNT; ++k)
      ref[k][i] = r.raw[k];
  }
  for (bool skip : {false, true})
  {
    SynthOptions o;
    o.polar_optimization = false;
    o.skip_forbidden = skip;
    const auto raw = synth_points(prep, pts, 2, o);
    for (int k = 0; k < F_COUNT; ++k)
      EXPECT_LT(max_rel(raw.f[k], ref[k]), 1e-10) << "skip " << skip << " field " << k;
  }
}

TEST(Synthesis, OptionsAgree)
{
  const int nmax = 500;
  const Model model = random_model(nmax, 5);
  const auto prep = prepare(model);
  const Points pts = random_points(1001, 9, 80.0);
  SynthOptions base;
  base.polar_optimization = false;
  base.skip_forbidden = false;
  base.rotate_trig = false;
  base.sort_latitude = false;
  base.threads = 1;
  base.tile_vectors = 1;
  base.parallel_orders = 0;
  const auto ref = synth_points(prep, pts, 1, base);
  const auto compare = [&](SynthOptions o, const char *label)
  {
    const auto raw = synth_points(prep, pts, 1, o);
    for (int k = 0; k < 4; ++k)
    {
      std::vector<double> b(ref.f[k].begin(), ref.f[k].end());
      EXPECT_LT(max_rel(raw.f[k], b), 1e-12) << label << " field " << k;
    }
  };
  for (int b : {1, 2, 4, 8, 16})
  {
    SynthOptions o = base;
    o.tile_vectors = b;
    o.threads = 0;
    compare(o, "tile");
  }
  SynthOptions o = base;
  o.skip_forbidden = true;
  compare(o, "skip");
  o = base;
  o.rotate_trig = true;
  compare(o, "rot");
  o = base;
  o.sort_latitude = true;
  o.threads = 0;
  compare(o, "sort");
  o = base;
  o.parallel_orders = 1;
  o.threads = 4;
  compare(o, "by order");
}

TEST(Synthesis, TileResultIndependentOfThreads)
{
  const int nmax = 200;
  const Model model = random_model(nmax, 5);
  const auto prep = prepare(model);
  const Points pts = random_points(500, 21);
  SynthOptions a;
  a.threads = 1;
  a.parallel_orders = 0;
  SynthOptions b = a;
  b.threads = 8;
  const auto ra = synth_points(prep, pts, 2, a), rb = synth_points(prep, pts, 2, b);
  for (int k = 0; k < F_COUNT; ++k)
    for (std::int64_t i = 0; i < pts.size(); ++i)
      ASSERT_EQ(ra.f[k][i], rb.f[k][i]);
}

TEST(Synthesis, Nmin)
{
  const int nmax = 50;
  Model model = random_model(nmax, 13);
  const auto prep = prepare(model, 10, nmax);
  Model cut = model;
  for (int m = 0; m <= nmax; ++m)
    for (int n = m; n < 10; ++n)
    {
      cut.C[idx_tri(n, m, nmax)] = 0;
      cut.S[idx_tri(n, m, nmax)] = 0;
    }
  const auto prep2 = prepare(cut);
  const Points pts = random_points(40, 3);
  SynthOptions o;
  o.parallel_orders = 0; // parallelism over m changes the order of the sums between runs
  const auto a = synth_points(prep, pts, 0, o), b = synth_points(prep2, pts, 0, o);
  for (std::int64_t i = 0; i < pts.size(); ++i)
    ASSERT_EQ(a.f[0][i], b.f[0][i]);
}

TEST(Oracle, LaplaceAndHessianConsistency)
{
  // Zero trace, and the raw-sum form equals Bucha and Janák (2019) away from the poles.
  const int nmax = 120;
  const Model model = random_model(nmax, 17);
  const auto coefs = legendre::make_coefs(nmax);
  for (double latd : {-70.0, -12.0, 3.0, 45.0, 88.0})
  {
    const double lat = deg2rad(latd), lon = 1.1, r = 6378136.3 + 3000.0;
    const auto R = reference_point(model, coefs, lat, lon, r);
    const double scale = std::abs(R.vzz);
    EXPECT_LT(std::abs(R.vxx + R.vyy + R.vzz), 1e-11 * scale);
    const double k = model.GM / (r * r * r), cl = std::cos(lat), tl = std::tan(lat);
    const double *S = R.raw;
    EXPECT_NEAR(R.vxx, k * (S[F_SPP] - S[F_SR]), 1e-11 * scale);
    EXPECT_NEAR(R.vyy, k * (-S[F_SLL] / (cl * cl) - S[F_SR] - tl * S[F_SP]), 2e-10 * scale);
    EXPECT_NEAR(R.vzz, k * S[F_SRR], 1e-11 * scale);
    EXPECT_NEAR(R.vxz, -k * (S[F_SRP] + S[F_SP]), 1e-11 * scale);
    EXPECT_NEAR(R.vxy, -k * (S[F_SPL] + tl * S[F_SL]) / cl, 2e-10 * scale);
    EXPECT_NEAR(R.vyz, k * (S[F_SRL] + S[F_SL]) / cl, 2e-10 * scale);
  }
}

TEST(Oracle, TensorAgainstFiniteDifferences)
{
  // Rotated Cartesian Hessian from finite differences (check 2 of Bucha and Janák (2019)).
  const int nmax = 30;
  const Model model = random_model(nmax, 23);
  const auto coefs = legendre::make_coefs(nmax);
  const double lat0 = 0.6, lon0 = 0.9, r0 = 6500000.0;
  const auto V = [&](double x, double y, double z)
  {
    const double r = std::sqrt(x * x + y * y + z * z);
    const auto rp = reference_point(model, coefs, std::asin(z / r), std::atan2(y, x), r);
    return model.GM / r * rp.raw[F_S0];
  };
  const double cl = std::cos(lat0), sl = std::sin(lat0), co = std::cos(lon0), so = std::sin(lon0);
  const double X[3] = {r0 * cl * co, r0 * cl * so, r0 * sl};
  // Axes of the LNOF frame: x north, y west, z up.
  const double ex[3] = {-sl * co, -sl * so, cl}, ey[3] = {so, -co, 0}, ez[3] = {cl * co, cl * so, sl};
  const double *axes[3] = {ex, ey, ez};
  const double h = 300.0;
  const auto second = [&](const double *a, const double *b)
  {
    const auto f = [&](double s, double t)
    { return V(X[0] + s * a[0] + t * b[0], X[1] + s * a[1] + t * b[1], X[2] + s * a[2] + t * b[2]); };
    return (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h);
  };
  const auto R = reference_point(model, coefs, lat0, lon0, r0);
  const double vals[3][3] = {{R.vxx, R.vxy, R.vxz}, {R.vxy, R.vyy, R.vyz}, {R.vxz, R.vyz, R.vzz}};
  for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j)
    {
      const double fd = second(axes[i], axes[j]);
      EXPECT_NEAR(vals[i][j], fd, 2e-5 * std::abs(R.vzz) + 1e-11) << i << "," << j;
    }
}
