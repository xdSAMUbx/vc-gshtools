// Grid synthesis must equal point-by-point synthesis at the same locations (grid_to_points).
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include "grid.hpp"
#include "kernels.hpp"
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
  void compare(const RawFields &grid, const RawFields &points, double tol, const char *label)
  {
    for (std::size_t k = 0; k < grid.f.size(); ++k)
    {
      double d = 0, s = 0;
      for (std::size_t i = 0; i < grid.f[k].size(); ++i)
      {
        d = std::max(d, std::abs(grid.f[k][i] - points.f[k][i]));
        s = std::max(s, std::abs(points.f[k][i]));
      }
      EXPECT_LT(d, tol * std::max(s, 1e-300)) << label << " field " << k;
    }
  }
}

class GridVsPoints : public ::testing::TestWithParam<int>
{
};

TEST_P(GridVsPoints, GaussLegendreAllGrads)
{
  const int nmax = GetParam();
  const Model model = random_model(nmax, 1);
  const auto prep = prepare(model);
  const Grid g = Grid::gauss_legendre(nmax, model.R * 1.001);
  const Points pts = grid_to_points(g);
  SynthOptions o;
  o.polar_optimization = false;
  for (int grad = 0; grad <= 2; ++grad)
    compare(synth_grid(prep, g, grad, o), synth_points(prep, pts, grad, o), 2e-11, "GL");
}
INSTANTIATE_TEST_SUITE_P(Degrees, GridVsPoints, ::testing::Values(1, 2, 3, 16, 61));

TEST(Grid, DriscollHealyWithPoleRow)
{
  const int nmax = 24;
  const Model model = random_model(nmax, 2);
  const auto prep = prepare(model);
  SynthOptions o;
  o.polar_optimization = false;
  const auto d1 = Grid::driscoll_healy1(nmax, model.R);
  compare(synth_grid(prep, d1, 0, o), synth_points(prep, grid_to_points(d1), 0, o), 1e-11, "DH1");
  const auto d2 = Grid::driscoll_healy2(nmax, model.R);
  compare(synth_grid(prep, d2, 0, o), synth_points(prep, grid_to_points(d2), 0, o), 1e-11, "DH2");

  // With derivatives, the pole row comes from the oracle and must be finite.
  const auto raw = synth_grid(prep, d1, 2, o);
  for (const auto &v : raw.f)
    for (double x : v)
      ASSERT_TRUE(std::isfinite(x));
  const auto coefs = legendre::make_coefs(nmax);
  const auto ref = reference_point(model, coefs, PI / 2, d1.dlon * 3, model.R);
  for (int k = 0; k < F_COUNT; ++k)
    EXPECT_NEAR(raw.f[k][3], ref.raw[k], 1e-10 * (std::abs(ref.raw[k]) + 1.0)) << k;
}

TEST(Grid, RegionalUsesDirectSynthesisAndAsymmetricLatitudes)
{
  const int nmax = 40;
  const Model model = random_model(nmax, 3);
  const auto prep = prepare(model);
  SynthOptions o;
  o.polar_optimization = false;
  const Grid g = Grid::regular(45.0, 30.0, 7, 250.0, 260.0, 11, model.R + 5000);
  EXPECT_FALSE(g.full_circle());
  compare(synth_grid(prep, g, 2, o), synth_points(prep, grid_to_points(g), 2, o), 1e-11, "regional");
}

TEST(Grid, FewLongitudesAliasing)
{
  // nlon small compared with nmax: the frequency folding must reproduce the direct evaluation.
  const int nmax = 30;
  const Model model = random_model(nmax, 4);
  const auto prep = prepare(model);
  Grid g = Grid::gauss_legendre(nmax, model.R);
  g.nlon = 16;
  g.dlon = 2 * PI / 16;
  g.lon0 = 0.3;
  SynthOptions o;
  o.polar_optimization = false;
  ASSERT_TRUE(g.full_circle());
  compare(synth_grid(prep, g, 1, o), synth_points(prep, grid_to_points(g), 1, o), 1e-11, "alias");
}

TEST(Grid, RadiusPerLatitudeBreaksSymmetrySafely)
{
  const int nmax = 30;
  const Model model = random_model(nmax, 5);
  const auto prep = prepare(model);
  Grid g = Grid::gauss_legendre(nmax, model.R);
  for (std::int64_t i = 0; i < g.nlat(); ++i)
    g.r[i] = model.R + 1000.0 * i;
  SynthOptions o;
  o.polar_optimization = false;
  compare(synth_grid(prep, g, 1, o), synth_points(prep, grid_to_points(g), 1, o), 1e-11, "radios");
}
