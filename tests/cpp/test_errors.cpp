// test_errors.cpp: commission and omission errors (errors.hpp).
// Omission model: Tscherning and Rapp (1974); functional weights: Bucha and Janák (2019).
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <string>
#include "errors.hpp"
#include "functionals.hpp"
#include "legendre.hpp"
#include "shs.hpp"

using namespace gsh;

namespace
{
  constexpr double kGM = 3.986004415e14, kR = 6378136.3;

  // Sigma model with sigma_C = sigma_S = amp / n^1.5 for degrees in [n0, n1] (S sigma of m = 0 is 0).
  Model make_sigma(int n0, int n1, int nmax, double amp)
  {
    Model s = Model::zeros(kGM, kR, nmax);
    for (int n = n0; n <= n1; ++n)
      for (int m = 0; m <= n; ++m)
      {
        s.C[idx_tri(n, m, nmax)] = amp / std::pow(n, 1.5);
        s.S[idx_tri(n, m, nmax)] = m == 0 ? 0.0 : amp * (1.0 + 0.3 * (m % 3)) / std::pow(n, 1.5);
      }
    return s;
  }
}

TEST(Errors, SingleCoefficientMatchesHandFormula)
{
  const int nmax = 8;
  Model s = Model::zeros(kGM, kR, nmax);
  s.C[idx_tri(5, 3, nmax)] = 1e-9;
  s.S[idx_tri(5, 3, nmax)] = 2e-9;
  const double lat = 0.4, lon = 1.1, r = kR + 2.0e5;
  const Points pts = Points::from_spherical({lat * 180 / PI}, {lon * 180 / PI}, {r});
  const auto c = legendre::make_coefs(nmax);
  std::vector<double> P(tri_size(nmax));
  legendre::sfc(std::sin(lat), std::cos(lat), c, P);
  const double p = P[idx_tri(5, 3, nmax)];
  const double sc = std::sqrt(1e-18 * std::pow(std::cos(3 * lon), 2) + 4e-18 * std::pow(std::sin(3 * lon), 2));
  const double q5 = std::pow(kR / r, 5);

  EXPECT_NEAR(commission_error(s, pts, ErrorFunctional::T)[0], kGM / r * q5 * std::fabs(p) * sc, 1e-12 * kGM / r * q5 * std::fabs(p) * sc + 1e-30);
  EXPECT_NEAR(commission_error(s, pts, ErrorFunctional::GravityAnomaly)[0], kGM / (r * r) * 4 * q5 * std::fabs(p) * sc, 1e-18);
  EXPECT_NEAR(commission_error(s, pts, ErrorFunctional::GravityDisturbance)[0], kGM / (r * r) * 6 * q5 * std::fabs(p) * sc, 1e-18);
}

TEST(Errors, CommissionMatchesMonteCarloSynthesis)
{
  // Random coefficients with the given standard deviations, synthesised with compute_functionals (nmin = 21 keeps
  // the ellipsoidal zonal terms of the disturbing model out of the picture). The sample standard deviation of
  // the functional must match the propagated commission error.
  const int nmax = 40, n0 = 21;
  const Model sigma = make_sigma(n0, nmax, nmax, 1e-8);
  const Points pts = Points::from_spherical({10.0, 60.0, -35.0}, {20.0, -100.0, 170.0}, {6.9e6, 6.9e6, 6.9e6});
  ErrorOptions eo;
  eo.nmin = n0;
  const std::vector<std::pair<ErrorFunctional, Functional>> cases = {
      {ErrorFunctional::T, Functional::T}, {ErrorFunctional::ZetaEll, Functional::zeta_ell}, {ErrorFunctional::GravityAnomaly, Functional::Dg_sa}};

  const int K = 3000;
  std::mt19937 rng(12345);
  std::normal_distribution<double> gauss(0.0, 1.0);
  std::vector<std::vector<std::vector<double>>> draws(cases.size(), std::vector<std::vector<double>>(pts.size()));
  FunctionalRequest req;
  req.nmin = n0;
  for (auto &c : cases)
    req.outputs.push_back(c.second);
  for (int k = 0; k < K; ++k)
  {
    Model m = Model::zeros(kGM, kR, nmax);
    for (std::size_t i = 0; i < m.C.size(); ++i)
    {
      m.C[i] = sigma.C[i] * gauss(rng);
      m.S[i] = sigma.S[i] * gauss(rng);
    }
    const FunctionalResult res = compute_functionals(m, pts, req);
    for (std::size_t c = 0; c < cases.size(); ++c)
      for (std::int64_t i = 0; i < pts.size(); ++i)
        draws[c][i].push_back(res.values.at(cases[c].second)[i]);
  }
  for (std::size_t c = 0; c < cases.size(); ++c)
  {
    const auto ce = commission_error(sigma, pts, cases[c].first, eo);
    for (std::int64_t i = 0; i < pts.size(); ++i)
    {
      const auto &v = draws[c][i];
      double mean = 0;
      for (double x : v)
        mean += x;
      mean /= K;
      double var = 0;
      for (double x : v)
        var += (x - mean) * (x - mean);
      const double sd = std::sqrt(var / (K - 1));
      EXPECT_NEAR(sd / ce[i], 1.0, 0.06) << error_functional_name(cases[c].first) << " point " << i;
    }
  }
}

TEST(Errors, GlobalCommissionMatchesSpectralFormula)
{
  const int nmax = 30;
  const Model s = make_sigma(2, nmax, nmax, 1e-9);
  const double a = kR, r = kR + 1.0e5;
  // CE^2(T) = (GM/a)^2 sum_n (a/r)^(2n+2) sum_m (sC^2 + sS^2)   (formula supplied by the project owner)
  double sum = 0.0;
  for (int n = 2; n <= nmax; ++n)
  {
    double s2 = 0.0;
    for (int m = 0; m <= n; ++m)
      s2 += s.c(n, m) * s.c(n, m) + s.s(n, m) * s.s(n, m);
    sum += std::pow(a / r, 2 * n + 2) * s2;
  }
  const double ce_t = kGM / a * std::sqrt(sum);
  EXPECT_NEAR(commission_global(s, ErrorFunctional::T, r, 0.0, 2, nmax), ce_t, 1e-12 * ce_t);

  // CE(N) = CE(T)/γ ~ R sqrt(sum sigma^2) at r = R (the equality uses γ ~ GM/R^2).
  double s2all = 0.0;
  for (int n = 2; n <= nmax; ++n)
    for (int m = 0; m <= n; ++m)
      s2all += s.c(n, m) * s.c(n, m) + s.s(n, m) * s.s(n, m);
  const double gamma = Ellipsoid::GRS80().normal_gravity(0.3, kR).gamma;
  EXPECT_NEAR(commission_global(s, ErrorFunctional::ZetaEll, kR, gamma, 2, nmax) / (kR * std::sqrt(s2all)), 1.0, 0.01);
}

TEST(Errors, PointCommissionAveragesToGlobalValue)
{
  // The area mean of the squared point error equals the global value (Parseval on the sphere).
  const int nmax = 24;
  const Model s = make_sigma(2, nmax, nmax, 1e-9);
  const double r = kR + 3.0e5;
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> u(-1.0, 1.0), l(-180.0, 180.0);
  std::vector<double> lat, lon, rr;
  for (int i = 0; i < 6000; ++i)
  {
    lat.push_back(std::asin(u(rng)) * 180.0 / PI);
    lon.push_back(l(rng));
    rr.push_back(r);
  }
  const auto ce = commission_error(s, Points::from_spherical(lat, lon, rr), ErrorFunctional::T);
  double ms = 0;
  for (double x : ce)
    ms += x * x;
  const double rms = std::sqrt(ms / ce.size());
  EXPECT_NEAR(rms / commission_global(s, ErrorFunctional::T, r, 0.0, 0, nmax), 1.0, 0.05);
}

TEST(Errors, OmissionTscherningRappConsistency)
{
  const double gamma = 9.81;
  TscherningRapp p; // defaults: s = 0.999617 (R_B = 6369.78 km), R_E = 6371 km
  // Direct formula versus the generic degree-variance route (model sphere R = R_E).
  const double GM = kGM, R = p.R_E;
  const int nmax = 360, ntail = 10800;
  const std::vector<double> sig2 = tscherning_rapp_variances(GM, R, 3, ntail, p);
  const double generic = omission_error(sig2, nmax, ErrorFunctional::ZetaEll, GM, R, R, gamma);
  const double direct = omission_tscherning_rapp_undulation(nmax, ntail, gamma, p);
  EXPECT_NEAR(generic, direct, 1e-9 * direct);

  // Anomaly omission at r = R equals the plain sum of c_n (units check mGal^2 -> (m/s^2)^2).
  const double dg = omission_error(sig2, nmax, ErrorFunctional::GravityAnomaly, GM, R, R, gamma);
  double sum = 0;
  for (int n = nmax + 1; n <= ntail; ++n)
    sum += p.A * (n - 1.0) / ((n - 2.0) * (n + p.B)) * std::pow(p.s(), n + 2.0);
  EXPECT_NEAR(dg, std::sqrt(sum) * 1e-5, 1e-9 * dg);

  // The omission error decreases with the truncation degree, and Kaula's rule is the same order of magnitude.
  EXPECT_GT(omission_tscherning_rapp_undulation(180, ntail, gamma, p), direct);
  EXPECT_GT(direct, omission_tscherning_rapp_undulation(720, ntail, gamma, p));
  EXPECT_NEAR(omission_kaula_undulation(6371000.0, 360), 6371000.0 * 1e-5 / 361.0, 1e-15);
  const double kaula = omission_kaula_undulation(6371000.0, nmax);
  EXPECT_LT(direct / kaula, 5.0);
  EXPECT_GT(direct / kaula, 0.2);
}

TEST(Errors, OmissionEqualsMeanSquareOfHighDegreeSignal)
{
  // Parseval: the global mean square of the signal of degrees n0..nmax equals sum w_n^2 sigma_n^2 with the
  // degree variances of the very same coefficients. Checked by sampling the sphere at random points.
  const int n0 = 21, nmax = 40;
  std::mt19937 rng(99);
  std::normal_distribution<double> gauss(0.0, 1.0);
  Model m = Model::zeros(kGM, kR, nmax);
  for (int n = n0; n <= nmax; ++n)
    for (int mm = 0; mm <= n; ++mm)
    {
      m.C[idx_tri(n, mm, nmax)] = 1e-8 / n * gauss(rng);
      m.S[idx_tri(n, mm, nmax)] = mm == 0 ? 0.0 : 1e-8 / n * gauss(rng);
    }
  const std::vector<double> sig2 = degree_variances(m);
  const double r = kR; // radius of the model sphere: no attenuation
  const double expected = omission_error(sig2, n0 - 1, ErrorFunctional::T, kGM, kR, r, 0.0);

  std::uniform_real_distribution<double> u(-1.0, 1.0), l(-180.0, 180.0);
  std::vector<double> lat, lon, rr;
  const int N = 30000;
  for (int i = 0; i < N; ++i)
  {
    lat.push_back(std::asin(u(rng)) * 180.0 / PI);
    lon.push_back(l(rng));
    rr.push_back(r);
  }
  const Points pts = Points::from_spherical(lat, lon, rr);
  const RawFields raw = synth_points(prepare(m, n0, nmax), pts, 0);
  double ms = 0;
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const double T = kGM / r * raw.f[0][i];
    ms += T * T;
  }
  EXPECT_NEAR(std::sqrt(ms / N) / expected, 1.0, 0.05);
}

TEST(Errors, NumericHelpers)
{
  EXPECT_NEAR(relative_error(1.0 + 1e-12, 1.0), 1e-12, 1e-15);
  EXPECT_DOUBLE_EQ(relative_error(2e-3, 0.0), 2e-3);
  EXPECT_NEAR(digit_loss_bound(2190, PI / 2), std::log10(2190.0) + 1.0, 1e-12); // theta cot theta = 0
  EXPECT_NEAR(digit_loss_bound(2190, 0.0), std::log10(2190.0) + std::log10(12.0), 1e-12);
  EXPECT_GT(digit_loss_bound(2190, 3.0), digit_loss_bound(2190, PI / 2)); // near the south pole, theta cot theta is large
  EXPECT_THROW(omission_tscherning_rapp_undulation(1, 100, 9.81), std::invalid_argument);
  EXPECT_THROW(tscherning_rapp_variances(kGM, kR, 2, 10), std::invalid_argument);
}

TEST(Errors, Egm2008SigmasIfPresent)
{
  const std::string path = std::string(GSH_REPO_DATA) + "/cache/modelos/EGM2008/EGM2008.gfc";
  FILE *f = std::fopen(path.c_str(), "rb");
  if (!f)
    GTEST_SKIP() << "EGM2008.gfc not found";
  std::fclose(f);
  const Model s = load_gfc_sigma(path, 360);
  EXPECT_NEAR(s.c(2, 0), 0.7481239490e-11, 1e-20);
  EXPECT_NEAR(s.s(2, 1), 0.7348347201e-11, 1e-20);
  const double gamma = 9.81;
  const double ce = commission_global(s, ErrorFunctional::ZetaEll, kR, gamma, 2, 360);
  const double oe = omission_tscherning_rapp_undulation(360, 10800, gamma);
  std::printf("EGM2008 n=360 undulation: commission %.4f m, omission (TR) %.4f m, Kaula %.4f m\n", ce, oe, omission_kaula_undulation(6371000.0, 360));
  EXPECT_GT(ce, 0.0);
  EXPECT_LT(ce, 0.5);
}
