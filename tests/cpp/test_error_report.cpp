// test_error_report.cpp: consolidated error report (error_report.hpp).
// Checks the Kaula value at n = 2190 against the quoted 2.9 cm, and internal consistency and formats for the rest.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include "error_report.hpp"
#include "errors.hpp"

using namespace gsh;

namespace
{
  constexpr double kGM = 3.986004415e14, kR = 6378136.3;

  // Sigma model (C and S hold standard deviations) with degrees 2..nmax.
  Model make_sigma(int nmax, double amp)
  {
    Model s = Model::zeros(kGM, kR, nmax);
    for (int n = 2; n <= nmax; ++n)
      for (int m = 0; m <= n; ++m)
      {
        s.C[idx_tri(n, m, nmax)] = amp / std::pow(n, 1.5);
        s.S[idx_tri(n, m, nmax)] = m == 0 ? 0.0 : amp * (1.0 + 0.3 * (m % 3)) / std::pow(n, 1.5);
      }
    s.name = "test model";
    return s;
  }

  std::size_t count_char(const std::string &s, char c)
  {
    return static_cast<std::size_t>(std::count(s.begin(), s.end(), c));
  }
}

TEST(ErrorReport, KaulaOmissionAtN2190IsAboutTwoPointNineCm)
{
  // Quoted by the project owner: R 1e-5 / (2190 + 1) = 2.9 cm.
  const Model s = make_sigma(10, 1e-8);
  const ErrorReport rep = make_error_report(s, 2190, kR);
  EXPECT_NEAR(rep.omission.kaula, 0.029, 0.0005);
  EXPECT_NEAR(rep.omission.kaula, kR * 1e-5 / 2191.0, 1e-15);
  EXPECT_GT(rep.omission.tr, 0.0);
  EXPECT_FALSE(rep.omission.has_tr_model); // run nmax is above the model degree: no model variances
  EXPECT_FALSE(rep.numerical.available);
  EXPECT_FALSE(rep.commission.has_points);
  EXPECT_EQ(rep.commission.nmax, 10); // clipped to the sigma model degree
}

TEST(ErrorReport, NumericalErrorFromPairsAndPrecomputed)
{
  const Model s = make_sigma(10, 1e-8);
  const std::vector<std::pair<double, double>> pairs = {{1.0 + 1e-13, 1.0}, {2.0, 2.0 + 2e-15}};
  ErrorReport rep = make_error_report(s, 360, kR, 9.81, pairs);
  ASSERT_TRUE(rep.numerical.available);
  EXPECT_NEAR(rep.numerical.max_relative_error, 1e-13, 1e-15);
  EXPECT_TRUE(rep.numerical.within_1e12);
  EXPECT_NEAR(rep.numerical.digits, 13.0, 0.01);
  EXPECT_NEAR(rep.numerical.digit_loss_bound, digit_loss_bound(360, PI / 2), 1e-12);

  set_relative_error(rep, 1e-10);
  EXPECT_FALSE(rep.numerical.within_1e12);
  EXPECT_NEAR(rep.numerical.digits, 10.0, 1e-9);

  set_relative_error(rep, 1e-3);
  EXPECT_NEAR(rep.numerical.digits, 3.0, 1e-9);

  // Digits are capped at 16 for exact agreement and for tiny errors.
  set_relative_error(rep, 0.0);
  EXPECT_NEAR(rep.numerical.digits, 16.0, 1e-12);
  set_relative_error(rep, 1e-20);
  EXPECT_NEAR(rep.numerical.digits, 16.0, 1e-12);
  EXPECT_TRUE(rep.numerical.within_1e12);

  const ErrorReport exact = make_error_report(s, 360, kR, 9.81, {{3.0, 3.0}});
  EXPECT_NEAR(exact.numerical.digits, 16.0, 1e-12);
}

TEST(ErrorReport, CommissionPointStatisticsAndGlobalValue)
{
  const int nmax = 30;
  const Model s = make_sigma(nmax, 1e-9);
  const double r = kR + 1.0e5;
  const Points pts = Points::from_spherical({10.0, 60.0, -35.0}, {20.0, -100.0, 170.0}, {r, r, r});
  const ErrorReport rep = make_error_report(s, nmax, r, 9.81, {}, &pts);

  ASSERT_TRUE(rep.commission.has_points);
  EXPECT_EQ(rep.commission.count, 3u);
  ErrorOptions eo;
  eo.nmin = 2;
  eo.nmax = nmax;
  const std::vector<double> ce = commission_error(s, pts, ErrorFunctional::ZetaEll, eo);
  const double lo = *std::min_element(ce.begin(), ce.end());
  const double hi = *std::max_element(ce.begin(), ce.end());
  const double mean = (ce[0] + ce[1] + ce[2]) / 3.0;
  EXPECT_NEAR(rep.commission.min, lo, 1e-12 * hi);
  EXPECT_NEAR(rep.commission.max, hi, 1e-12 * hi);
  EXPECT_NEAR(rep.commission.mean, mean, 1e-12 * hi);
  EXPECT_LE(rep.commission.min, rep.commission.mean);
  EXPECT_LE(rep.commission.mean, rep.commission.max);

  const double global = commission_global(s, ErrorFunctional::ZetaEll, r, 9.81, 2, nmax);
  EXPECT_NEAR(rep.commission.global_N, global, 1e-12 * global);

  // A run degree above the sigma model degree is clipped: the global value equals the one at the model degree.
  const ErrorReport clipped = make_error_report(s, 360, r, 9.81);
  EXPECT_EQ(clipped.commission.nmax, nmax);
  EXPECT_NEAR(clipped.commission.global_N, global, 1e-12 * global);
  EXPECT_FALSE(clipped.commission.has_points);
}

TEST(ErrorReport, ModelVarianceTailMatchesIndependentSum)
{
  // Coefficient model to degree 40, run at nmax = 20: OE^2 = sum_{21..40} (GM/(R gamma))^2 sigma_n^2
  //                                                    + sum_{41..n_tail} R^2 c_n / (gamma^2 (n-1)^2).
  const int nmax = 20, dmax = 40, ntail = 10800;
  const double gamma = 9.81;
  std::mt19937 rng(5);
  std::normal_distribution<double> gauss(0.0, 1.0);
  Model coef = Model::zeros(kGM, kR, dmax);
  for (int n = 2; n <= dmax; ++n)
    for (int m = 0; m <= n; ++m)
    {
      coef.C[idx_tri(n, m, dmax)] = 1e-8 / n * gauss(rng);
      coef.S[idx_tri(n, m, dmax)] = m == 0 ? 0.0 : 1e-8 / n * gauss(rng);
    }
  const Model sigma = make_sigma(dmax, 1e-9);

  ReportOptions opt;
  opt.coefficients = &coef;
  opt.n_tail = ntail;
  const ErrorReport rep = make_error_report(sigma, nmax, kR, gamma, {}, nullptr, opt);
  ASSERT_TRUE(rep.omission.has_tr_model);
  EXPECT_EQ(rep.omission.coefficient_nmax, dmax);

  const TscherningRapp p;
  double sum = 0.0;
  for (int n = nmax + 1; n <= dmax; ++n)
  {
    double s2 = 0.0;
    for (int m = 0; m <= n; ++m)
      s2 += coef.c(n, m) * coef.c(n, m) + coef.s(n, m) * coef.s(n, m);
    sum += std::pow(kGM / (kR * gamma), 2) * s2;
  }
  for (int n = dmax + 1; n <= ntail; ++n)
  {
    const double cn = p.A * (n - 1.0) / ((n - 2.0) * (n + p.B)) * std::pow(p.s(), n + 2.0) * 1e-10; // (m/s^2)^2
    sum += kR * kR * cn / (gamma * gamma * (n - 1.0) * (n - 1.0));
  }
  EXPECT_NEAR(rep.omission.tr_model, std::sqrt(sum), 1e-9 * std::sqrt(sum));

  // Without a coefficient model, or with the run degree at or above it, the model value is not computed.
  EXPECT_FALSE(make_error_report(sigma, nmax, kR, gamma).omission.has_tr_model);
  EXPECT_FALSE(make_error_report(sigma, dmax, kR, gamma, {}, nullptr, opt).omission.has_tr_model);
}

TEST(ErrorReport, InvalidArgumentsThrow)
{
  const Model s = make_sigma(10, 1e-8);
  EXPECT_THROW(make_error_report(s, 1, kR), std::invalid_argument);
  ReportOptions bad;
  bad.n_tail = 100;
  EXPECT_THROW(make_error_report(s, 360, kR, 9.81, {}, nullptr, bad), std::invalid_argument);
}

TEST(ErrorReport, TextAndJsonFormats)
{
  Model s = make_sigma(10, 1e-8);
  s.name = "EGM \"x\"\\";
  const Points pts = Points::from_spherical({0.0, 45.0}, {0.0, 90.0}, {kR, kR});
  ErrorReport rep = make_error_report(s, 360, kR, 9.81, {{1.0 + 1e-13, 1.0}}, &pts);

  const std::string txt = format_error_report(rep);
  EXPECT_NE(txt.find("Numerical error"), std::string::npos);
  EXPECT_NE(txt.find("digits retained"), std::string::npos);
  EXPECT_NE(txt.find("Kaula rule"), std::string::npos);
  EXPECT_NE(txt.find("Tscherning-Rapp"), std::string::npos);
  EXPECT_NE(txt.find("CE(N) global"), std::string::npos);
  EXPECT_NE(txt.find("points (2)"), std::string::npos);
  EXPECT_NE(txt.find("not computed"), std::string::npos); // model-variance TR is not available here
  EXPECT_EQ(txt.find("nan"), std::string::npos);

  const std::string json = error_report_json(rep);
  EXPECT_EQ(json.front(), '{');
  EXPECT_EQ(json.back(), '}');
  EXPECT_EQ(count_char(json, '{'), count_char(json, '}'));
  EXPECT_EQ(count_char(json, '['), count_char(json, ']'));
  EXPECT_NE(json.find("\"kaula_m\":"), std::string::npos);
  EXPECT_NE(json.find("\"within_1e-12\":true"), std::string::npos);
  EXPECT_NE(json.find("\"tr_model_m\":null"), std::string::npos);
  EXPECT_NE(json.find("\"name\":\"EGM \\\"x\\\"\\\\\""), std::string::npos); // escaped quote and backslash
  EXPECT_NE(json.find("\"points\":{\"count\":2"), std::string::npos);
  EXPECT_EQ(json.find("nan"), std::string::npos);
  EXPECT_EQ(json.find("inf"), std::string::npos);

  // Unavailable parts are written as null / not computed.
  const ErrorReport bare = make_error_report(s, 360, kR);
  const std::string json2 = error_report_json(bare);
  EXPECT_NE(json2.find("\"max_relative_error\":null"), std::string::npos);
  EXPECT_NE(json2.find("\"within_1e-12\":null"), std::string::npos);
  EXPECT_NE(json2.find("\"points\":null"), std::string::npos);
  EXPECT_NE(format_error_report(bare).find("E_rel: not computed"), std::string::npos);
}
