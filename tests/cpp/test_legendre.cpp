// Acceptance tests for the Legendre functions, following Holmes and Featherstone (2002): closed forms,
// sum rule, orthonormality, polar boundary values, SFC / MSFC / ERA agreement and latitude derivatives.
#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include "legendre.hpp"
using namespace gsh;
using legendre::Coefs;

namespace
{
  std::vector<double> alf(double t, double u, const Coefs &c, bool modified = false)
  {
    std::vector<double> P(tri_size(c.nmax));
    if (modified)
      legendre::msfc(t, u, c, P);
    else
      legendre::sfc(t, u, c, P);
    return P;
  }
  // Maximum relative error of the sum rule sum_m P̄_nm² = 2n+1 over all degrees n <= nmax (Holmes and Featherstone (2002), item 3).
  double sum_rule_error(const std::vector<double> &P, int nmax)
  {
    double worst = 0.0;
    for (int n = 0; n <= nmax; ++n)
    {
      double s = 0.0;
      for (int m = 0; m <= n; ++m)
      {
        const double v = P[idx_tri(n, m, nmax)];
        s += v * v;
      }
      worst = std::max(worst, std::abs(s / (2.0 * n + 1.0) - 1.0));
    }
    return worst;
  }
  // Gauss-Legendre nodes and weights by Newton iteration, for the orthonormality test (Holmes and Featherstone (2002), item 4).
  void gauss(int N, std::vector<double> &x, std::vector<double> &w)
  {
    x.resize(N);
    w.resize(N);
    for (int i = 0; i < N; ++i)
    {
      double z = std::cos(PI * (i + 0.75) / (N + 0.5)), dp = 0;
      for (int it = 0; it < 100; ++it)
      {
        double p0 = 1, p1 = z;
        for (int k = 2; k <= N; ++k)
        {
          const double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = N * (z * p1 - p0) / (z * z - 1);
        const double dz = p1 / dp;
        z -= dz;
        if (std::abs(dz) < 1e-15)
          break;
      }
      x[i] = z;
      w[i] = 2 / ((1 - z * z) * dp * dp);
    }
  }
}

TEST(Legendre, ClosedForms)
{
  const auto c = legendre::make_coefs(2);
  for (double phi : {-1.2, -0.3, 0.0, 0.5, 1.4})
  {
    const double t = std::sin(phi), u = std::cos(phi);
    for (bool mod : {false, true})
    {
      const auto P = alf(t, u, c, mod);
      EXPECT_NEAR(P[idx_tri(0, 0, 2)], 1.0, 1e-14);
      EXPECT_NEAR(P[idx_tri(1, 0, 2)], std::sqrt(3.0) * t, 1e-14);
      EXPECT_NEAR(P[idx_tri(1, 1, 2)], std::sqrt(3.0) * u, 1e-14);
      EXPECT_NEAR(P[idx_tri(2, 0, 2)], std::sqrt(5.0) / 2 * (3 * t * t - 1), 1e-14);
      EXPECT_NEAR(P[idx_tri(2, 1, 2)], std::sqrt(15.0) * t * u, 1e-14);
      EXPECT_NEAR(P[idx_tri(2, 2, 2)], std::sqrt(15.0) / 2 * u * u, 1e-14);
    }
  }
}

TEST(Legendre, SmallDegrees)
{
  for (int nmax : {0, 1, 2})
  {
    const auto c = legendre::make_coefs(nmax);
    for (bool mod : {false, true})
    {
      const auto P = alf(0.3, std::sqrt(1 - 0.09), c, mod);
      EXPECT_NEAR(P[0], 1.0, 1e-15);
      EXPECT_LT(sum_rule_error(P, nmax), 1e-14);
    }
  }
}

TEST(Legendre, MatchesScipy)
{
  std::ifstream in(std::string(GSH_TEST_DATA) + "/scipy_alf_n60.txt");
  ASSERT_TRUE(in.good());
  const auto c = legendre::make_coefs(60);
  std::string line;
  std::getline(in, line);
  double last_t = 2, worst = 0;
  std::vector<double> P, Q;
  long count = 0;
  while (std::getline(in, line))
  {
    std::istringstream is(line);
    double t, ref;
    int n, m;
    is >> t >> n >> m >> ref;
    if (t != last_t)
    {
      P = alf(t, std::sqrt(1 - t * t), c);
      Q = alf(t, std::sqrt(1 - t * t), c, true);
      last_t = t;
    }
    const double scale = std::max(std::abs(ref), 1.0);
    worst = std::max(worst, std::abs(P[idx_tri(n, m, 60)] - ref) / scale);
    worst = std::max(worst, std::abs(Q[idx_tri(n, m, 60)] - ref) / scale);
    ++count;
  }
  EXPECT_EQ(count, 3 * 61 * 62 / 2);
  EXPECT_LT(worst, 1e-12);
}

TEST(Legendre, SumRule360)
{
  const auto c = legendre::make_coefs(360);
  for (double phi : {-1.5, -0.7, 0.01, 0.34, 1.1, 1.55})
    EXPECT_LT(sum_rule_error(alf(std::sin(phi), std::cos(phi), c), 360), 1e-12) << phi;
}

TEST(Legendre, Orthonormality)
{
  const int nmax = 40, N = 64;
  std::vector<double> x, w;
  gauss(N, x, w);
  const auto c = legendre::make_coefs(nmax);
  std::vector<std::vector<double>> P;
  for (double t : x)
    P.push_back(alf(t, std::sqrt(1 - t * t), c));
  double worst = 0;
  for (int m = 0; m <= nmax; ++m)
    for (int n = m; n <= nmax; ++n)
      for (int k = n; k <= nmax; ++k)
      {
        double s = 0;
        for (int i = 0; i < N; ++i)
          s += w[i] * P[i][idx_tri(n, m, nmax)] * P[i][idx_tri(k, m, nmax)];
        const double exact = (n == k) ? (m == 0 ? 2.0 : 4.0) : 0.0;
        worst = std::max(worst, std::abs(s - exact));
      }
  EXPECT_LT(worst, 1e-12);
}

TEST(Legendre, Poles)
{
  const int nmax = 400;
  const auto c = legendre::make_coefs(nmax);
  for (double t : {1.0, -1.0})
    for (bool mod : {false, true})
    {
      const auto P = alf(t, 0.0, c, mod);
      for (int m = 0; m <= nmax; ++m)
        for (int n = m; n <= nmax; ++n)
        {
          const double v = P[idx_tri(n, m, nmax)];
          ASSERT_TRUE(std::isfinite(v));
          if (m > 0)
            ASSERT_EQ(v, 0.0);
          else
            ASSERT_NEAR(std::abs(v), std::sqrt(2.0 * n + 1.0), 1e-11);
        }
    }
}

TEST(Legendre, MsfcHardCase2700)
{
  const int nmax = 2700;
  const auto c = legendre::make_coefs(nmax);
  const double u = 0.34, t = std::sqrt(1 - u * u);
  EXPECT_LT(sum_rule_error(alf(t, u, c, true), nmax), 1e-11);
  // SFC must fail here. This documents why MSFC is needed (Holmes and Featherstone 2002).
  EXPECT_GT(sum_rule_error(alf(t, u, c, false), nmax), 1e-3);
}

TEST(Legendre, MsfcEqualsSfcAt1800)
{
  const int nmax = 1800;
  const auto c = legendre::make_coefs(nmax);
  const double t = std::sin(0.8), u = std::cos(0.8);
  const auto A = alf(t, u, c), B = alf(t, u, c, true);
  double worst = 0;
  for (std::size_t i = 0; i < A.size(); ++i)
    worst = std::max(worst, std::abs(A[i] - B[i]));
  EXPECT_LT(worst, 1e-9);
}

// ERA (literal extended-range arithmetic) as the precision reference for MSFC, with no scale factors.
namespace
{
  std::vector<double> era_alf(double t, double u, const Coefs &c)
  {
    std::vector<double> P(tri_size(c.nmax));
    legendre::era(t, u, c, P);
    return P;
  }
  // Maximum difference between two triangles, relative to the scale sqrt(2n+1) of each degree.
  double scaled_diff(const std::vector<double> &A, const std::vector<double> &B, int nmax)
  {
    double worst = 0.0;
    for (int m = 0; m <= nmax; ++m)
      for (int n = m; n <= nmax; ++n)
      {
        const auto i = idx_tri(n, m, nmax);
        worst = std::max(worst, std::abs(A[i] - B[i]) / std::sqrt(2.0 * n + 1.0));
      }
    return worst;
  }
}

TEST(Legendre, EraMatchesSfcAtLowDegree)
{
  const int nmax = 360;
  const auto c = legendre::make_coefs(nmax);
  for (const double lat : {0.0, 0.7, 1.4, 1.5707})
    EXPECT_LT(scaled_diff(era_alf(std::sin(lat), std::cos(lat), c), alf(std::sin(lat), std::cos(lat), c), nmax), 1e-11) << lat;
}

TEST(Legendre, EraHardCase2700)
{
  const int nmax = 2700;
  const auto c = legendre::make_coefs(nmax);
  const double u = 0.34, t = std::sqrt(1 - u * u);
  const auto E = era_alf(t, u, c);
  EXPECT_LT(sum_rule_error(E, nmax), 1e-11);
  EXPECT_LT(scaled_diff(E, alf(t, u, c, true), nmax), 1e-12);
}

TEST(Legendre, EraVsMsfcNearPole)
{
  // Latitudes 89.99° and 89.9999° (u = cos φ ≈ 1.7e-4 and 1.7e-6) with nmax = 2190: IEEE underflow destroys SFC,
  // so MSFC must agree with ERA here.
  const int nmax = 2190;
  const auto c = legendre::make_coefs(nmax);
  for (const double u : {1.7e-4, 1.7e-6})
  {
    const double t = std::sqrt(1 - u * u);
    const auto E = era_alf(t, u, c);
    EXPECT_LT(sum_rule_error(E, nmax), 1e-9) << u;
    EXPECT_LT(scaled_diff(E, alf(t, u, c, true), nmax), 1e-9) << u;
  }
}

TEST(Legendre, Derivatives)
{
  const int nmax = 120;
  const auto c = legendre::make_coefs(nmax);
  const double h = 1e-4;
  for (double phi : {-1.1, 0.2, 0.9})
  {
    const auto P = alf(std::sin(phi), std::cos(phi), c);
    const auto Pp = alf(std::sin(phi + h), std::cos(phi + h), c), Pm = alf(std::sin(phi - h), std::cos(phi - h), c);
    const auto Ppp = alf(std::sin(phi + 2 * h), std::cos(phi + 2 * h), c), Pmm = alf(std::sin(phi - 2 * h), std::cos(phi - 2 * h), c);
    std::vector<double> dP(P.size()), ddP(P.size());
    legendre::dphi(nmax, P, dP);
    legendre::dphi(nmax, dP, ddP);
    double err1 = 0, err2 = 0;
    for (std::size_t i = 0; i < P.size(); ++i)
    {
      const double fd1 = (8 * (Pp[i] - Pm[i]) - (Ppp[i] - Pmm[i])) / (12 * h);
      const double fd2 = (16 * (Pp[i] + Pm[i]) - (Ppp[i] + Pmm[i]) - 30 * P[i]) / (12 * h * h);
      err1 = std::max(err1, std::abs(dP[i] - fd1));
      err2 = std::max(err2, std::abs(ddP[i] - fd2));
    }
    EXPECT_LT(err1, 1e-4);
    EXPECT_LT(err2, 2e-2);
  }
  // Exact pole: all derivatives must be finite (the identity never divides by cos φ).
  const auto Pn = alf(1.0, 0.0, c);
  std::vector<double> dP(Pn.size()), ddP(Pn.size());
  legendre::dphi(nmax, Pn, dP);
  legendre::dphi(nmax, dP, ddP);
  for (double v : ddP)
    ASSERT_TRUE(std::isfinite(v));
}

TEST(Legendre, ColumnDerivativeMatchesIdentity)
{
  // The column-wise derivative used by the kernels, (f P(n-1,m) - n t P(n,m)) / u (Fukushima 2012, II), must agree with the
  // non-singular identity of Bucha and Janák (2019), which is what dphi() implements.
  const int nmax = 300;
  const auto c = legendre::make_coefs(nmax);
  const double phi = 0.7, t = std::sin(phi), u = std::cos(phi);
  const auto P = alf(t, u, c);
  std::vector<double> dP(P.size());
  legendre::dphi(nmax, P, dP);
  double worst = 0;
  for (int m = 0; m <= nmax; ++m)
    for (int n = m; n <= nmax; ++n)
    {
      const auto i = idx_tri(n, m, nmax);
      const double prev = n > m ? P[idx_tri(n - 1, m, nmax)] : 0.0;
      const double v = (c.f[i] * prev - n * t * P[i]) / u;
      worst = std::max(worst, std::abs(v - dP[i]) / (n + 1.0));
    }
  EXPECT_LT(worst, 1e-11);
}

// SFC at nmax = 2190 near latitude 60 deg: the sectorial value underflows into the denormal range, and the forward
// recurrence used to amplify its rounding error to ~1e23 where the true values are ~1e-36. SFC is not accurate there
// (Holmes and Featherstone 2002), but it must underflow to zero, never return values larger than the true ones.
TEST(Legendre, SfcUnderflowsToZeroInsteadOfBlowingUp)
{
  const int nmax = 2190;
  const auto c = legendre::make_coefs(nmax);
  for (const double lat_deg : {-60.2, -59.8, -59.3, 59.3, 59.8, 60.2})
  {
    const double lat = lat_deg * PI / 180.0, t = std::sin(lat), u = std::cos(lat);
    const auto P = alf(t, u, c, false);
    // |P̄_nm| <= sqrt(2n+1) always holds (sum rule), so nothing may exceed sqrt(2 nmax + 1) by more than rounding.
    const double bound = std::sqrt(2.0 * nmax + 1.0) * (1.0 + 1e-6);
    for (const double v : P)
      ASSERT_LE(std::fabs(v), bound) << "latitude " << lat_deg;
  }
}
