// Acceptance protocol of this project, with the derivative formulas of Fukushima (2012, II) and Bosch (2000) and the
// functional formulas of Bucha and Janák (2013, 2019). CPU only.
// The "[METRIC]" lines feed docs/validacion.md.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <random>
#include <string>
#include <vector>
#include "functionals.hpp"
#include "height_model.hpp"
#include "kernels.hpp"
#include "legendre.hpp"
using namespace gsh;
using F = Functional;
namespace fs = std::filesystem;

namespace
{
  void metric(const std::string &name, double v) { std::printf("[METRIC] %s %.3e\n", name.c_str(), v); }

  // ---------- Legendre: sum identities (the Holmes and Featherstone 2002 identities) ----------
  enum class Alg
  {
    Sfc,
    Msfc,
    Era
  };
  struct Sums
  {
    double d[4];
  };

  // Relative error of I^(0..3) against the closed-form values of the protocol (sum identities of Holmes and Featherstone 2002), at the point (t, u).
  Sums sum_rules(int M, double t, double u, Alg alg)
  {
    const auto c = legendre::make_coefs(M);
    const auto n = tri_size(M);
    std::vector<double> P(n), D1(n), D2(n), D3(n);
    switch (alg)
    {
    case Alg::Sfc:
      legendre::sfc(t, u, c, P);
      break;
    case Alg::Msfc:
      legendre::msfc(t, u, c, P);
      break;
    case Alg::Era:
      legendre::era(t, u, c, P);
      break;
    }
    legendre::dphi(M, P, D1);
    legendre::dphi(M, D1, D2);
    legendre::dphi(M, D2, D3);
    const std::vector<double> *a[4] = {&P, &D1, &D2, &D3};
    long double s[4] = {0, 0, 0, 0};
    for (int k = 0; k < 4; ++k)
      for (int m = 0; m <= M; ++m)
      {
        double col = 0.0;
        const auto off = col_offset(m, M);
        for (int i = 0; i <= M - m; ++i)
          col += (*a[k])[off + i] * (*a[k])[off + i];
        s[k] += col;
      }
    const long double x = M;
    const long double i0 = (x + 1) * (x + 1), i1 = x * (x + 2) * (x + 1) * (x + 1) / 4;
    const long double i2 = (x * x + 2 * x - 1) / 2 * i1;
    const long double i3 = (5 * x * x * x * x + 20 * x * x * x + 5 * x * x - 30 * x + 16) / 16 * i1;
    const long double th[4] = {i0, i1, i2, i3};
    Sums out{};
    for (int k = 0; k < 4; ++k)
      out.d[k] = static_cast<double>(std::abs(s[k] / th[k] - 1.0L));
    return out;
  }
  Sums sum_rules_lat(int M, double lat_deg, Alg alg)
  {
    const double phi = deg2rad(lat_deg);
    return sum_rules(M, std::sin(phi), std::cos(phi), alg);
  }
  double worst(const Sums &s) { return *std::max_element(s.d, s.d + 4); }
}

TEST(ProtocolLegendre, SumRulesDoublePrecision360)
{
  for (double lat : {0.0, 30.0, 60.0, 85.0, 89.0})
  {
    const double w = worst(sum_rules_lat(360, lat, Alg::Sfc));
    metric("sumrule_M360_lat" + std::to_string(static_cast<int>(lat)), w);
    EXPECT_LT(w, 1e-12) << "lat " << lat; // acceptance criterion of this project: < 1e-12 for M <= 360 in double precision
  }
}

TEST(ProtocolLegendre, SumRulesExtendedRange2700)
{
  for (double lat : {0.0, 30.0, 60.0, 85.0, 89.0})
  {
    const auto e = sum_rules_lat(2700, lat, Alg::Era), m = sum_rules_lat(2700, lat, Alg::Msfc);
    metric("sumrule_M2700_era_lat" + std::to_string(static_cast<int>(lat)), worst(e));
    EXPECT_LT(worst(e), 1e-11) << "lat " << lat; // acceptance criterion of this project: < 1e-11 in extended-range arithmetic
    for (int k = 0; k < 4; ++k)
      EXPECT_NEAR(e.d[k], m.d[k], 1e-13) << "ERA and MSFC must agree, k=" << k << " lat " << lat;
  }
}

TEST(ProtocolLegendre, SumRulesExtendedRange5400)
{
  for (double lat : {45.0, 80.0})
  {
    const double w = worst(sum_rules_lat(5400, lat, Alg::Era));
    metric("sumrule_M5400_era_lat" + std::to_string(static_cast<int>(lat)), w);
    EXPECT_LT(w, 1e-11) << "lat " << lat;
  }
}

// Near the pole the error is set by the resolution of t = sin φ in double precision, not by the algorithm:
// changing t by 1 ulp changes the identity by about 4e-10 at M = 2700 (the protocol asks for 1e-11, which is
// unreachable there).
TEST(ProtocolLegendre, NearPoleErrorComesFromInputResolution)
{
  const int M = 2700;
  const double phi = deg2rad(89.99), t = std::sin(phi), u = std::cos(phi);
  const double base = sum_rules(M, t, u, Alg::Era).d[0];
  const double up = sum_rules(M, std::nextafter(t, 2.0), u, Alg::Era).d[0];
  const double down = sum_rules(M, std::nextafter(t, 0.0), u, Alg::Era).d[0];
  metric("sumrule_M2700_lat89.99_base", base);
  metric("sumrule_M2700_lat89.99_plus1ulp", up);
  metric("sumrule_M2700_lat89.99_minus1ulp", down);
  EXPECT_LT(base, 2e-9);
  EXPECT_LT(up, 5e-9);
  EXPECT_LT(down, 5e-9);
  EXPECT_GT(std::max(up, down), 1e-10); // the input resolution, not the exponent range, dominates the error
}

// ---------- Legendre derivatives from the derivative formulas, implemented independently here ----------
TEST(ProtocolLegendre, DerivativeFormulasOfTheDocument)
{
  const int nmax = 360;
  const auto c = legendre::make_coefs(nmax);
  const auto size = tri_size(nmax);
  for (double lat : {30.0, 60.0, -45.0, 85.0})
  {
    const double phi = deg2rad(lat), t = std::sin(phi), u = std::cos(phi);
    std::vector<double> P(size), D1(size), D2(size), D3(size);
    legendre::sfc(t, u, c, P);
    legendre::dphi(nmax, P, D1);
    legendre::dphi(nmax, D1, D2);
    legendre::dphi(nmax, D2, D3);
    double w1 = 0, w2 = 0, w3 = 0;
    for (int m = 0; m <= nmax; ++m)
      for (int n = m; n <= nmax; ++n)
      {
        const double N = n, M = m;
        const double Pn = P[idx_tri(n, m, nmax)];
        const double Pm1 = (n > m) ? P[idx_tri(n - 1, m, nmax)] : 0.0;
        const double e = (n > m) ? std::sqrt((2 * N + 1) * (N - M) * (N + M) / (2 * N - 1)) : 0.0;
        // Derivative formula (written for m >= 1), applied with the same sign also for m = 0 (see DocErrataZonalSign).
        const double p1 = -N * t / u * Pn + e / u * Pm1;
        const double p2 = t / u * p1 + (M * M / (u * u) - N * (N + 1)) * Pn;
        const double p3 = t / u * p2 + ((M * M + 1) / (u * u) - N * (N + 1)) * p1 + 2 * M * M * t / (u * u * u) * Pn;
        const double s1 = std::sqrt(2 * N + 1) * (N + 1), s2 = s1 * (N + 1), s3 = s2 * (N + 1);
        const auto i = idx_tri(n, m, nmax);
        w1 = std::max(w1, std::abs(D1[i] - p1) / s1);
        w2 = std::max(w2, std::abs(D2[i] - p2) / s2);
        w3 = std::max(w3, std::abs(D3[i] - p3) / s3);
      }
    metric("deriv_formula_lat" + std::to_string(static_cast<int>(lat)) + "_d1", w1);
    metric("deriv_formula_lat" + std::to_string(static_cast<int>(lat)) + "_d2", w2);
    metric("deriv_formula_lat" + std::to_string(static_cast<int>(lat)) + "_d3", w3);
    EXPECT_LT(w1, 1e-12) << lat;
    EXPECT_LT(w2, 1e-12) << lat;
    EXPECT_LT(w3, 1e-11) << lat;
  }
}

// Sign check: the m >= 1 derivative formula, applied literally at m = 0, has the opposite sign to the closed-form value
// P^(1)_{1,0} = +sqrt(3) u. The code uses the m >= 1 form, which reproduces the closed-form value.
TEST(ProtocolLegendre, DocErrataZonalSign)
{
  const int nmax = 3;
  const auto c = legendre::make_coefs(nmax);
  const double phi = deg2rad(40.0), t = std::sin(phi), u = std::cos(phi);
  std::vector<double> P(tri_size(nmax)), D1(tri_size(nmax));
  legendre::sfc(t, u, c, P);
  legendre::dphi(nmax, P, D1);
  EXPECT_NEAR(D1[idx_tri(1, 0, nmax)], std::sqrt(3.0) * u, 1e-14); // value stated by the document and produced by the code
  const double doc_line = 1.0 * t / u * P[idx_tri(1, 0, nmax)] - std::sqrt(3.0) / u * P[idx_tri(0, 0, nmax)];
  EXPECT_NEAR(doc_line, -std::sqrt(3.0) * u, 1e-14); // the formula as written gives the opposite sign
}

// ---------- Functionals (Bucha and Janák 2013, 2019) ----------
namespace
{
  Model random_model(int nmax, unsigned seed, double amp0 = 1e-6)
  {
    Model m = Model::zeros(3.986004415e14, 6378136.3, nmax);
    std::mt19937_64 g(seed);
    std::normal_distribution<double> nd;
    m.C[idx_tri(0, 0, nmax)] = 1.0;
    m.C[idx_tri(2, 0, nmax)] = -4.84e-4;
    for (int mm = 0; mm <= nmax; ++mm)
      for (int n = std::max(mm, 2); n <= nmax; ++n)
      {
        const double amp = amp0 / (n * n);
        m.C[idx_tri(n, mm, nmax)] += amp * nd(g);
        if (mm > 0)
          m.S[idx_tri(n, mm, nmax)] = amp * nd(g);
      }
    return m;
  }
  Points random_points(int count, int seed, double lat_max = 75.0)
  {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> lat(-lat_max, lat_max), lon(-180, 180), h(0, 5000);
    Points p;
    for (int i = 0; i < count; ++i)
    {
      p.lat.push_back(deg2rad(lat(g)));
      p.lon.push_back(deg2rad(lon(g)));
      p.r.push_back(6378000.0 + h(g));
    }
    return p;
  }
  AlignedVec<double> eval(const Model &m, const Points &p, F f)
  {
    FunctionalRequest rq;
    rq.outputs = {f};
    return compute_functionals(m, p, rq).values.at(f);
  }
}

// ξ = -(1/(rγ)) ∂T/∂φ, η = -(1/(rγ cosφ)) ∂T/∂λ and ζ_ell = T/γ, checked against central finite differences of T
// (Bucha and Janák 2019).
TEST(ProtocolFunctionals, DeflectionsAndHeightAnomalyAgainstFiniteDifferences)
{
  const Model model = random_model(60, 11);
  const Points base = random_points(12, 4);
  const auto ell = Ellipsoid::GRS80();
  FunctionalRequest rq;
  rq.outputs = {F::xi, F::eta, F::Theta, F::zeta_ell, F::T};
  const auto res = compute_functionals(model, base, rq);
  const double d = 1e-6;
  const auto shifted = [&](double dlat, double dlon)
  {
    Points p = base;
    for (std::int64_t i = 0; i < p.size(); ++i)
    {
      p.lat[i] += dlat;
      p.lon[i] += dlon;
    }
    return eval(model, p, F::T);
  };
  const auto Tn = shifted(d, 0), Ts = shifted(-d, 0), Te = shifted(0, d), Tw = shifted(0, -d);
  double scale = 0;
  for (std::int64_t i = 0; i < base.size(); ++i)
    scale = std::max(scale, std::abs(res.values.at(F::xi)[i]));
  double wxi = 0, weta = 0, wze = 0, wth = 0;
  for (std::int64_t i = 0; i < base.size(); ++i)
  {
    const double r = base.r[i], gamma = ell.normal_gravity(base.lat[i], r).gamma, cl = std::cos(base.lat[i]);
    const double xi = -(Tn[i] - Ts[i]) / (2 * d) / (r * gamma);
    const double eta = -(Te[i] - Tw[i]) / (2 * d) / (r * gamma * cl);
    wxi = std::max(wxi, std::abs(res.values.at(F::xi)[i] - xi));
    weta = std::max(weta, std::abs(res.values.at(F::eta)[i] - eta));
    wze = std::max(wze, std::abs(res.values.at(F::zeta_ell)[i] - res.values.at(F::T)[i] / gamma));
    wth = std::max(wth, std::abs(res.values.at(F::Theta)[i] - std::hypot(res.values.at(F::xi)[i], res.values.at(F::eta)[i])));
  }
  metric("xi_vs_fd_relative", wxi / scale);
  metric("eta_vs_fd_relative", weta / scale);
  EXPECT_LT(wxi, 1e-6 * scale);
  EXPECT_LT(weta, 1e-6 * scale);
  EXPECT_LT(wze, 1e-12);
  EXPECT_LT(wth, 1e-18);
}

// N = (T - 2πGρH²)/γ and ζ = T/γ - δg_sa (H+N)/γ (Bucha and Janák 2019, PDF 75-78), assembled here independently from a
// harmonic topography model.
TEST(ProtocolFunctionals, GeoidUndulationAndHeightAnomalyAssembly)
{
  const Model model = random_model(60, 21);
  const int hmax = 20;
  Model hm = Model::zeros(1.0, 6378136.3, hmax);
  std::mt19937_64 g(5);
  std::normal_distribution<double> nd;
  for (int m = 0; m <= hmax; ++m)
    for (int n = std::max(m, 1); n <= hmax; ++n)
    {
      hm.C[idx_tri(n, m, hmax)] = 300.0 / (n * n) * nd(g);
      if (m > 0)
        hm.S[idx_tri(n, m, hmax)] = 300.0 / (n * n) * nd(g);
    }
  const Points pts = random_points(10, 8);
  const auto ell = Ellipsoid::GRS80();

  FunctionalRequest rq;
  rq.outputs = {F::N, F::zeta};
  rq.height = height_from_model(hm);
  const auto res = compute_functionals(model, pts, rq);

  // Foot points on the ellipsoid, and H by direct summation with the ERA Legendre functions.
  Points foot = pts;
  for (std::int64_t i = 0; i < pts.size(); ++i)
    foot.r[i] = ell.radius(pts.lat[i]);
  const auto T = eval(model, foot, F::T), dg = eval(model, foot, F::dg_sa);
  const auto coefs = legendre::make_coefs(hmax);
  double wN = 0, wZ = 0;
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    std::vector<double> P(tri_size(hmax));
    legendre::era(std::sin(pts.lat[i]), std::cos(pts.lat[i]), coefs, P);
    double H = 0;
    for (int m = 0; m <= hmax; ++m)
      for (int n = m; n <= hmax; ++n)
        H += (hm.c(n, m) * std::cos(m * pts.lon[i]) + hm.s(n, m) * std::sin(m * pts.lon[i])) * P[idx_tri(n, m, hmax)];
    const double gamma = ell.normal_gravity(pts.lat[i], foot.r[i]).gamma;
    const double N = (T[i] - 2.0 * PI * 6.67259e-11 * 2670.0 * H * H) / gamma;
    const double zeta = T[i] / gamma - dg[i] * (H + N) / gamma;
    wN = std::max(wN, std::abs(res.values.at(F::N)[i] - N));
    wZ = std::max(wZ, std::abs(res.values.at(F::zeta)[i] - zeta));
  }
  metric("N_assembly_abs_m", wN);
  metric("zeta_assembly_abs_m", wZ);
  EXPECT_LT(wN, 1e-9);
  EXPECT_LT(wZ, 1e-9);
}

// ---------- Real models (EGM2008, XGM2019e), degree 2190, CPU ----------
namespace
{
  std::string model_path(const std::string &file)
  {
    // Extracted models live in data/cache/modelos/<file stem>/<file>.
    const std::string stem = file.substr(0, file.find_last_of('.'));
    const std::string p = std::string(GSH_REPO_DATA) + "/cache/modelos/" + stem + "/" + file;
    return fs::exists(p) ? p : std::string{};
  }
  // Neumaier compensated summation, so that the reference sum does not limit the accuracy being measured.
  struct Kahan
  {
    double s = 0, c = 0;
    void add(double v)
    {
      const double t = s + v;
      c += (std::abs(s) >= std::abs(v)) ? (s - t) + v : (v - t) + s;
      s = t;
    }
    double value() const { return s + c; }
  };
}

class RealModel : public ::testing::TestWithParam<const char *>
{
};

// Raw potential (degree >= 2) from the fast kernel, against a direct sum with extended-range Legendre functions.
TEST_P(RealModel, FastKernelMatchesExtendedRangeSum)
{
  const std::string path = model_path(GetParam());
  if (path.empty())
    GTEST_SKIP() << GetParam() << " is not in data/";
  const int nmax = 2190;
  const Model model = load_gfc(path, nmax);
  ASSERT_EQ(model.nmax, nmax);
  const auto coefs = legendre::make_coefs(nmax);
  const auto prep = prepare(model, 2, nmax);

  struct Case
  {
    double lat, tol;
  };
  const std::vector<Case> cases = {{0.0, 1e-12}, {30.5, 1e-12}, {45.0, 1e-12}, {80.0, 1e-12}, {85.0, 1e-12},
                                   {-85.0, 1e-12}, {89.9, 1e-10}, {89.99, 2e-9}, {-89.99, 2e-9}};
  Points pts;
  for (const auto &c : cases)
    for (double lon : {13.7, 200.3})
    {
      pts.lat.push_back(deg2rad(c.lat));
      pts.lon.push_back(deg2rad(lon));
      pts.r.push_back(model.R + 1500.0);
    }
  const auto raw = synth_points(prep, pts, 0);

  std::vector<double> ref(pts.size());
  double scale = 0;
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    std::vector<double> P(tri_size(nmax));
    legendre::era(std::sin(pts.lat[i]), std::cos(pts.lat[i]), coefs, P);
    const double q = model.R / pts.r[i];
    Kahan acc;
    for (int m = 0; m <= nmax; ++m)
    {
      const double cm = std::cos(m * pts.lon[i]), sm = std::sin(m * pts.lon[i]);
      for (int n = std::max(m, 2); n <= nmax; ++n)
        acc.add(std::pow(q, n) * (model.c(n, m) * cm + model.s(n, m) * sm) * P[idx_tri(n, m, nmax)]);
    }
    ref[i] = acc.value();
    scale = std::max(scale, std::abs(ref[i]));
  }
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const double rel = std::abs(raw.f[kernels::F_S0][i] - ref[i]) / scale;
    const double lat = cases[i / 2].lat;
    metric(std::string(GetParam()) + "_S0_lat" + std::to_string(lat), rel);
    EXPECT_LT(rel, cases[i / 2].tol) << GetParam() << " lat " << lat << " lon " << rad2deg(pts.lon[i]);
  }
}

// With real coefficients and n up to 2190, the trace of the disturbing tensor and of V must vanish (Laplace).
TEST_P(RealModel, TensorTraceVanishesAtUltraHighDegree)
{
  const std::string path = model_path(GetParam());
  if (path.empty())
    GTEST_SKIP() << GetParam() << " is not in data/";
  const Model model = load_gfc(path, 2190);
  Points pts;
  for (double lat : {0.0, 38.2, -61.0, 80.0, 89.0})
  {
    pts.lat.push_back(deg2rad(lat));
    pts.lon.push_back(deg2rad(254.4));
    pts.r.push_back(model.R + 5376.0);
  }
  FunctionalRequest rq;
  rq.outputs = {F::V_xx, F::V_yy, F::V_zz, F::T_xx, F::T_yy, F::T_zz, F::Dg_sa, F::dg_sa, F::T};
  const auto res = compute_functionals(model, pts, rq);
  double wV = 0, wT = 0, wRel = 0;
  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    const double sV = std::abs(res.values.at(F::V_zz)[i]);
    wV = std::max(wV, std::abs(res.values.at(F::V_xx)[i] + res.values.at(F::V_yy)[i] + res.values.at(F::V_zz)[i]) / sV);
    const double sT = std::max({std::abs(res.values.at(F::T_xx)[i]), std::abs(res.values.at(F::T_yy)[i]), std::abs(res.values.at(F::T_zz)[i])});
    wT = std::max(wT, std::abs(res.values.at(F::T_xx)[i] + res.values.at(F::T_yy)[i] + res.values.at(F::T_zz)[i]) / sT);
    const double lhs = res.values.at(F::Dg_sa)[i], rhs = res.values.at(F::dg_sa)[i] - 2.0 * res.values.at(F::T)[i] / pts.r[i];
    wRel = std::max(wRel, std::abs(lhs - rhs) / std::abs(res.values.at(F::dg_sa)[i]));
  }
  metric(std::string(GetParam()) + "_trace_V_rel", wV);
  metric(std::string(GetParam()) + "_trace_T_rel", wT);
  metric(std::string(GetParam()) + "_Dg_identity_rel", wRel);
  EXPECT_LT(wV, 1e-9);
  EXPECT_LT(wT, 1e-9);
  EXPECT_LT(wRel, 1e-10);
}

INSTANTIATE_TEST_SUITE_P(Models, RealModel, ::testing::Values("EGM2008.gfc", "XGM2019e_2159.gfc"));
