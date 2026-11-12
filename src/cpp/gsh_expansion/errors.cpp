// errors.cpp: implementation of errors.hpp (commission and omission errors).
#include "errors.hpp"
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include "legendre.hpp"

namespace gsh
{
  namespace
  {
    // Degree weight w_n of the functional at radius r, without the attenuation (R/r)^n (applied by the callers).
    double base_weight(ErrorFunctional f, int n, double GM, double r, double gamma)
    {
      switch (f)
      {
      case ErrorFunctional::T:
        return GM / r;
      case ErrorFunctional::ZetaEll:
        return GM / (r * gamma);
      case ErrorFunctional::GravityAnomaly:
        return GM / (r * r) * (n - 1);
      case ErrorFunctional::GravityDisturbance:
        return GM / (r * r) * (n + 1);
      }
      throw std::invalid_argument("unknown ErrorFunctional");
    }

    // Tscherning-Rapp anomaly degree variance c_n in (m/s^2)^2.
    double tr_anomaly_variance(int n, const TscherningRapp &p)
    {
      const double mgal2 = p.A * (n - 1.0) / ((n - 2.0) * (n + p.B)) * std::pow(p.s(), n + 2.0);
      return mgal2 * 1e-10; // (1 mGal = 1e-5 m/s^2)
    }
  }

  const char *error_functional_name(ErrorFunctional f)
  {
    switch (f)
    {
    case ErrorFunctional::T:
      return "T";
    case ErrorFunctional::ZetaEll:
      return "zeta_ell";
    case ErrorFunctional::GravityAnomaly:
      return "Dg_sa";
    case ErrorFunctional::GravityDisturbance:
      return "dg_radial";
    }
    return "?";
  }

  std::vector<double> commission_error(const Model &sigma, const Points &pts, ErrorFunctional f, const ErrorOptions &opts)
  {
    const int nmax = opts.nmax < 0 ? sigma.nmax : opts.nmax;
    if (nmax > sigma.nmax)
      throw std::invalid_argument("commission_error: nmax is larger than the degree of the sigma model");
    if (opts.nmin < 0 || opts.nmin > nmax)
      throw std::invalid_argument("commission_error: nmin must be in [0, nmax]");

    const std::int64_t np = pts.size();
    std::vector<double> out(static_cast<std::size_t>(np), 0.0);
    if (np == 0)
      return out;
    const legendre::Coefs coefs = legendre::make_coefs(nmax);
    const int nthreads = opts.threads > 0 ? opts.threads : omp_get_max_threads();

#pragma omp parallel num_threads(nthreads)
    {
      AlignedVec<double> P(tri_size(nmax));
      std::vector<double> cs(nmax + 1), sn(nmax + 1);
#pragma omp for schedule(dynamic, 4)
      for (std::int64_t i = 0; i < np; ++i)
      {
        const double lat = pts.lat[i], lon = pts.lon[i], r = pts.r[i];
        legendre::compute(std::sin(lat), std::cos(lat), coefs, P);
        for (int m = 0; m <= nmax; ++m)
        {
          cs[m] = std::cos(m * lon);
          sn[m] = std::sin(m * lon);
        }
        const double gamma = f == ErrorFunctional::ZetaEll ? opts.ellipsoid.normal_gravity(lat, r).gamma : 0.0;
        const double q = sigma.R / r;
        double qn = std::pow(q, opts.nmin);
        double var = 0.0;
        for (int n = opts.nmin; n <= nmax; ++n, qn *= q)
        {
          const double w = base_weight(f, n, sigma.GM, r, gamma) * qn;
          double sum = 0.0;
          for (int m = 0; m <= n; ++m)
          {
            const double p = P[idx_tri(n, m, nmax)];
            const double sc = sigma.c(n, m), ss = sigma.s(n, m);
            sum += p * p * (sc * sc * cs[m] * cs[m] + ss * ss * sn[m] * sn[m]);
          }
          var += w * w * sum;
        }
        out[static_cast<std::size_t>(i)] = std::sqrt(var);
      }
    }
    return out;
  }

  std::vector<double> commission_per_degree(const Model &sigma, ErrorFunctional f, double r, double gamma, int nmin, int nmax)
  {
    if (nmax > sigma.nmax)
      throw std::invalid_argument("commission_per_degree: nmax is larger than the degree of the sigma model");
    if (nmin < 0 || nmin > nmax)
      throw std::invalid_argument("commission_per_degree: nmin must be in [0, nmax]");
    std::vector<double> out(static_cast<std::size_t>(nmax) + 1, 0.0);
    const double q = sigma.R / r;
    for (int n = nmin; n <= nmax; ++n)
    {
      double s2 = 0.0;
      for (int m = 0; m <= n; ++m)
        s2 += sigma.c(n, m) * sigma.c(n, m) + sigma.s(n, m) * sigma.s(n, m);
      out[static_cast<std::size_t>(n)] = std::fabs(base_weight(f, n, sigma.GM, r, gamma)) * std::pow(q, n) * std::sqrt(s2);
    }
    return out;
  }

  std::vector<double> degree_variances(const Model &model)
  {
    std::vector<double> v(static_cast<std::size_t>(model.nmax) + 1, 0.0);
    for (int n = 0; n <= model.nmax; ++n)
    {
      double s2 = 0.0;
      for (int m = 0; m <= n; ++m)
        s2 += model.c(n, m) * model.c(n, m) + model.s(n, m) * model.s(n, m);
      v[static_cast<std::size_t>(n)] = s2;
    }
    return v;
  }

  std::vector<double> tscherning_rapp_variances(double GM, double R, int n_first, int n_last, const TscherningRapp &p)
  {
    if (n_first < 3 || n_last < n_first)
      throw std::invalid_argument("tscherning_rapp_variances: need 3 <= n_first <= n_last");
    std::vector<double> v(static_cast<std::size_t>(n_last) + 1, 0.0);
    const double g0 = GM / (R * R);
    for (int n = n_first; n <= n_last; ++n)
    {
      const double a = (n - 1.0) * g0; // c_n = a^2 sigma_n^2 at r = R
      v[static_cast<std::size_t>(n)] = tr_anomaly_variance(n, p) / (a * a);
    }
    return v;
  }

  std::vector<double> degree_variances_with_tail(const Model &model, int n_end, const TscherningRapp &p)
  {
    if (n_end <= model.nmax)
      throw std::invalid_argument("degree_variances_with_tail: n_end must exceed the degree of the model");
    std::vector<double> v = degree_variances(model);
    const std::vector<double> tail = tscherning_rapp_variances(model.GM, model.R, model.nmax + 1, n_end, p);
    v.resize(static_cast<std::size_t>(n_end) + 1, 0.0);
    for (int n = model.nmax + 1; n <= n_end; ++n)
      v[static_cast<std::size_t>(n)] = tail[static_cast<std::size_t>(n)];
    return v;
  }

  double omission_error(const std::vector<double> &sigma2, int nmax, ErrorFunctional f, double GM, double R, double r, double gamma)
  {
    if (nmax < 0)
      throw std::invalid_argument("omission_error: nmax must be >= 0");
    const int N = static_cast<int>(sigma2.size()) - 1;
    const double q = R / r;
    double sum = 0.0;
    for (int n = nmax + 1; n <= N; ++n)
    {
      const double w = base_weight(f, n, GM, r, gamma) * std::pow(q, n);
      sum += w * w * sigma2[static_cast<std::size_t>(n)];
    }
    return std::sqrt(sum);
  }
}

namespace gsh
{
  double commission_global(const Model &sigma, ErrorFunctional f, double r, double gamma, int nmin, int nmax)
  {
    const std::vector<double> d = commission_per_degree(sigma, f, r, gamma, nmin, nmax);
    double sum = 0.0;
    for (double x : d)
      sum += x * x;
    return std::sqrt(sum);
  }

  double omission_kaula_undulation(double R, int nmax)
  {
    if (nmax < 0)
      throw std::invalid_argument("omission_kaula_undulation: nmax must be >= 0");
    return R * 1e-5 / (nmax + 1.0);
  }

  double omission_tscherning_rapp_undulation(int nmax, int n_tail, double gamma, const TscherningRapp &p)
  {
    if (nmax < 2 || n_tail <= nmax)
      throw std::invalid_argument("omission_tscherning_rapp_undulation: need 2 <= nmax < n_tail");
    double sum = 0.0;
    for (int n = nmax + 1; n <= n_tail; ++n)
      sum += tr_anomaly_variance(n, p) / ((n - 1.0) * (n - 1.0));
    return p.R_E / gamma * std::sqrt(sum);
  }

  double relative_error(double calc, double ref)
  {
    const double d = std::fabs(calc - ref);
    return ref == 0.0 ? d : d / std::fabs(ref);
  }

  double digit_loss_bound(int M, double theta)
  {
    if (M < 1)
      throw std::invalid_argument("digit_loss_bound: M must be >= 1");
    const double tc = std::fabs(theta) < 1e-12 ? 1.0 : theta / std::tan(theta);
    return std::log10(static_cast<double>(M)) + std::log10(10.0 + 2.0 * std::fabs(tc));
  }
}
