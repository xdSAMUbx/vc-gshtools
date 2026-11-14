#include "legendre.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace gsh::legendre
{
  Coefs make_coefs(int nmax)
  {
    if (nmax < 0)
      throw std::invalid_argument("make_coefs: nmax negativo");
    Coefs c;
    c.nmax = nmax;
    const auto size = tri_size(nmax);
    c.a.assign(size, 0.0);
    c.b.assign(size, 0.0);
    c.f.assign(size, 0.0);
    c.d.assign(nmax + 1, 1.0);
    for (int m = 1; m <= nmax; ++m)
      c.d[m] = (m == 1) ? std::sqrt(3.0) : std::sqrt((2.0 * m + 1.0) / (2.0 * m));

// Each column m writes only its own packed slice, so the outer loop is race-free.
#pragma omp parallel for schedule(dynamic, 16)
    for (int m = 0; m <= nmax; ++m)
    {
      const double M = m;
      const auto off = col_offset(m, nmax);
      // a: sqrt((2N+1)(2N-1)/((N-M)(N+M))); f: sqrt((N-M)(N+M)(2N+1)/(2N-1)); b: only for n >= m+2 (Holmes and Featherstone 2002).
      for (int n = m + 1; n <= nmax; ++n)
      {
        const double N = n;
        const auto i = off + (n - m);
        const double nmm = N - M, npm = N + M;
        c.a[i] = std::sqrt((2.0 * N + 1.0) * (2.0 * N - 1.0) / (nmm * npm));
        c.f[i] = std::sqrt(nmm * npm * (2.0 * N + 1.0) / (2.0 * N - 1.0));
        if (n >= m + 2)
          c.b[i] = std::sqrt((2.0 * N + 1.0) * (npm - 1.0) * (nmm - 1.0) / (nmm * npm * (2.0 * N - 3.0)));
      }
    }
    return c;
  }

  namespace
  {
    // Requires a buffer P of at least tri_size(nmax) elements.
    void check(const Coefs &c, std::span<const double> P)
    {
      if (static_cast<std::int64_t>(P.size()) < tri_size(c.nmax))
        throw std::invalid_argument("legendre: buffer P insuficiente");
    }

    // Analytic pole shortcut shared by sfc(), msfc() and era(). At t = ±1 (u = 0) the fully normalised functions are
    //   P̄_n0(±1) = (±1)^n sqrt(2n+1)   and   P̄_nm(±1) = 0 for m >= 1.
    // The zonal column is then exact. Evaluating it with the three-term recurrence at t = ±1 is not: the
    // characteristic equation has a double root there, so the rounding error of each step is amplified linearly in
    // the remaining steps and accumulates as O(M² ε) (3.5e-11 in the sum rule at M = 1800, 9.0e-11 at M = 2700).
    // Applicability: |(|t| - 1)| < 1e-15 as requested, AND u < 1e-15. The second guard matters: |t| = 1 - 1e-15
    // corresponds to u = 4.5e-8 (about 0.3 m from the pole), where P̄_n1 is still of order 1e-3 at n = 2700 and
    // zeroing it would be an error; the pole value is used only when u itself is negligible (cos(π/2) = 6e-17 passes).
    // Returns true if P was filled.
    bool pole_shortcut(double t, double u, const Coefs &c, std::span<double> P)
    {
      if (!(std::abs(std::abs(t) - 1.0) < 1e-15 && std::abs(u) < 1e-15))
        return false;
      std::fill(P.begin(), P.begin() + tri_size(c.nmax), 0.0);
      double sign = 1.0; // (±1)^n: alternates only at the south pole
      for (int n = 0; n <= c.nmax; ++n)
      {
        P[idx_tri(n, 0, c.nmax)] = sign * std::sqrt(2.0 * n + 1.0);
        if (t < 0.0)
          sign = -sign;
      }
      return true;
    }
  }

  void sfc(double t, double u, const Coefs &c, std::span<double> P)
  {
    check(c, P);
    if (pole_shortcut(t, u, c, P))
      return;
    const int nmax = c.nmax;
    // sect = P(m,m) = prod_{k=1..m} d_k u (Holmes and Featherstone 2002). It is the value that underflows for large m and u.
    double sect = 1.0, lsum = 0.0; // lsum = log(sect), carried only to rebuild sect in the denormal range
    for (int m = 0; m <= nmax; ++m)
    {
      if (m > 0)
      {
        const double step = c.d[m] * u;
        lsum += std::log(step);
        sect *= step;
        // In the denormal range the repeated product stalls: when the factor d_m u is a little above 0.5 (latitude
        // near 60 deg) the product of the smallest denormal rounds back up to itself and never reaches zero. Every
        // later column is then seeded with 4.9e-324 instead of ~0, and the forward recurrence amplifies that to ~1e23
        // (true values ~1e-36 at nmax = 2190). Below the smallest normal number the value is taken from the log
        // sum instead, which rounds to the correctly rounded denormal and to zero when u^m underflows completely, as a
        // direct power does (Holmes and Featherstone 2002).
        if (sect < std::numeric_limits<double>::min())
          sect = std::exp(lsum);
      }
      const auto off = col_offset(m, nmax);
      const double *ca = c.a.data() + off, *cb = c.b.data() + off;
      double *col = P.data() + off;
      // col[j] holds P(n, m) with n = m + j. Start at the sectorial term, then the first subdiagonal,
      // then the three-term recurrence P(n,m) = a t P(n-1,m) - b P(n-2,m).
      col[0] = sect;
      if (m < nmax)
        col[1] = ca[1] * t * sect;
      const int len = nmax - m;
      for (int j = 2; j <= len; ++j)
        col[j] = ca[j] * t * col[j - 1] - cb[j] * col[j - 2];
    }
  }

  void msfc(double t, double u, const Coefs &c, std::span<double> P)
  {
    check(c, P);
    if (pole_shortcut(t, u, c, P))
      return;
    const int nmax = c.nmax;
    // Sectorial mantissa xs and its exponent ks (value = xs * R^-ks). ks is shared by the whole sectorial chain.
    double xs = 1.0;
    int ks = 0;
    // Reconstructs x * R^-k as a plain double: k = 0 is exact, k >= 3 is treated as underflow to 0.
    const auto rebuild = [](double x, int k)
    {
      if (k == 0)
        return x;
      if (k >= 3)
        return 0.0;
      return std::ldexp(x, -RBITS * k);
    };
    for (int m = 0; m <= nmax; ++m)
    {
      if (m > 0)
      {
        // Sectorial step xs <- d_m u xs; when it falls below R^-1, scale it up by R and count one exponent step.
        xs *= c.d[m] * u;
        while (xs != 0.0 && std::abs(xs) < R_INV)
        {
          xs *= R_BIG;
          ++ks;
        }
      }
      const auto off = col_offset(m, nmax);
      const double *ca = c.a.data() + off, *cb = c.b.data() + off;
      double *col = P.data() + off;
      // Column recurrence on mantissas: same coefficients and same exponent k for the whole column.
      int k = ks;
      double y2 = xs;
      col[0] = rebuild(y2, k);
      if (m == nmax)
        break;
      double y1 = ca[1] * t * y2;
      col[1] = rebuild(y1, k);
      const int len = nmax - m;
      for (int j = 2; j <= len; ++j)
      {
        double y = ca[j] * t * y1 - cb[j] * y2;
        // Keep the mantissa bounded: if the new term exceeds R, divide the last two terms by R and decrement k.
        // The current value is unchanged because y, y1 and k all move together.
        if (k > 0 && std::abs(y) > R_BIG)
        {
          y *= R_INV;
          y1 *= R_INV;
          --k;
        }
        col[j] = rebuild(y, k);
        y2 = y1;
        y1 = y;
      }
    }
  }

  namespace
  {
    // Extended-range number: value = m * 2^e with 0.5 <= |m| < 1 (m = 0 represents zero).
    // The exponent e is a 64-bit integer, so it cannot overflow for any practical degree.
    struct XNum
    {
      double m = 0.0;
      long long e = 0;
    };
    // Renormalises a mantissa to [0.5, 1) and folds the binary exponent into e.
    inline XNum xnorm(double m, long long e)
    {
      if (m == 0.0)
        return {0.0, 0};
      int k;
      m = std::frexp(m, &k);
      return {m, e + k};
    }
    inline XNum xmul(const XNum &a, double s) { return xnorm(a.m * s, a.e); }
    // Returns a*s1 - b*s2, aligning exponents first. An addend more than 2^1100 times smaller than the other
    // one is dropped, because it lies below double-precision resolution of the sum.
    inline XNum xfma2(const XNum &a, double s1, const XNum &b, double s2)
    {
      const XNum x = xmul(a, s1), y = xmul(b, -s2);
      if (x.m == 0.0)
        return y;
      if (y.m == 0.0)
        return x;
      const XNum &hi = (x.e >= y.e) ? x : y;
      const XNum &lo = (x.e >= y.e) ? y : x;
      const long long d = hi.e - lo.e;
      return xnorm(hi.m + (d > 1100 ? 0.0 : std::ldexp(lo.m, -static_cast<int>(d))), hi.e);
    }
    // Converts back to a double; values below 2^-1200 are returned as 0 (they are far below the IEEE range).
    inline double xvalue(const XNum &x)
    {
      if (x.m == 0.0 || x.e < -1200)
        return 0.0;
      return std::ldexp(x.m, static_cast<int>(x.e));
    }
  }

  void era(double t, double u, const Coefs &c, std::span<double> P)
  {
    check(c, P);
    if (pole_shortcut(t, u, c, P))
      return;
    const int nmax = c.nmax;
    // Sectorial seeds: seed[m] = P(m,m) in extended range, built serially with the same d_m u product.
    std::vector<XNum> seed(nmax + 1);
    seed[0] = {0.5, 1};
    for (int m = 1; m <= nmax; ++m)
      seed[m] = xmul(seed[m - 1], c.d[m] * u);

    // Columns are independent once their seeds exist. Dynamic scheduling because column m has nmax - m + 1 entries.
#pragma omp parallel for schedule(dynamic, 8) if (nmax > 64)
    for (int m = 0; m <= nmax; ++m)
    {
      const auto off = col_offset(m, nmax);
      const double *ca = c.a.data() + off, *cb = c.b.data() + off;
      double *col = P.data() + off;
      XNum y2 = seed[m];
      col[0] = xvalue(y2);
      if (m == nmax)
        continue;
      XNum y1 = xmul(y2, ca[1] * t);
      col[1] = xvalue(y1);
      const int len = nmax - m;
      for (int j = 2; j <= len; ++j)
      {
        const XNum y = xfma2(y1, ca[j] * t, y2, cb[j]);
        col[j] = xvalue(y);
        y2 = y1;
        y1 = y;
      }
    }
  }

  void compute(double t, double u, const Coefs &c, std::span<double> P)
  {
    // 1800 is the SFC validity limit measured in this project (benchmark/alf); Holmes and Featherstone (2002) discuss
    // the underflow of the sectorial term that causes it.
    if (c.nmax <= 1800)
      sfc(t, u, c, P);
    else
      msfc(t, u, c, P);
  }

  void dphi(int nmax, std::span<const double> P, std::span<double> dP)
  {
    if (static_cast<std::int64_t>(P.size()) < tri_size(nmax) || static_cast<std::int64_t>(dP.size()) < tri_size(nmax))
      throw std::invalid_argument("dphi: buffer insuficiente");
    // Each (n, m) entry reads only the same degree n (neighbouring orders), so the loop over m is parallel.
#pragma omp parallel for schedule(dynamic, 16) if (nmax > 256)
    for (int m = 0; m <= nmax; ++m)
    {
      const double M = m;
      // Factor sqrt(1 + δ_m1) of the identity: sqrt(2) for m = 1, 1 otherwise.
      const double sqrt_m1 = (m == 1) ? std::sqrt(2.0) : 1.0;
      for (int n = m; n <= nmax; ++n)
      {
        const double N = n;
        double v;
        if (m == 0)
          // Zonal case: dP(n,0)/dφ = sqrt(n(n+1)/2) P(n,1); degree 0 has zero derivative.
          v = (n >= 1) ? std::sqrt(N * (N + 1.0) / 2.0) * P[idx_tri(n, 1, nmax)] : 0.0;
        else
        {
          // Non-zonal case: 1/2 [ -sqrt((n+m)(n-m+1)) sqrt(1+δ_m1) P(n,m-1) + sqrt((n-m)(n+m+1)) P(n,m+1) ].
          v = -std::sqrt((N + M) * (N - M + 1.0)) * sqrt_m1 * P[idx_tri(n, m - 1, nmax)];
          if (m < n)
            v += std::sqrt((N - M) * (N + M + 1.0)) * P[idx_tri(n, m + 1, nmax)];
          v *= 0.5;
        }
        dP[idx_tri(n, m, nmax)] = v;
      }
    }
  }

  namespace
  {
    // Team size used when the caller passes threads = 0.
    int default_team()
    {
#ifdef _OPENMP
      return std::max(1, omp_get_max_threads());
#else
      return 1;
#endif
    }
  }

  Extracted extract(std::span<const double> lat_rad, const Coefs &c, const Selection &sel, int threads)
  {
    if (c.nmax < 0)
      throw std::invalid_argument("extract: coefs are empty (call make_coefs first)");
    const int cn = c.nmax;
    const int snmax = sel.nmax < 0 ? cn : sel.nmax;
    if (sel.nmin < 0 || sel.nmin > snmax || snmax > cn)
      throw std::invalid_argument("extract: require 0 <= nmin <= nmax <= coefs.nmax (nmin=" + std::to_string(sel.nmin) +
                                  ", nmax=" + std::to_string(snmax) + ", coefs.nmax=" + std::to_string(cn) + ")");
    if (threads < 0)
      throw std::invalid_argument("extract: threads must be >= 0 (0 = OpenMP default)");

    // Orders: empty = all 0..nmax. Given orders are sorted, so the output order does not depend on the input order.
    std::vector<int> orders = sel.orders;
    if (orders.empty())
    {
      orders.resize(static_cast<std::size_t>(snmax) + 1);
      std::iota(orders.begin(), orders.end(), 0);
    }
    else
    {
      std::sort(orders.begin(), orders.end());
      for (const int m : orders)
        if (m < 0 || m > snmax)
          throw std::invalid_argument("extract: order m=" + std::to_string(m) + " outside 0..nmax=" + std::to_string(snmax));
      if (const auto dup = std::adjacent_find(orders.begin(), orders.end()); dup != orders.end())
        throw std::invalid_argument("extract: duplicate order m=" + std::to_string(*dup));
    }

    Extracted out;
    out.count = 0;
    for (const int m : orders)
      out.count += snmax - std::max(m, sel.nmin) + 1;
    const std::int64_t count = out.count;

    // Output size check, before any allocation. Overflow-safe: np * count * per_point > limit is tested as a division.
    const auto np = static_cast<std::ptrdiff_t>(lat_rad.size());
    const std::int64_t per_point = sel.derivative ? 2 : 1;
    if (np > 0 && count > 0 && static_cast<std::int64_t>(np) > extract_max_doubles / (count * per_point))
      throw std::length_error("extract: output of " + std::to_string(np) + " points x " + std::to_string(count) +
                              " columns x " + std::to_string(per_point) + " array(s) exceeds the limit of 2^31 doubles (16 GB). " +
                              "Request fewer points or orders, a lower nmax, or derivative=false (memory = 8 * npoints * count bytes).");

    // Column table: one entry per stored element, plus the packed index into the full triangle of degree cn.
    std::vector<std::int64_t> src;
    src.reserve(static_cast<std::size_t>(count));
    out.n.reserve(static_cast<std::size_t>(count));
    out.m.reserve(static_cast<std::size_t>(count));
    for (const int m : orders)
      for (int n = std::max(m, sel.nmin); n <= snmax; ++n)
      {
        out.n.push_back(n);
        out.m.push_back(m);
        src.push_back(idx_tri(n, m, cn));
      }

    const auto total = static_cast<std::size_t>(np) * static_cast<std::size_t>(count);
    out.P.assign(total, 0.0);
    if (sel.derivative)
      out.dP.assign(total, 0.0);

    const int team = threads > 0 ? threads : default_team();
    const bool want_dp = sel.derivative;
    double *const outP = out.P.data();
    double *const outDP = out.dP.data();
    const auto ucount = static_cast<std::size_t>(count);
    const std::size_t nsrc = src.size();
    const std::int64_t *const srcp = src.data();

    // One buffer per thread (allocated once per call, reused for every point of that thread). Each point writes
    // only its own row of the output, so there are no races and the result does not depend on the thread count.
#pragma omp parallel num_threads(team)
    {
      AlignedVec<double> P(static_cast<std::size_t>(tri_size(cn)));
      AlignedVec<double> dP(want_dp ? static_cast<std::size_t>(tri_size(cn)) : 0);
#pragma omp for schedule(dynamic, 4)
      for (std::ptrdiff_t i = 0; i < np; ++i)
      {
        const double phi = lat_rad[static_cast<std::size_t>(i)];
        compute(std::sin(phi), std::cos(phi), c, std::span<double>(P.data(), P.size()));
        const auto row = static_cast<std::size_t>(i) * ucount;
        for (std::size_t k = 0; k < nsrc; ++k)
          outP[row + k] = P[static_cast<std::size_t>(srcp[k])];
        if (want_dp)
        {
          dphi(cn, std::span<const double>(P.data(), P.size()), std::span<double>(dP.data(), dP.size()));
          for (std::size_t k = 0; k < nsrc; ++k)
            outDP[row + k] = dP[static_cast<std::size_t>(srcp[k])];
        }
      }
    }
    return out;
  }
}
