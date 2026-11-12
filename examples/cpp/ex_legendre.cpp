// ex_legendre.cpp: fully normalised associated Legendre functions (legendre.hpp) at one latitude.
//
// Shows: the three recurrences sfc(), msfc() and era() and how they agree; the dispatcher compute(); the sum-rule
// identity sum_{m=0}^{n} Pbar_nm^2 = 2n+1; and the latitude derivative dphi(), checked against a central finite
// difference of sfc(). No model file is needed.
//
// Usage: ex_legendre [lat_deg=45] [nmax=200]
//   lat_deg: GEOCENTRIC (spherical) latitude in degrees, -90..90.
//   nmax:    maximum degree, 2..1800. This is the range where sfc() is valid, so all recurrences can be compared.
//
// Output (all values dimensionless, fully normalised with 4 pi, no Condon-Shortley phase):
//   max |a - b| is the largest absolute difference and max rel is the largest relative difference over the
//   entries where the reference is not below 1e-100 (msfc and era differ from sfc only by rounding).
//   The sum rule should give a deviation near 1e-13; it holds at every latitude.
//   dP/dphi is the derivative with respect to latitude in radians; the finite difference uses a step of 1e-6 rad.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "legendre.hpp"

namespace
{
  constexpr int kMaxSfcDegree = 1800; // sfc() is accurate up to this degree (legendre.hpp)

  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "ex_legendre: %s\nusage: ex_legendre [lat_deg=45] [nmax=200] (nmax in 2..%d)\n", why, kMaxSfcDegree);
    std::exit(2);
  }

  double parse_double(const char *s)
  {
    char *end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s || *end != '\0')
      usage("latitude must be a number");
    return v;
  }

  int parse_int(const char *s)
  {
    char *end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0')
      usage("nmax must be an integer");
    return static_cast<int>(v);
  }

  struct Diff
  {
    double abs = 0.0, rel = 0.0;
  };

  // Largest absolute and relative difference between a and the reference ref.
  Diff compare(const std::vector<double> &a, const std::vector<double> &ref)
  {
    Diff d;
    for (std::size_t i = 0; i < ref.size(); ++i)
    {
      const double e = std::fabs(a[i] - ref[i]);
      d.abs = std::max(d.abs, e);
      if (std::fabs(ref[i]) > 1e-100)
        d.rel = std::max(d.rel, e / std::fabs(ref[i]));
    }
    return d;
  }
} // namespace

int main(int argc, char **argv)
{
  if (argc > 3)
    usage("too many arguments");
  const double lat_deg = argc > 1 ? parse_double(argv[1]) : 45.0;
  const int nmax = argc > 2 ? parse_int(argv[2]) : 200;
  if (lat_deg < -90.0 || lat_deg > 90.0)
    usage("latitude must be in [-90, 90] degrees");
  if (nmax < 2 || nmax > kMaxSfcDegree)
    usage("nmax must be in 2..1800");

  const double phi = gsh::deg2rad(lat_deg);
  const double t = std::sin(phi), u = std::cos(phi); // t = sin(phi), u = cos(phi)
  const gsh::legendre::Coefs c = gsh::legendre::make_coefs(nmax);
  const std::size_t size = static_cast<std::size_t>(gsh::tri_size(nmax));

  std::vector<double> P_sfc(size), P_msfc(size), P_era(size), P_auto(size);
  gsh::legendre::sfc(t, u, c, P_sfc);
  gsh::legendre::msfc(t, u, c, P_msfc);
  gsh::legendre::era(t, u, c, P_era);
  gsh::legendre::compute(t, u, c, P_auto);

  const Diff d_msfc = compare(P_msfc, P_sfc);
  const Diff d_era = compare(P_era, P_sfc);
  const Diff d_auto = compare(P_auto, P_sfc);

  // Sum rule: for every n, sum_m Pbar_nm^2 = 2n+1 on the unit sphere (t^2 + u^2 = 1).
  double worst_sum = 0.0;
  for (int n = 0; n <= nmax; ++n)
  {
    double s = 0.0;
    for (int m = 0; m <= n; ++m)
    {
      const double p = P_sfc[gsh::idx_tri(n, m, nmax)];
      s += p * p;
    }
    worst_sum = std::max(worst_sum, std::fabs(s / (2 * n + 1) - 1.0));
  }

  // Latitude derivative (analytic, dphi) against a central finite difference of sfc().
  std::vector<double> dP(size), P_plus(size), P_minus(size);
  gsh::legendre::dphi(nmax, P_sfc, dP);
  const double h = 1e-6;
  gsh::legendre::sfc(std::sin(phi + h), std::cos(phi + h), c, P_plus);
  gsh::legendre::sfc(std::sin(phi - h), std::cos(phi - h), c, P_minus);
  double worst_fd = 0.0;
  for (std::size_t i = 0; i < size; ++i)
    worst_fd = std::max(worst_fd, std::fabs(dP[i] - (P_plus[i] - P_minus[i]) / (2.0 * h)));

  const int n_show = nmax, m_show = 2;
  const std::size_t i_show = static_cast<std::size_t>(gsh::idx_tri(n_show, m_show, nmax));

  std::printf("ex_legendre: lat=%.6f deg (geocentric), nmax=%d, fully normalised (4 pi)\n", lat_deg, nmax);
  std::printf("  sfc vs msfc:      max |a-b| = %.3e, max rel = %.3e\n", d_msfc.abs, d_msfc.rel);
  std::printf("  sfc vs era:       max |a-b| = %.3e, max rel = %.3e\n", d_era.abs, d_era.rel);
  std::printf("  compute vs sfc:   max |a-b| = %.3e (compute uses sfc up to degree %d)\n", d_auto.abs, kMaxSfcDegree);
  std::printf("  sum rule:         max |sum_m Pbar_nm^2/(2n+1) - 1| = %.3e over n = 0..%d\n", worst_sum, nmax);
  std::printf("  (n,m) = (%d,%d):  P = %.12e, dP/dphi = %.12e, finite difference = %.12e\n", n_show, m_show,
              P_sfc[i_show], dP[i_show], (P_plus[i_show] - P_minus[i_show]) / (2.0 * h));
  std::printf("  dP/dphi vs finite difference over the whole triangle: max |diff| = %.3e\n", worst_fd);
  return 0;
}
