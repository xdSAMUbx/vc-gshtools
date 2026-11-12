// ex_errors.cpp: commission and omission errors of the height anomaly and of the undulation, and the error report.
//
// Shows: the commission error from the formal standard deviations of the coefficients (load_gfc_sigma,
// commission_error, commission_global), the omission error with Kaula's rule (omission_kaula_undulation) and with the
// Tscherning-Rapp degree-variance model (omission_tscherning_rapp_undulation), and the consolidated report
// (make_error_report, format_error_report). EGM2008 is a calibrated model (its file header says "errors calibrated"),
// so its sigma columns are realistic errors (Pavlis et al. 2012). Numerical error needs reference values, so the
// report shows it as not available here.
//
// Usage: ex_errors <model.gfc> [nmax=2190]
//   nmax: truncation degree, >= 2 and not above the degree of the file. Lower values give quick runs.
//
// Output: commission error of zeta_ell at three points [m], the global commission of N, the Kaula and Tscherning-Rapp
// omission errors of N [m], and the full report. All global values are one standard deviation.

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>
#include "ellipsoid.hpp"
#include "error_report.hpp"
#include "errors.hpp"
#include "model.hpp"
#include "synth_types.hpp"

namespace
{
  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "ex_errors: %s\nusage: ex_errors <model.gfc> [nmax=2190]\n", why);
    std::exit(2);
  }

  int parse_int(const char *s)
  {
    char *end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0')
      usage("nmax must be an integer");
    return static_cast<int>(v);
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    if (argc < 2 || argc > 3)
      usage("expected 1 or 2 arguments");
    const int nmax = argc > 2 ? parse_int(argv[2]) : 2190;
    if (nmax < 2)
      usage("nmax must be at least 2");

    const gsh::Ellipsoid ell = gsh::Ellipsoid::GRS80();
    const double gamma = 9.81; // mean normal gravity [m/s^2], the value used by the report
    const gsh::Model sigma = gsh::load_gfc_sigma(argv[1], nmax); // C and S hold the formal standard deviations
    const gsh::Model coef = gsh::load_gfc(argv[1], nmax);        // the coefficients themselves (for the TR tail)

    // Three points on the ellipsoid (height 0): equator, mid-latitude and a high-latitude point.
    const gsh::Points pts = gsh::Points::from_geodetic({0.0, 45.0, 80.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, ell.a, ell.e2);
    gsh::ErrorOptions opts;
    opts.nmin = 2;
    opts.nmax = nmax;
    const std::vector<double> ce = gsh::commission_error(sigma, pts, gsh::ErrorFunctional::ZetaEll, opts);

    std::printf("Commission error of zeta_ell (degrees 2..%d, one sigma):\n", nmax);
    const double lat[] = {0.0, 45.0, 80.0};
    for (std::size_t i = 0; i < ce.size(); ++i)
      std::printf("  lat %5.1f deg (geodetic, h = 0): %.6f m\n", lat[i], ce[i]);

    const double ce_global_N = gsh::commission_global(sigma, gsh::ErrorFunctional::ZetaEll, sigma.R, gamma, 2, nmax);
    std::printf("Global commission error CE(N), degrees 2..%d: %.6f m\n", nmax, ce_global_N);

    const double oe_kaula = gsh::omission_kaula_undulation(sigma.R, nmax);
    const double oe_tr = gsh::omission_tscherning_rapp_undulation(nmax, 10800, gamma);
    std::printf("Omission error of N above degree %d: Kaula %.6f m, Tscherning-Rapp (degrees to 10800) %.6f m\n", nmax,
                oe_kaula, oe_tr);

    gsh::ReportOptions opt;
    opt.coefficients = &coef; // enables the TR value with the model's own variances above nmax
    const gsh::ErrorReport report = gsh::make_error_report(sigma, nmax, sigma.R, gamma, {}, &pts, opt);
    std::printf("\n%s", gsh::format_error_report(report).c_str());
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "ex_errors: %s\n", e.what());
    return 1;
  }
}
