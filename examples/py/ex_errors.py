"""ex_errors.py: commission and omission errors of the height anomaly and of the undulation, and the error report.

Python equivalent of examples/cpp/ex_errors.cpp. Same arguments, same printed quantities.

Shows: the commission error from the formal standard deviations of the coefficients (load_gfc_sigma,
commission_error, commission_global), the omission error with Kaula's rule (omission_kaula_undulation) and with the
Tscherning-Rapp degree-variance model (omission_tscherning_rapp_undulation), and the consolidated report
(make_error_report, format_error_report). EGM2008 is a calibrated model (its file header says "errors calibrated"),
so its sigma columns are realistic errors (Pavlis et al. 2012). Numerical error needs reference values, so the
report shows it as not available here.

Usage: python examples/py/ex_errors.py <model.gfc> [nmax=2190]
  nmax: truncation degree, >= 2 and not above the degree of the file. Lower values give quick runs.

Output: commission error of zeta_ell at three points [m], the global commission of N, the Kaula and Tscherning-Rapp
omission errors of N [m], and the full report. All global values are one standard deviation.
"""

import argparse
import sys

import gshtools as gsh

GAMMA = 9.81  # mean normal gravity [m/s^2], the value used by the report
TR_N_TAIL = 10800  # upper degree of the Tscherning-Rapp sums


def main():
    parser = argparse.ArgumentParser(description="Commission and omission errors of N and zeta (gshtools).")
    parser.add_argument("model", help="fully normalised ICGEM .gfc model with formal errors (sigma columns)")
    parser.add_argument("nmax", nargs="?", type=int, default=2190, help="truncation degree, >= 2 (default 2190)")
    args = parser.parse_args()
    if args.nmax < 2:
        parser.error("nmax must be at least 2")
    nmax = args.nmax
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    ell = gsh.Ellipsoid.GRS80()
    sigma = gsh.load_gfc_sigma(args.model, nmax)  # C and S hold the formal standard deviations
    coef = gsh.load_gfc(args.model, nmax)         # the coefficients themselves (for the TR tail)

    # Three points on the ellipsoid (height 0): equator, mid-latitude and a high-latitude point.
    lat_geodetic = [0.0, 45.0, 80.0]
    pts = gsh.Points.from_geodetic(lat_geodetic, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ell.a, ell.e2)
    opts = gsh.ErrorOptions()
    opts.nmin = 2
    opts.nmax = nmax
    ce = gsh.commission_error(sigma, pts, gsh.ErrorFunctional.ZetaEll, opts)

    print(f"Commission error of zeta_ell (degrees 2..{nmax}, one sigma):")
    for lat, value in zip(lat_geodetic, ce):
        print(f"  lat {lat:5.1f} deg (geodetic, h = 0): {value:.6f} m")

    ce_global_N = gsh.commission_global(sigma, gsh.ErrorFunctional.ZetaEll, sigma.R, GAMMA, 2, nmax)
    print(f"Global commission error CE(N), degrees 2..{nmax}: {ce_global_N:.6f} m")

    oe_kaula = gsh.omission_kaula_undulation(sigma.R, nmax)
    oe_tr = gsh.omission_tscherning_rapp_undulation(nmax, TR_N_TAIL, GAMMA)
    print(f"Omission error of N above degree {nmax}: Kaula {oe_kaula:.6f} m, "
          f"Tscherning-Rapp (degrees to {TR_N_TAIL}) {oe_tr:.6f} m")

    opt = gsh.ReportOptions()  # the Tscherning-Rapp value uses the model's own variances above nmax (coefficients=)
    report = gsh.make_error_report(sigma, nmax, sigma.R, GAMMA, [], pts, opt, coef)
    print()
    print(gsh.format_error_report(report), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
