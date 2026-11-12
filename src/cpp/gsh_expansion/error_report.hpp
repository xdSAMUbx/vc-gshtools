// error_report.hpp: consolidated error report of a run (numerical, omission and commission errors).
//
// Builds on errors.hpp. Everything here is for the undulation N N (metres) unless stated otherwise:
//   - Numerical error (NE): max relative error of (calc, ref) pairs, or a precomputed E_rel, with the digits
//     retained and the digit_loss_bound of the Legendre recurrence at the run degree and a colatitude.
//   - Omission error (OE): Kaula rule and Tscherning-Rapp (default degree variances) for degrees above nmax up
//     to n_tail. Optionally the Tscherning-Rapp value with the degree variances of a coefficient model.
//   - Commission error (CE): global CE(N) of the sigma model (commission_global, ZetaEll) and optional point
//     values (commission_error, ZetaEll) with min/mean/max.
#pragma once
#include <cstddef>
#include <string>
#include <utility>
#include <vector>
#include "errors.hpp"
#include "model.hpp"
#include "synth_types.hpp"

namespace gsh
{
  // Numerical error. The digit-loss bound is always filled; the E_rel fields only when available is true.
  struct NumericalError
  {
    bool available = false;                    // true when E_rel was given (pairs or precomputed)
    double max_relative_error = 0.0;           // E_rel = max |calc - ref| / |ref|; NaN when not available
    bool within_1e12 = false;                  // E_rel <= 1e-12
    double digits = 0.0;                       // -log10(E_rel), capped at 16; NaN when not available
    double colatitude = 0.0;                   // colatitude [rad] used for the digit-loss bound
    double digit_loss_bound = 0.0;             // digit_loss_bound(nmax, colatitude) [decimal digits]
  };

  // Omission error of the undulation N (global, metres), for degrees nmax+1 .. n_tail.
  struct OmissionError
  {
    int nmax = 0;                              // run truncation degree
    int n_tail = 0;                            // upper degree of the Tscherning-Rapp sums
    double kaula = 0.0;                        // R 1e-5 / (nmax+1), R = sigma model radius
    double tr = 0.0;                           // Tscherning-Rapp with the default TR degree variances
    bool has_tr_model = false;                 // true when tr_model was computed
    double tr_model = 0.0;                     // TR tail plus the coefficient model's own variances above nmax
    int coefficient_nmax = 0;                  // degree of the coefficient model used for tr_model
  };

  // Commission error of the undulation N (metres).
  struct CommissionError
  {
    int nmax = 0;                              // degrees used: min(run nmax, sigma model nmax)
    double global_N = 0.0;                     // CE(N) = commission_global(ZetaEll, degrees 2..nmax)
    bool has_points = false;                   // true when point values were computed
    std::size_t count = 0;                     // number of points
    double min = 0.0, mean = 0.0, max = 0.0;   // statistics of the point commission errors
  };

  struct ErrorReport
  {
    std::string model_name;                    // header name of the sigma model (may be empty)
    double GM = 0.0;                           // m^3/s^2
    double R = 0.0;                            // m
    int model_nmax = 0;                        // degree of the sigma model
    int nmax = 0;                              // run truncation degree
    double r = 0.0;                            // m, radius of the global commission value
    double gamma = 0.0;                        // m/s^2, mean normal gravity used by the global values
    NumericalError numerical;
    OmissionError omission;
    CommissionError commission;
  };

  // Optional settings of make_error_report.
  struct ReportOptions
  {
    double colatitude = PI / 2;                // colatitude [rad] for digit_loss_bound (default: equator)
    int n_tail = 10800;                        // upper degree of the Tscherning-Rapp sums (> nmax)
    const Model *coefficients = nullptr;       // coefficient model (load_gfc, not the sigma model); enables tr_model
    TscherningRapp tr = {};                    // Tscherning-Rapp parameters
  };

  // Builds the report of a run.
  //   sigma_model: model read with load_gfc_sigma (its C and S hold the sigmas). Its GM, R and name are reported.
  //   nmax: run truncation degree (>= 2). Values above sigma_model.nmax are clipped for the commission part.
  //   r: radius [m] of the global commission value. Undulation omission values are on the model sphere R.
  //   gamma: mean normal gravity [m/s^2] for the global values (9.81 is typical).
  //   calc_ref: optional (calc, ref) pairs for the numerical error. Empty means not computed.
  //   pts: optional points for the point commission statistics (ZetaEll with per-point normal gravity).
  //   opt: colatitude, n_tail, coefficient model and TR parameters.
  // Throws std::invalid_argument if nmax < 2, n_tail <= nmax, or sigma_model.nmax < 2.
  ErrorReport make_error_report(const Model &sigma_model, int nmax, double r, double gamma = 9.81,
                                const std::vector<std::pair<double, double>> &calc_ref = {}, const Points *pts = nullptr,
                                const ReportOptions &opt = {});

  // Sets the numerical part from a precomputed max relative error (digits, within_1e12 and available follow).
  void set_relative_error(ErrorReport &report, double max_relative_error);

  // Compact human-readable text, several lines, units in the labels.
  std::string format_error_report(const ErrorReport &report);

  // Single-line JSON object. Values that are not available are written as null.
  std::string error_report_json(const ErrorReport &report);
}
