// errors.hpp: commission and omission errors of spherical-harmonic functionals.
//
// Commission and omission errors: Tscherning and Rapp (1974) for the omission model; the commission weights follow
// Bucha and Janák (2019).
//
// All functionals handled here are linear in the coefficients, F = sum_n w_n sum_m Pbar_nm(sin φ) (C cos mλ + S sin mλ),
// with a degree weight w_n that depends on the functional (see ErrorFunctional). They use the spherical
// approximation of the references (no ellipsoidal corrections), which is the usual one for formal errors.
//
//   Functional            w_n                                       Unit
//   T                     GM/r (R/r)^n                              m^2/s^2
//   ZetaEll (T/γ)         GM/(r γ) (R/r)^n                          m
//   GravityAnomaly        GM/r^2 (n-1) (R/r)^n                      m/s^2   (spherical Δg, same as Functional::Dg_sa)
//   GravityDisturbance    GM/r^2 (n+1) (R/r)^n                      m/s^2   (radial disturbance -dT/dr)
//
// Commission error: the formal error of F caused by the standard deviations of the coefficients. The .gfc
// format only carries sigma C and sigma S, so the covariance matrix of the coefficients is taken as diagonal
// (correlations between coefficients are ignored). Then
//   sigma_F^2 = sum_n w_n^2 sum_m Pbar_nm^2 (sigma_C,nm^2 cos^2 mλ + sigma_S,nm^2 sin^2 mλ).
// Whether the result is an absolute error depends on the sigmas of the model: "calibrated" models (EGM2008) give
// realistic errors, "formal" ones give only relative geographic precision (Bucha and Janák 2019).
//
// Omission error: the signal of the degrees above the truncation degree nmax. Its global mean square follows from
// the degree variances sigma_n^2 = sum_m (C_nm^2 + S_nm^2) (Bucha and Janák 2019):
//   e^2 = sum_{n>nmax} w_n^2 sigma_n^2.
// This is a global (area-mean) value, not a point value.
#pragma once
#include <vector>
#include "core.hpp"
#include "ellipsoid.hpp"
#include "model.hpp"
#include "synth_types.hpp"

namespace gsh
{
  enum class ErrorFunctional
  {
    T,
    ZetaEll,
    GravityAnomaly,
    GravityDisturbance
  };
  // Name of the functional, e.g. "T".
  const char *error_functional_name(ErrorFunctional f);

  struct ErrorOptions
  {
    int nmin = 0;                             // lowest degree included
    int nmax = -1;                            // highest degree included; < 0 means the degree of the model
    Ellipsoid ellipsoid = Ellipsoid::GRS80(); // supplies the normal gravity γ for ZetaEll
    int threads = 0;                          // OpenMP threads; 0 = the OpenMP maximum
  };

  // Commission error at every point (one standard deviation, in the unit of the table above).
  //   sigma: model whose C and S hold the standard deviations (load_gfc_sigma). Its GM and R scale the result.
  //   pts: geocentric latitude and longitude in radians and radius in metres (Points).
  // Cost: one Legendre triangle per point, O(nmax^2) each; parallel over points. Handles any degree supported by
  // legendre::compute (including 2190 and above).
  // Throws std::invalid_argument if opts.nmax > sigma.nmax or nmin is outside [0, nmax].
  std::vector<double> commission_error(const Model &sigma, const Points &pts, ErrorFunctional f, const ErrorOptions &opts = {});

  // Commission error per degree: for each n in [nmin, nmax], the area-mean contribution
  // sqrt(w_n^2 sum_m (sigma_C^2 + sigma_S^2)) at radius r (an unweighted global value, no Legendre functions needed).
  // The vector has nmax + 1 entries; degrees outside [nmin, nmax] are 0. gamma is only used by ZetaEll.
  std::vector<double> commission_per_degree(const Model &sigma, ErrorFunctional f, double r, double gamma, int nmin, int nmax);

  // Degree variances sigma_n^2 = sum_{m=0}^{n} (C_nm^2 + S_nm^2) of a coefficient model, for n = 0..model.nmax.
  std::vector<double> degree_variances(const Model &model);

  // Parameters of the Tscherning-Rapp anomaly degree-variance model (Tscherning and Rapp 1974, OSU Report 208, Table 7):
  //   c_n(Δg) = A (n-1) / ((n-2)(n+B)) s^(n+2)   [mGal^2], n >= 3, s = R_B^2 / R_E^2.
  // Defaults: A = 425.28 mGal^2, B = 24, s = 0.999617 with R_E = 6371 km (mean Earth radius), which gives the
  // Bjerhammar sphere radius R_B = R_E sqrt(s) = 6369.78 km (note: not 6370 km, which would give s = 0.999686).
  struct TscherningRapp
  {
    double A = 425.28;     // mGal^2
    double B = 24.0;
    double R_B = 6369.78e3; // m, R_E sqrt(0.999617)
    double R_E = 6371.0e3; // m
    double s() const { return (R_B * R_B) / (R_E * R_E); }
  };

  // Global (area-mean) commission error of the model: sqrt( sum_{n=nmin}^{nmax} w_n^2 sum_m (sigma_C^2 + sigma_S^2) ),
  // with the weights of the table above at radius r. For T and nmin = 2 this is
  //   CE(T) = (GM/a) sqrt( sum_n (a/r)^(2n+2) sum_m (sigma_C^2 + sigma_S^2) );  CE(N) = CE(T)/γ.
  // gamma is only used by ZetaEll. It is the quadrature sum of commission_per_degree().
  double commission_global(const Model &sigma, ErrorFunctional f, double r, double gamma, int nmin, int nmax);

  // Kaula-rule estimate of the omission error of the undulation N (global, in metres):
  //   OE_Kaula(N) ~ R 1e-5 / (nmax + 1).
  double omission_kaula_undulation(double R, int nmax);

  // Tscherning-Rapp omission error of the undulation N (global, in metres):
  //   OE^2(N) = (R_E/γ)^2 sum_{n=nmax+1}^{n_tail} c_n(Δg) / (n-1)^2,   c_n in (m/s^2)^2 after converting from mGal^2.
  // Typical n_tail is about 10800. gamma is a mean normal gravity (about 9.81 m/s^2).
  // Throws std::invalid_argument if nmax < 2 or n_tail <= nmax.
  double omission_tscherning_rapp_undulation(int nmax, int n_tail, double gamma, const TscherningRapp &p = {});

  // Relative error |calc - ref| / |ref| (Bucha and Janák 2019). Returns |calc - ref| if ref is 0.
  double relative_error(double calc, double ref);

  // Upper bound of the decimal digits lost by the Legendre recurrence at degree M and colatitude theta [rad]
  // (Fukushima 2012): log10(M) + log10(10 + 2 |theta cot theta|). At theta = 0 the limit
  // theta cot theta = 1 is used.
  double digit_loss_bound(int M, double theta);

  // Coefficient degree variances sigma_n^2 (same meaning as degree_variances) implied by the TR model for the
  // degrees n_first..n_last, at the model sphere of radius R: sigma_n^2 = c_n / ((n-1) GM/R^2)^2, with c_n
  // converted from mGal^2 to (m/s^2)^2. The result has n_last + 1 entries; entries below n_first are 0.
  // Throws std::invalid_argument if n_first < 3 or n_last < n_first.
  std::vector<double> tscherning_rapp_variances(double GM, double R, int n_first, int n_last, const TscherningRapp &p = {});

  // Degree variances of `model` up to model.nmax, extended with the Tscherning-Rapp model up to n_end, which is
  // the usual way to bound the omission error of a model beyond its own maximum degree.
  std::vector<double> degree_variances_with_tail(const Model &model, int n_end, const TscherningRapp &p = {});

  // Omission error: sqrt(sum_{n=nmax+1}^{N} w_n^2 sigma_n^2), where N = sigma2.size() - 1 and sigma2[n] is the
  // degree variance at the model sphere of radius R. r is the evaluation radius (upper-degree attenuation
  // (R/r)^n is included) and gamma is only used by ZetaEll. Returns 0 when nmax >= N.
  // Throws std::invalid_argument if nmax < 0.
  double omission_error(const std::vector<double> &sigma2, int nmax, ErrorFunctional f, double GM, double R, double r, double gamma);
}
