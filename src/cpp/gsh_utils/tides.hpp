// tides.hpp: permanent-tide systems (tide-free NT / zero-tide ZT / mean-tide MT): conversion of the
// geopotential coefficient C̄20, of the geoid height N, of the ellipsoidal height h, of the gravity
// value g, and of the separation surface S used by the IHRF.
//
// Conventions:
//   - Geodetic latitude in degrees; heights in metres; gravity in mGal (1 mGal = 1e-5 m/s^2);
//     C̄20 fully normalized and dimensionless.
//   - Every shift is defined as "value in the target system minus value in the source system".
//   - Vectorized overloads (std::vector<double>) run element by element with OpenMP.
//
// Sources (see the comment of each function in tides.cpp for the formula it implements):
//   - Tide systems, IHRF conventions and conversion formulas: Makinen, J. (2021), "The permanent tide and the
//     International Height Reference Frame IHRF", Journal of Geodesy 95, 106 (arXiv:2006.08440); IAG Resolution No. 1 of 2015.
//   - IERS Conventions (2010), IERS Technical Note 36 (Petit and Luzum, eds.): the zero-tide tidal coefficients.
//   - GRS80 constants (gamma_e, k, e^2): Moritz, H. (2000), "Geodetic Reference System 1980",
//     Journal of Geodesy 74, 128-133.
//   - EGM2008 constants (GM, R): values of the header of the EGM2008 .gfc file (Pavlis et al. 2012, JGR 117, B04406).
//   - H0 = -0.31460 m (height of the permanent tide) is Eq. (16) of Makinen (2021), from Cartwright and Tayler (1971)
//     and Cartwright and Edden (1973). A0 = 4.4228e-8 m^-1 is Eq. (6.8c) of the IERS Conventions (2010), and the
//     permanent shift of C20 is A0 H0 k20 (Eq. 6.14); the same chapter gives -4.1736e-9 as the zero-tide minus
//     tide-free difference of C20 for EGM2008.
//   - k20 = 0.3 is the default degree-2 Love number (the nominal value of Makinen 2021); the same paper uses
//     k20 = 0.30190 for the conventional IERS value. The default is kept at 0.3 and can be overridden per call.
#pragma once
#include <string>
#include <vector>
#include "model.hpp"
#include "tide_system.hpp"

namespace gsh::tides
{
  // Numerical constants of the permanent-tide model. Changing them changes every conversion.
  inline constexpr double A0 = 4.4228e-8;          // amplitude factor of the permanent tidal potential [-]
  inline constexpr double H0 = -0.31460;           // height of the permanent tide [m]
  inline constexpr double GM = 0.3986004415e15;    // gravitational constant [m^3/s^2], EGM2008
  inline constexpr double R = 0.63781363e7;        // reference radius [m], EGM2008
  inline constexpr double GAMMA_E = 9.780326771;   // normal gravity at the equator [m/s^2], GRS80 (Somigliana)
  inline constexpr double GRS80_K = 0.001931851353;  // Somigliana constant k of GRS80 [-]
  inline constexpr double GRS80_E2 = 0.00669438002290; // first eccentricity squared of GRS80 [-]
  inline constexpr double K20_DEFAULT = 0.3;       // degree-2 Love number k20 [-] (see the file header)

  // Error message raised when C̄20 is requested for the mean-tide system. Its text is part of the
  // observable behaviour: keep it unchanged. The quoted phrase comes from arXiv:2006.08440 (see the file header).
  inline constexpr const char *NO_C20_MT =
      "There is no C̄20 in the mean-tide system (arXiv:2006.08440: 'there is no such thing as "
      "mean-tide C20'). C̄20 can only be converted between tide_free and zero_tide.";

  // ---- Scalar functions ----------------------------------------------------------------------
  // Change of C̄20 when going from `from` to `to`: offset(to) - offset(from), with offset(TF) = 0 and
  // offset(ZT) = k2 A0 H0 [-]. Throws std::invalid_argument if either system is MeanTide, because
  // mean-tide has no C̄20 of its own (NO_C20_MT).
  double c20_shift(TideSystem from, TideSystem to, double k2 = K20_DEFAULT);

  // Amount to ADD to the geoid height N [m] when going from `from` to `to`, at geodetic latitude
  // lat_deg. Does not apply to h or zeta (see convert_ellipsoidal_height and ihrf_separation_surface).
  double height_shift(double lat_deg, TideSystem from, TideSystem to, double k = K20_DEFAULT);
  // Geoid height N [m] at lat_deg converted from `from` to `to`: N + height_shift(...).
  double convert_geoid_height(double N, double lat_deg, TideSystem from, TideSystem to,
                              double k = K20_DEFAULT);

  // Vertical crust displacement h_T(phi) = (60.34 - 179.01 sin^2 phi - 1.82 sin^4 phi) mm, in metres.
  // Source: Makinen (2021).
  double crust_height_shift_m(double lat_deg);
  // Mean-minus-zero geoid separation H_T(phi) = (99.40 - 295.41 sin^2 phi - 0.42 sin^4 phi) mm, in
  // metres. Source: Makinen (2021).
  double geoid_tide_separation_m(double lat_deg);
  // Direct tidal gravity term g_T(phi) = (-30.49 + 90.95 sin^2 phi + 0.31 sin^4 phi) microGal, in mGal.
  // Source: Makinen (2021).
  double gravity_tide_mgal(double lat_deg);
  // Amount to ADD to the observed gravity [mGal] when going from `from` to `to`. Only the direct term
  // g_T is applied, and only for the mean-tide system.
  double gravity_tide_shift_mgal(double lat_deg, TideSystem from, TideSystem to);

  // Ellipsoidal height h [m] converted between systems: h(MT) = h(ZT) = h(NT) + h_T. ITRF/GNSS
  // coordinates are tide-free, so the usual input is `from` = TideFree.
  double convert_ellipsoidal_height(double h, double lat_deg, TideSystem from, TideSystem to);

  // Separation surface S(phi, lambda) = zeta_ZT - h_T(phi) + H_T(phi) [m], used to obtain the IHRF
  // normal height from the ellipsoidal height. Source: Makinen (2021)
  // (conversion to IHRF heights).
  double ihrf_separation_surface(double zeta_zero_tide, double lat_deg);

  // Copy of `model` with C̄20 moved to system `to`. The source system is model.tide_system. If that
  // string is empty or unknown, std::invalid_argument is thrown. If the source already equals `to`,
  // an unchanged copy is returned.
  Model convert_model_c20(const Model &model, TideSystem to, double k2 = K20_DEFAULT);
  // Copy of `model` with C̄20 converted from `from` to `to` (explicit source, not model.tide_system).
  // Returns an unchanged copy if from == to. Throws std::invalid_argument if one system is mean-tide,
  // or if the model degree is below 2. The tide_system field of the copy is set to `to`.
  Model convert_c20_from(const Model &model, TideSystem from, TideSystem to, double k2 = K20_DEFAULT);

  // Parses a system name (same as parse_tide in tide_system.hpp). Exposed here for convenience.
  inline TideSystem tide_from_name(const std::string &name) { return parse_tide(name); }

  // ---- Vectorized functions (OpenMP) ----------------------------------------------------------
  // Same computation as the scalar functions, element by element. The functions that take two or
  // more vectors throw std::invalid_argument if the sizes differ. Output has the size of the first
  // (or only) input vector.
  std::vector<double> crust_height_shift_m(const std::vector<double> &lat_deg);
  std::vector<double> geoid_tide_separation_m(const std::vector<double> &lat_deg);
  std::vector<double> gravity_tide_mgal(const std::vector<double> &lat_deg);
  std::vector<double> height_shift(const std::vector<double> &lat_deg, TideSystem from, TideSystem to,
                                   double k = K20_DEFAULT);
  std::vector<double> convert_geoid_height(const std::vector<double> &N, const std::vector<double> &lat_deg,
                                           TideSystem from, TideSystem to, double k = K20_DEFAULT);
  std::vector<double> gravity_tide_shift_mgal(const std::vector<double> &lat_deg, TideSystem from,
                                              TideSystem to);
  std::vector<double> convert_ellipsoidal_height(const std::vector<double> &h,
                                                 const std::vector<double> &lat_deg, TideSystem from,
                                                 TideSystem to);
  std::vector<double> ihrf_separation_surface(const std::vector<double> &zeta_zero_tide,
                                              const std::vector<double> &lat_deg);
}
