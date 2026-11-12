// functionals.hpp: geopotential functionals, as defined for GrafLab.
// Equations: Bucha & Janák (2019), "Definition of functionals of the geopotential used in GrafLab software"; "PDF n"
// in the comments below refers to the page numbering of that document.
// Local frame (LNOF): x north, y west, z up. Latitude: spherical (geocentric), in radians. Units: SI, with
// the angle-valued functionals (ξ, η, Θ) in radians.
// Bibliography: Bucha & Janák (2013), Computers & Geosciences 56,
// 186-196 (GrafLab software); Barthelmes (2013), GFZ STR09/02 (definition of the functionals); Petrovskaya &
// Vershkov (2006), J. Geodesy 80, 117-127 (non-singular gravity-gradient expressions; GrafLab uses a slightly
// modified form of them, Bucha & Janák 2013);
// Moritz (2000), J. Geodesy 74, 128-133 (GRS80 reference system).
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "ellipsoid.hpp"
#include "height.hpp"
#include "shs.hpp"

namespace gsh
{
  // Units of the functionals (SI):
  //   V, W, T, dV_phi, dV_lambda: m² s⁻² (per radian for the angular derivatives)
  //   dV_r, g_x, g_y, g_z, g, g_sa, dg_sa, Dg_sa, dg: m s⁻²
  //   V_*, W_rr, T_*: second derivatives or tensor components, s⁻² (1 Eötvös = 1e-9 s⁻²)
  //   xi, eta, Theta: radians
  //   zeta_ell, N, zeta: metres
  enum class Functional
  {
    // Gravitational potential V and its first derivatives (V: PDF 1; first derivatives: PDF 41)
    V, dV_r, dV_phi, dV_lambda,
    // Gravity potential W = V + centrifugal term, gravity vector in LNOF, and magnitude (PDF 36-48)
    W, g_x, g_y, g_z, g, g_sa, W_rr,
    // Tensor in spherical coordinates, as scaled second derivatives (PDF 2-8)
    V_rr, V_rphi, V_rlambda, V_phiphi, V_philambda, V_lambdalambda,
    // Tensor in the local LNOF frame, non-singular form (PDF 9-35)
    V_xx, V_xy, V_xz, V_yy, V_yz, V_zz,
    // Disturbing potential T, disturbance and anomaly of gravity, exact gravity disturbance (PDF 49-57, 54)
    T, dg_sa, Dg_sa, T_rr, dg,
    // Disturbing tensor, spherical coordinates and LNOF (same formulas as V with ΔC̄, ΔS̄; PDF 58-71)
    T_rphi, T_rlambda, T_phiphi, T_philambda, T_lambdalambda,
    T_xx, T_xy, T_xz, T_yy, T_yz, T_zz,
    // Deflections of the vertical, height anomaly, undulation and height anomaly on the surface (PDF 72-78)
    xi, eta, Theta, zeta_ell, N, zeta,
    Count
  };
  // Name as a string.
  const char *functional_name(Functional f);
  // All 45 functionals, in enum order.
  std::vector<Functional> all_functionals();

  // What to compute and how. Only the functionals listed in outputs are evaluated; the minimum number of synthesis
  // passes is derived from them (see compute_functionals).
  struct FunctionalRequest
  {
    std::vector<Functional> outputs;
    int nmin = 0;                                   // lowest degree included (n_min)
    int nmax = -1;                                  // highest degree; < 0 means the degree of the model
    Ellipsoid ellipsoid = Ellipsoid::GRS80();       // reference ellipsoid: r_ell, γ, ΔC̄ (Bucha and Janák 2019)
    std::shared_ptr<HeightSource> height;           // required for N and zeta: topography H(φ, λ), metres
    double density = 2670.0;                        // topographic density ρ, kg/m³ (Bucha and Janák 2019)
    double G = 6.67259e-11;                         // gravitational constant, m³ kg⁻¹ s⁻² (Bucha and Janák 2019)
    // Below this cos φ, dP and ddP are inaccurate and the point is evaluated with the oracle (reference_point).
    // 5e-3 is a choice of this project: the derivatives of P carry a factor 1/cos φ, so their accuracy degrades toward
    // the poles (see benchmark/meridian).
    double pole_cos_threshold = 5e-3;
    SynthOptions synth;                             // threads, ISA and kernel switches
  };

  struct FunctionalResult
  {
    // One array per requested functional, one entry per input point.
    std::map<Functional, AlignedVec<double>> values;
    // Points where the height source returned NaN; N and zeta are NaN there.
    std::int64_t points_without_height = 0;
    // Points handled by the oracle because cos φ < pole_cos_threshold.
    std::int64_t polar_points = 0;
    // Wall-clock seconds per stage, keyed by name: "model_synthesis", "disturbing_synthesis", "height" (height
    // evaluation) and "ellipsoid_synthesis" (synthesis at the ellipsoid foot points).
    std::map<std::string, double> timing_s;
  };

  // Computes the requested functionals at every point of pts (spherical latitude and longitude in radians, radius in
  // metres). The synthesis is planned to use the fewest passes: grad 0/1/2 on the model C̄, S̄, grad 0/1/2 on the
  // disturbing coefficients ΔC̄, ΔS̄, and T evaluated at the ellipsoid radius r_ell for N and ζ.
  // Points with cos φ < pole_cos_threshold use the oracle. Functionals that contain 1/cos φ (the singular set in
  // functionals.cpp) are undefined at the exact poles and are rejected there.
  // Throws std::invalid_argument if nmax > model.nmax, if N or zeta is requested without a height source, or if a
  // functional undefined at the pole (for example eta, or V_lambdalambda) is requested for a point at exactly ±π/2.
  // Cost: each synthesis pass is O(points × nmax²) operations, parallel over points; the oracle adds O(nmax²) per point.
  FunctionalResult compute_functionals(const Model &model, const Points &pts, const FunctionalRequest &req);
}
