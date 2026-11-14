#include "functionals.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include "kernels.hpp"

namespace gsh
{
  using namespace kernels;

  // Names of the functionals as strings. The static_assert below ties the array to the enum order.
  const char *functional_name(Functional f)
  {
    static const char *const names[] = {
        "V", "dV_r", "dV_phi", "dV_lambda",
        "W", "g_x", "g_y", "g_z", "g", "g_sa", "W_rr",
        "V_rr", "V_rphi", "V_rlambda", "V_phiphi", "V_philambda", "V_lambdalambda",
        "V_xx", "V_xy", "V_xz", "V_yy", "V_yz", "V_zz",
        "T", "dg_sa", "Dg_sa", "T_rr", "dg",
        "T_rphi", "T_rlambda", "T_phiphi", "T_philambda", "T_lambdalambda",
        "T_xx", "T_xy", "T_xz", "T_yy", "T_yz", "T_zz",
        "xi", "eta", "Theta", "zeta_ell", "N", "zeta"};
    static_assert(sizeof(names) / sizeof(names[0]) == static_cast<int>(Functional::Count));
    return names[static_cast<int>(f)];
  }

  // Lists the enum values in declaration order (the same order as functional_name).
  std::vector<Functional> all_functionals()
  {
    std::vector<Functional> v;
    for (int i = 0; i < static_cast<int>(Functional::Count); ++i)
      v.push_back(static_cast<Functional>(i));
    return v;
  }

  namespace
  {
    using Clock = std::chrono::steady_clock;
    double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

    // Gradient order that each functional needs from the model (model: C̄, S̄) and from the disturbing model (disturb:
    // ΔC̄, ΔS̄); -1 means not needed. normal: needs the normal gravity γ. ellipsoid_points: needs the N / ζ path.
    struct Need
    {
      int model = -1, disturb = -1;
      bool normal = false, ellipsoid_points = false;
    };
    // Maps one functional to the passes it requires. The grouping follows the sums required by each functional (Bucha and Janák 2019).
    Need need_of(Functional f)
    {
      using F = Functional;
      Need n;
      switch (f)
      {
      case F::V:
        n.model = 0;
        break;
      case F::W:
        n.model = 0;
        break;
      case F::dV_r:
      case F::dV_phi:
      case F::dV_lambda:
      case F::g_x:
      case F::g_y:
      case F::g_z:
      case F::g:
      case F::g_sa:
        n.model = 1;
        break;
      case F::V_rr:
      case F::W_rr:
        n.model = 2;
        break;
      case F::V_rphi:
      case F::V_rlambda:
      case F::V_phiphi:
      case F::V_philambda:
      case F::V_lambdalambda:
      case F::V_xx:
      case F::V_xy:
      case F::V_xz:
      case F::V_yy:
      case F::V_yz:
      case F::V_zz:
        n.model = 2;
        break;
      case F::T:
        n.disturb = 0;
        break;
      case F::dg_sa:
      case F::Dg_sa:
        n.disturb = 1;
        break;
      case F::T_rr:
      case F::T_rphi:
      case F::T_rlambda:
      case F::T_phiphi:
      case F::T_philambda:
      case F::T_lambdalambda:
      case F::T_xx:
      case F::T_xy:
      case F::T_xz:
      case F::T_yy:
      case F::T_yz:
      case F::T_zz:
        n.disturb = 2;
        break;
      case F::dg:
        n.disturb = 1;
        n.normal = true;
        break;
      case F::xi:
      case F::eta:
      case F::Theta:
        n.disturb = 1;
        n.normal = true;
        break;
      case F::zeta_ell:
        n.disturb = 0;
        n.normal = true;
        break;
      case F::N:
      case F::zeta:
        n.ellipsoid_points = true;
        break;
      default:
        break;
      }
      return n;
    }

    // Local (LNOF) tensor of V or T at one point, from its raw sums S and with k3 = GM/r³ and the latitude φ.
    // The result is the covariant Hessian of the potential in the LNOF frame. In terms of partial derivatives of V
    // (V_r = ∂V/∂r, V_φλ = ∂²V/∂φ∂λ, ...), the formulas coded below are:
    //   V_xx = V_φφ/r² + V_r/r               V_yy = V_λλ/(r² cos²φ) + V_r/r - tanφ V_φ/r²
    //   V_zz = V_rr                           V_xz = V_rφ/r - V_φ/r²
    //   V_xy = -(V_φλ + tanφ V_λ)/(r² cosφ)  V_yz = -(V_rλ/r - V_λ/r²)/cosφ
    // Note the 1/cosφ factors in V_xy and V_yz. These are the same tensor as the non-singular form of PDF 9-35 that
    // reference_point() evaluates; test_shs.cpp checks both (LaplaceAndHessianConsistency, TensorAgainstFiniteDifferences).
    // Units: s⁻² (the caller passes k3 = GM/r³).
    struct Tensor
    {
      double xx, xy, xz, yy, yz, zz;
    };
    // Computes the six LNOF components of the tensor for one point; the formulas are documented on struct Tensor.
    inline Tensor local_tensor(const double *S, double k3, double lat)
    {
      const double c = std::cos(lat), tn = std::sin(lat) / c;
      Tensor t;
      t.xx = k3 * (S[F_SPP] - S[F_SR]);
      t.xy = -k3 * (S[F_SPL] + tn * S[F_SL]) / c;
      t.xz = -k3 * (S[F_SRP] + S[F_SP]);
      t.yy = k3 * (-S[F_SLL] / (c * c) - S[F_SR] - tn * S[F_SP]);
      t.yz = k3 * (S[F_SRL] + S[F_SL]) / c;
      t.zz = k3 * S[F_SRR];
      return t;
    }

    // Relative cost of one synthesis pass per point, by gradient order (grad 0 = 1). Only used to decide whether fusing the
    // model and disturbing passes pays off.
    double pass_cost(int grad) { return grad == 0 ? 1.0 : (grad == 1 ? 3.5 : 6.0); }
    // Cost of the zonal normal-field pass (at most 231 terms of degree <= 20, against about nmax^2/2 terms for a full pass),
    // about 1e-4 of a full pass; rounded up to 0.01 to stay conservative.
    constexpr double kZonalPassCost = 0.01;
  }

  // Evaluates the requested functionals. The structure is: (1) decide the passes from the requested outputs,
  // (2) run the model and disturbing syntheses, (3) replace near-polar points with the oracle, (4) fill each group
  // of functionals from the raw sums, (5) for N and zeta, evaluate the height and the foot-point synthesis.
  FunctionalResult compute_functionals(const Model &model, const Points &pts, const FunctionalRequest &req)
  {
    const int nmax = req.nmax < 0 ? model.nmax : req.nmax;
    if (nmax > model.nmax)
      throw std::invalid_argument("compute_functionals: nmax is larger than the model degree");
    using F = Functional;
    const std::int64_t np = pts.size();
    FunctionalResult res;
    const std::set<F> want(req.outputs.begin(), req.outputs.end());
    const auto wants = [&](F f)
    { return want.count(f) > 0; };

    int gm = -1, gd = -1;
    bool need_normal = false, need_ell = false;
    for (F f : want)
    {
      const Need n = need_of(f);
      gm = std::max(gm, n.model);
      gd = std::max(gd, n.disturb);
      need_normal = need_normal || n.normal;
      need_ell = need_ell || n.ellipsoid_points;
    }
    // Allocates (zero-filled) the output array of f and returns a reference to it.
    const auto out = [&](F f) -> AlignedVec<double> &
    {
      auto &v = res.values[f];
      v.assign(np, 0.0);
      return v;
    };
    const Ellipsoid &ell = req.ellipsoid;
    const double w2 = ell.omega * ell.omega; // ω², Earth's angular velocity squared (rad² s⁻²)

    // Exact poles (Bucha and Janák 2019): functionals with 1/cos φ are undefined at φ = ±π/2 and are rejected there. The set holds
    // the functionals that contain ∂V/∂λ / cos φ directly (g_y, V_rλ, V_φλ, V_λλ, T_rλ, T_φλ, T_λλ, η) and those
    // built from them (g, dg, Θ).
    static const std::set<F> singular = {F::g_y, F::g, F::dg, F::V_rlambda, F::V_philambda, F::V_lambdalambda,
                                         F::T_rlambda, F::T_philambda, F::T_lambdalambda, F::eta, F::Theta};
    std::vector<std::int64_t> polar;
    for (std::int64_t i = 0; i < np; ++i)
    {
      if (std::cos(pts.lat[i]) < req.pole_cos_threshold)
        polar.push_back(i);
      if (std::abs(pts.lat[i]) == PI / 2)
        for (F f : want)
          if (singular.count(f))
            throw std::invalid_argument(std::string("The functional ") + functional_name(f) + " is not defined at an exact pole");
    }
    res.polar_points = static_cast<std::int64_t>(polar.size());

    // One recurrence table serves both the model and the disturbing synthesis (same degree).
    auto rec = make_recurrence(nmax);
    // Disturbing coefficients ΔC̄ = C̄ - C̄^Ell (GM^Ell/GM)(a/R)^n and ΔS̄ = S̄ (PDF 50-51). Built only when needed.
    const Model dist = (gd >= 0 || need_ell) ? disturbing_model(model, ell) : Model{};
    RawFields rm, rd;

    // Fused path: model = dist + zonal, coefficient by coefficient (normal_zonal_model), and synthesis is linear. So one pass on
    // the disturbing coefficients at grad max(gm, gd) gives rd, and the model sums are rm = rd + Z, where Z is the synthesis of
    // the zonal model (degree <= 20, grad gm). Only the sums are added, never subtracted, so T (from rd) keeps its accuracy.
    // Used when both sums are needed and the cost model says it is cheaper; otherwise the two full passes below run.
    const bool fuse = req.fuse_model_and_disturbing && gm >= 0 && gd >= 0 && nmax >= 2 &&
                      pass_cost(std::max(gm, gd)) + kZonalPassCost < pass_cost(gm) + pass_cost(gd);

    auto t0 = Clock::now();
    if (fuse)
    {
      rd = synth_points(prepare(dist, req.nmin, nmax, rec), pts, std::max(gm, gd), req.synth);
      res.timing_s["disturbing_synthesis"] = secs(t0);

      t0 = Clock::now();
      // The zonal model has degree min(20, nmax). If nmin is above it, all its terms are below nmin and Z = 0.
      const int nz = std::min(20, nmax);
      const bool zonal_terms = req.nmin <= nz;
      RawFields z;
      if (zonal_terms)
      {
        const Model zonal = normal_zonal_model(model, ell);
        const auto rec_z = nz == nmax ? rec : make_recurrence(nz);
        z = synth_points(prepare(zonal, req.nmin, nz, rec_z), pts, gm, req.synth);
      }
      rm.grad = gm;
      rm.f.assign(static_cast<std::size_t>(fields_for_grad(gm)), AlignedVec<double>(static_cast<std::size_t>(np)));
      for (std::size_t k = 0; k < rm.f.size(); ++k)
      {
        const double *zk = zonal_terms ? z.f[k].data() : nullptr;
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < np; ++i)
          rm.f[k][i] = zk ? rd.f[k][i] + zk[i] : rd.f[k][i];
      }
      // At polar points the model value S0 (grad 0) is taken from the model directly, as in the two-pass path. Other model
      // fields at polar points are replaced by reference_point() below.
      if (gm == 0 && !polar.empty())
      {
        Points sub;
        for (std::int64_t id : polar)
        {
          sub.lat.push_back(pts.lat[id]);
          sub.lon.push_back(pts.lon[id]);
          sub.r.push_back(pts.r[id]);
        }
        const RawFields ms = synth_points(prepare(model, req.nmin, nmax, rec), sub, 0, req.synth);
        for (std::size_t j = 0; j < polar.size(); ++j)
          rm.f[F_S0][polar[j]] = ms.f[F_S0][j];
      }
      res.timing_s["model_synthesis"] = secs(t0);
    }
    else
    {
      if (gm >= 0)
      {
        rm = synth_points(prepare(model, req.nmin, nmax, rec), pts, gm, req.synth);
        res.timing_s["model_synthesis"] = secs(t0);
      }
      t0 = Clock::now();
      if (gd >= 0)
      {
        rd = synth_points(prepare(dist, req.nmin, nmax, rec), pts, gd, req.synth);
        res.timing_s["disturbing_synthesis"] = secs(t0);
      }
    }

    // Near the poles, the column-wise dP and ddP (1/cos φ form) lose accuracy. For those points the raw sums are
    // replaced by reference_point(), which evaluates the full triangle with the non-singular identities of Bucha and Janák (2019)
    if (!polar.empty() && (gm >= 1 || gd >= 1))
    {
      const auto coefs = legendre::make_coefs(nmax);
      for (std::int64_t id : polar)
      {
        if (gm >= 1)
        {
          const auto ref = reference_point(model, coefs, pts.lat[id], pts.lon[id], pts.r[id], req.nmin, nmax);
          for (std::size_t k = 0; k < rm.f.size(); ++k)
            rm.f[k][id] = ref.raw[k];
        }
        if (gd >= 1)
        {
          const auto ref = reference_point(dist, coefs, pts.lat[id], pts.lon[id], pts.r[id], req.nmin, nmax);
          for (std::size_t k = 0; k < rd.f.size(); ++k)
            rd.f[k][id] = ref.raw[k];
        }
      }
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();

    const auto raw_at = [&](const RawFields &rf, std::int64_t i, double *S)
    {
      for (std::size_t k = 0; k < rf.f.size(); ++k)
        S[k] = rf.f[k][i];
    };

    // Group V: potential, gravity vector and second derivatives, from the model's raw sums rm.
    if (gm >= 0)
    {
      const double GM = model.GM;
      const auto fill = [&](F f, auto &&fn)
      {
        if (!wants(f))
          return;
        auto &v = out(f);
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < np; ++i)
          v[i] = fn(i);
      };
      const auto S = [&](int k, std::int64_t i)
      { return rm.f[k][i]; };
      const auto R = [&](std::int64_t i)
      { return pts.r[i]; };
      const auto cosl = [&](std::int64_t i)
      { return std::cos(pts.lat[i]); };
      const auto sinl = [&](std::int64_t i)
      { return std::sin(pts.lat[i]); };

      // V = GM/r S_0 (PDF 1). W = V + V_c with V_c = ½ ω² r² cos²φ (PDF 36-37).
      fill(F::V, [&](auto i)
           { return GM / R(i) * S(F_S0, i); });
      fill(F::W, [&](auto i)
           { return GM / R(i) * S(F_S0, i) + 0.5 * w2 * R(i) * R(i) * cosl(i) * cosl(i); });
      if (gm >= 1)
      {
        // First derivatives of V (PDF 41): ∂V/∂r = -GM/r² S_R, ∂V/∂φ = GM/r S_P, ∂V/∂λ = GM/r S_L.
        // Gravity vector in LNOF (PDF 37-45), with ∂V_c/∂r = ω² r cos²φ, ∂V_c/∂φ = -ω² r² cosφ sinφ, ∂V_c/∂λ = 0:
        //   g_x = (∂V/∂φ + ∂V_c/∂φ)/r,  g_y = -(∂V/∂λ)/(r cosφ),  g_z = ∂V/∂r + ∂V_c/∂r.
        //   g = |(g_x, g_y, g_z)|. g_sa = |g_z| is the spherical approximation (within PDF 37-45).
        // g_y and g contain 1/cosφ and are singular at the poles. The exact gravity disturbance dg is not computed here:
        // see the group T below (it is evaluated from the disturbing-potential derivatives).
        const auto dVr = [&](auto i)
        { return -GM / (R(i) * R(i)) * S(F_SR, i); };
        const auto dVp = [&](auto i)
        { return GM / R(i) * S(F_SP, i); };
        const auto dVl = [&](auto i)
        { return GM / R(i) * S(F_SL, i); };
        const auto gx = [&](auto i)
        { return (dVp(i) - w2 * R(i) * R(i) * cosl(i) * sinl(i)) / R(i); };
        const auto gy = [&](auto i)
        { return -dVl(i) / (R(i) * cosl(i)); };
        const auto gz = [&](auto i)
        { return dVr(i) + w2 * R(i) * cosl(i) * cosl(i); };
        const auto gmod = [&](auto i)
        { const double a = gx(i), b = gy(i), c = gz(i); return std::sqrt(a * a + b * b + c * c); };
        fill(F::dV_r, dVr);
        fill(F::dV_phi, dVp);
        fill(F::dV_lambda, dVl);
        fill(F::g_x, gx);
        fill(F::g_y, gy);
        fill(F::g_z, gz);
        fill(F::g, gmod);
        fill(F::g_sa, [&](auto i)
             { return std::abs(gz(i)); });
      }
      if (gm == 2)
      {
        // Spherical-coordinate tensor, scaled second derivatives (PDF 3-8), with k3 = GM/r³:
        //   V_rr = k3 S_RR;  W_rr = V_rr + ω² cos²φ;  V_rφ = -k3 S_RP;  V_rλ = -k3 S_RL / cosφ;
        //   V_φφ = k3 S_PP;  V_φλ = k3 S_PL / cosφ;   V_λλ = -k3 S_LL / cos²φ.
        const auto k3 = [&](auto i)
        { return GM / (R(i) * R(i) * R(i)); };
        fill(F::V_rr, [&](auto i)
             { return k3(i) * S(F_SRR, i); });
        fill(F::W_rr, [&](auto i)
             { return k3(i) * S(F_SRR, i) + w2 * cosl(i) * cosl(i); });
        fill(F::V_rphi, [&](auto i)
             { return -k3(i) * S(F_SRP, i); });
        fill(F::V_rlambda, [&](auto i)
             { return -k3(i) / cosl(i) * S(F_SRL, i); });
        fill(F::V_phiphi, [&](auto i)
             { return k3(i) * S(F_SPP, i); });
        fill(F::V_philambda, [&](auto i)
             { return k3(i) / cosl(i) * S(F_SPL, i); });
        fill(F::V_lambdalambda, [&](auto i)
             { return -k3(i) / (cosl(i) * cosl(i)) * S(F_SLL, i); });
      }
    }

    // Local tensor of V (or T) at every point, with the polar points taken from reference_point() (see local_tensor).
    const auto local = [&](const RawFields &rf, const Model &m, F xx, F xy, F xz, F yy, F yz, F zz)
    {
      if (!(wants(xx) || wants(xy) || wants(xz) || wants(yy) || wants(yz) || wants(zz)))
        return;
      auto &a = out(xx), &b = out(xy), &c = out(xz), &d = out(yy), &e = out(yz), &f = out(zz);
#pragma omp parallel for schedule(static)
      for (std::int64_t i = 0; i < np; ++i)
      {
        double S[10];
        raw_at(rf, i, S);
        const double r = pts.r[i];
        const Tensor t = local_tensor(S, m.GM / (r * r * r), pts.lat[i]);
        a[i] = t.xx;
        b[i] = t.xy;
        c[i] = t.xz;
        d[i] = t.yy;
        e[i] = t.yz;
        f[i] = t.zz;
      }
      if (!polar.empty())
      {
        const auto coefs = legendre::make_coefs(nmax);
        for (std::int64_t id : polar)
        {
          const auto ref = reference_point(m, coefs, pts.lat[id], pts.lon[id], pts.r[id], req.nmin, nmax);
          a[id] = ref.vxx;
          b[id] = ref.vxy;
          c[id] = ref.vxz;
          d[id] = ref.vyy;
          e[id] = ref.vyz;
          f[id] = ref.vzz;
        }
      }
    };
    if (gm == 2)
      local(rm, model, F::V_xx, F::V_xy, F::V_xz, F::V_yy, F::V_yz, F::V_zz);
    if (gd == 2)
      local(rd, dist, F::T_xx, F::T_xy, F::T_xz, F::T_yy, F::T_yz, F::T_zz);

    // Group T: disturbing potential and its derivatives, from the disturbing raw sums rd (ΔC̄, ΔS̄).
    if (gd >= 0)
    {
      const double GM = model.GM;
      const auto fill = [&](F f, auto &&fn)
      {
        if (!wants(f))
          return;
        auto &v = out(f);
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < np; ++i)
          v[i] = fn(i);
      };
      const auto S = [&](int k, std::int64_t i)
      { return rd.f[k][i]; };
      const auto R = [&](std::int64_t i)
      { return pts.r[i]; };
      const auto cosl = [&](std::int64_t i)
      { return std::cos(pts.lat[i]); };
      const auto gamma = [&](std::int64_t i)
      { return ell.normal_gravity(pts.lat[i], R(i)).gamma; };
      // T = GM/r S_0 (PDF 49). zeta_ell = T/γ is the height anomaly on the ellipsoid (PDF 77).
      fill(F::T, [&](auto i)
           { return GM / R(i) * S(F_S0, i); });
      fill(F::zeta_ell, [&](auto i)
           { return GM / R(i) * S(F_S0, i) / gamma(i); });
      if (gd >= 1)
      {
        // δg_sa = -∂T/∂r = GM/r² S_R (PDF 55). Δg_sa = δg_sa - 2T/r = GM/r² (S_R - 2 S_0) (PDF 56).
        fill(F::dg_sa, [&](auto i)
             { return GM / (R(i) * R(i)) * S(F_SR, i); });
        fill(F::Dg_sa, [&](auto i)
             { return GM / (R(i) * R(i)) * (S(F_SR, i) - 2.0 * S(F_S0, i)); });
        // Deflections of the vertical (PDF 72-74): ξ = -∂T/∂φ / (rγ) = -GM/(r²γ) S_P, η = -GM/(r²γ cosφ) S_L,
        // Θ = sqrt(ξ² + η²). ξ, η, Θ are in radians; η and Θ carry 1/cosφ.
        const auto xi_ = [&](auto i)
        { return -GM / (R(i) * R(i) * gamma(i)) * S(F_SP, i); };
        const auto eta_ = [&](auto i)
        { return -GM / (R(i) * R(i) * gamma(i) * cosl(i)) * S(F_SL, i); };
        // Exact gravity disturbance dg = |g| - γ (PDF 54), g = γ + δ, with γ the normal gravity vector (γ_x, 0, γ_z) of the
        // ellipsoid and δ = -∇T the disturbing gravity: δ_x = (∂T/∂φ)/r = GM/r² S_P, δ_y = -(∂T/∂λ)/(r cosφ) =
        // -GM/(r² cosφ) S_L, δ_z = ∂T/∂r = -GM/r² S_R (disturbing sums). Then
        //   dg = |γ + δ| - |γ| = (2 γ·δ + |δ|²) / (|γ + δ| + |γ|),
        // which is algebraically |g| - γ but has no difference of two numbers of 9.8 m/s² and no monopole term in the
        // sums, so it does not lose digits to cancellation (benchmark/meridian, comparison against binary128).
        // Singular at the poles (δ_y has 1/cosφ), like g.
        fill(F::dg, [&](auto i)
             {
               const auto nrm = ell.normal_gravity(pts.lat[i], R(i));
               const double k = GM / (R(i) * R(i));
               const double dx = k * S(F_SP, i);
               const double dy = -k * S(F_SL, i) / cosl(i);
               const double dz = -k * S(F_SR, i);
               const double gx = nrm.g_x + dx, gz = nrm.g_z + dz;
               const double num = 2.0 * (nrm.g_x * dx + nrm.g_z * dz) + dx * dx + dy * dy + dz * dz;
               return num / (std::sqrt(gx * gx + dy * dy + gz * gz) + nrm.gamma); });
        fill(F::xi, xi_);
        fill(F::eta, eta_);
        fill(F::Theta, [&](auto i)
             { return std::hypot(xi_(i), eta_(i)); });
      }
      if (gd == 2)
      {
        // Same formulas as the model tensor above, with the disturbing sums (PDF 57 for T_rr; PDF 58-71 for the rest).
        const auto k3 = [&](auto i)
        { return GM / (R(i) * R(i) * R(i)); };
        fill(F::T_rr, [&](auto i)
             { return k3(i) * S(F_SRR, i); });
        fill(F::T_rphi, [&](auto i)
             { return -k3(i) * S(F_SRP, i); });
        fill(F::T_rlambda, [&](auto i)
             { return -k3(i) / cosl(i) * S(F_SRL, i); });
        fill(F::T_phiphi, [&](auto i)
             { return k3(i) * S(F_SPP, i); });
        fill(F::T_philambda, [&](auto i)
             { return k3(i) / cosl(i) * S(F_SPL, i); });
        fill(F::T_lambdalambda, [&](auto i)
             { return -k3(i) / (cosl(i) * cosl(i)) * S(F_SLL, i); });
      }
    }

    // N and zeta (PDF 75-78). T and δg_sa are evaluated on the reference ellipsoid, at the same spherical latitude φ:
    // the foot point has r_ell(φ) = a sqrt((1-e²)/(1-e² cos²φ)) (ell.radius) and the height H is read from the height
    // source at the geodetic latitude of that point. Then N = (T - 2πGρH²)/γ and ζ = T/γ - δg_sa (H + N)/γ, with
    // γ = γ(r_ell, φ), taken from Ellipsoid::normal_gravity (the harmonic γ_SH of Bucha and Janák 2019). Points without height data (NaN) get NaN and are
    // counted in points_without_height.
    if (need_ell)
    {
      if (!req.height)
        throw std::invalid_argument("N and zeta require a height source (HeightSource)");
      t0 = Clock::now();
      Points foot = pts;
      AlignedVec<double> lat_geod(np);
      AlignedVec<double> gamma_ell(np);
#pragma omp parallel for schedule(static)
      for (std::int64_t i = 0; i < np; ++i)
      {
        foot.r[i] = ell.radius(pts.lat[i]);
        lat_geod[i] = rad2deg(ell.geodetic_of_spherical(pts.lat[i]));
        gamma_ell[i] = ell.normal_gravity(pts.lat[i], foot.r[i]).gamma;
      }
      std::vector<double> H;
      req.height->evaluate(foot, lat_geod, H, req.synth);
      res.timing_s["height"] = secs(t0);

      t0 = Clock::now();
      const auto rf = synth_points(prepare(dist, req.nmin, nmax, rec), foot, 1, req.synth);
      res.timing_s["ellipsoid_synthesis"] = secs(t0);

      AlignedVec<double> *pN = wants(F::N) ? &out(F::N) : nullptr;
      AlignedVec<double> *pZ = wants(F::zeta) ? &out(F::zeta) : nullptr;
      std::int64_t missing = 0;
      const double cterm = 2.0 * PI * req.G * req.density;
#pragma omp parallel for schedule(static) reduction(+ : missing)
      for (std::int64_t i = 0; i < np; ++i)
      {
        const double r = foot.r[i], g = gamma_ell[i], h = H[i];
        if (std::isnan(h))
        {
          ++missing;
          if (pN)
            (*pN)[i] = nan;
          if (pZ)
            (*pZ)[i] = nan;
          continue;
        }
        const double T = model.GM / r * rf.f[F_S0][i];
        const double N = (T - cterm * h * h) / g;
        if (pN)
          (*pN)[i] = N;
        if (pZ)
        {
          const double dgsa = model.GM / (r * r) * rf.f[F_SR][i];
          (*pZ)[i] = T / g - dgsa * (h + N) / g;
        }
      }
      res.points_without_height = missing;
    }
    return res;
  }
}
