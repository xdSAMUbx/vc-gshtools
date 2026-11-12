// reference.cpp: the sequential oracle. Evaluates the literal formulas of Bucha and Janák (2019) (raw sums and
// local tensor, PDF 9-35) on the full triangle of
// Legendre functions, with no singularity at the poles. Slow (O(nmax^2) per point); used for validation and for
// near-polar points in compute_functionals() and grid synthesis.
#include <cmath>
#include <vector>
#include "kernels.hpp"
#include "shs.hpp"

namespace gsh
{
  namespace
  {
    // Square root with the convention of Bucha and Janák (2019): a negative radicand gives 0 (it only happens where the
    // associated P̄ is zero anyway).
    inline double sq(double v) { return v > 0.0 ? std::sqrt(v) : 0.0; }
    // Kronecker delta δ_ab, used for the 1 + δ factors of the local-tensor coefficients.
    inline double kron(int a, int b) { return a == b ? 1.0 : 0.0; }
  }

  // Evaluates the oracle at one point. Steps: (1) Legendre P̄, dP̄/dφ, d²P̄/dφ² on the full triangle with MSFC and the
  // non-singular identities of Bucha and Janák (2019); (2) raw sums of Bucha and Janák (2019) over m and n; (3) local tensor of Bucha and Janák (2019), times GM/r³.
  ReferencePoint reference_point(const Model &model, const legendre::Coefs &coefs, double lat, double lon, double r, int nmin, int nmax)
  {
    if (nmax < 0)
      nmax = model.nmax;
    const int cn = coefs.nmax;
    std::vector<double> P(tri_size(cn)), dP(tri_size(cn)), ddP(tri_size(cn));
    const bool pole = std::fabs(lat) == PI / 2;
    const double t = pole ? (lat > 0 ? 1.0 : -1.0) : std::sin(lat);
    const double u = pole ? 0.0 : std::cos(lat);
    legendre::msfc(t, u, coefs, P);
    legendre::dphi(cn, P, dP);
    legendre::dphi(cn, dP, ddP);

    const auto Pb = [&](int n, int k) -> double
    {
      if (n < 0 || k < 0 || k > n || n > cn)
        return 0.0;
      return P[idx_tri(n, k, cn)];
    };

    ReferencePoint out{};
    using namespace kernels;
    const double q = model.R / r;
    std::vector<double> qn(nmax + 1);
    qn[0] = 1.0;
    for (int n = 1; n <= nmax; ++n)
      qn[n] = qn[n - 1] * q;

    double Vxx = 0, Vyy = 0, Vzz = 0, Vxz = 0, Vxy = 0, Vyz = 0;

    // Raw sums run over (nmax+1)(nmax+2)/2 terms (2.4e6 at 2190). A plain double accumulator loses about sqrt(N) ε
    // of relative accuracy (1e-5 m²/s² in V at degree 2190, measured against the kernel), which is too much for a
    // reference. Neumaier's compensated summation keeps the error near ε. Requires /fp:precise (no fast-math).
    struct Neumaier
    {
      double s = 0.0, c = 0.0;
      void add(double x)
      {
        const double t = s + x;
        c += (std::fabs(s) >= std::fabs(x)) ? (s - t) + x : (x - t) + s;
        s = t;
      }
      double value() const { return s + c; }
    };
    Neumaier ks[10];

    for (int m = 0; m <= nmax; ++m)
    {
      const double cm = std::cos(m * lon), sm = std::sin(m * lon);
      const double M = m;
      for (int n = (m > nmin ? m : nmin); n <= nmax; ++n)
      {
        const double N = n, C = model.c(n, m), S = model.s(n, m);
        const auto i = idx_tri(n, m, cn);
        const double base = C * cm + S * sm, cross = S * cm - C * sm;
        const double A = base * P[i], Ap = base * dP[i], App = base * ddP[i];
        const double Al = M * cross * P[i], Apl = M * cross * dP[i], All = M * M * base * P[i];
        const double w1 = N + 1, w2 = N + 2, Q = qn[n];
        ks[F_S0].add(Q * A);
        ks[F_SR].add(Q * w1 * A);
        ks[F_SP].add(Q * Ap);
        ks[F_SL].add(Q * Al);
        ks[F_SRR].add(Q * w1 * w2 * A);
        ks[F_SRP].add(Q * w1 * Ap);
        ks[F_SRL].add(Q * w1 * Al);
        ks[F_SPP].add(Q * App);
        ks[F_SPL].add(Q * Apl);
        ks[F_SLL].add(Q * All);
      }
    }
    for (int k = 0; k < 10; ++k)
      out.raw[k] = ks[k].value();

    // Local tensor (Bucha and Janák 2019): sum over the signed order ms in [-n, n]. For ms >= 0 the coefficient is C̄_n,ms and
    // the longitude factor is cos(mλ) (Q_m); for ms < 0 it is S̄_n,|ms| with sin(|m|λ). The coefficients a, b, c, d,
    // g, h, mu, nu, beta, gamma below are the closed forms of Bucha and Janák (2019) (PDF 9-35), with P̄ evaluated at n, n-1 and
    // the orders |m| - 2 ... |m| + 2 (Pb returns 0 outside the triangle).
    for (int n = nmin; n <= nmax; ++n)
    {
      const double N = n, Q = qn[n];
      const double kappa = n > 0 ? std::sqrt((2 * N + 1) / (2 * N - 1)) : 0.0;
      for (int ms = -n; ms <= n; ++ms)
      {
        const int am = std::abs(ms);
        const double A = am;
        const double cnm = ms >= 0 ? model.c(n, am) : model.s(n, am);
        const double qm = ms >= 0 ? std::cos(am * lon) : std::sin(am * lon);    // Q_m
        const double qmn = ms > 0 ? std::sin(am * lon) : std::cos(am * lon);    // Q_{-m}
        const double sgn = ms > 0 ? 1.0 : (ms < 0 ? -1.0 : 0.0);
        const double term = cnm * qm;

        double a = 0, b, c;
        if (am >= 2)
          a = sq(1 + kron(am, 2)) / 4 * sq(N * N - (A - 1) * (A - 1)) * sq(N + A) * sq(N - A + 2);
        if (am <= 1)
        {
          b = (N + A + 1) * (N + A + 2) / (2 * (A + 1));
          c = sq(1 + kron(am, 0)) / 4 * sq(N * N - (A + 1) * (A + 1)) * sq(N - A) * sq(N + A + 2);
        }
        else
        {
          b = (N * N + A * A + 3 * N + 2) / 2;
          c = 0.25 * sq(N * N - (A + 1) * (A + 1)) * sq(N - A) * sq(N + A + 2);
        }
        const double p0 = Pb(n, am), pm2 = Pb(n, am - 2), pp2 = Pb(n, am + 2);
        const double bracket = a * pm2 + b * p0 + c * pp2;
        Vxx += Q * term * (bracket - (N + 1) * (N + 2) * p0);
        Vyy -= Q * term * bracket;
        Vzz += Q * (N + 1) * (N + 2) * term * p0;

        double beta = 0, gamma;
        if (am >= 1)
        {
          beta = (N + 2) / 2 * sq(1 + kron(am, 1)) * sq(N + A) * sq(N - A + 1);
          gamma = -(N + 2) / 2 * sq(N - A) * sq(N + A + 1);
        }
        else
          gamma = -(N + 2) * sq(N * (N + 1) / 2);
        Vxz += Q * term * (beta * Pb(n, am - 1) + gamma * Pb(n, am + 1));

        if (ms != 0)
        {
          const double tn = cnm * qmn;
          double d = 0, g;
          if (am >= 2)
          {
            d = -sgn / 4 * kappa * sq(1 + kron(am, 2)) * sq(N * N - (A - 1) * (A - 1)) * sq(N + A) * sq(N + A - 2);
            g = static_cast<double>(ms) / 2 * kappa * sq(N + A) * sq(N - A);
          }
          else
            g = sgn / 4 * kappa * sq(N + 1) * sq(N - 1) * (N + 2);
          const double h = sgn / 4 * kappa * sq(N * N - (A + 1) * (A + 1)) * sq(N - A) * sq(N - A - 2);
          Vxy += Q * tn * (d * Pb(n - 1, am - 2) + g * Pb(n - 1, am) + h * Pb(n - 1, am + 2));
          const double mu = -sgn * (N + 2) / 2 * kappa * sq(1 + kron(am, 1)) * sq(N + A) * sq(N + A - 1);
          const double nu = -sgn * (N + 2) / 2 * kappa * sq(N - A) * sq(N - A - 1);
          Vyz += Q * tn * (mu * Pb(n - 1, am - 1) + nu * Pb(n - 1, am + 1));
        }
      }
    }
    const double s3 = model.GM / (r * r * r);
    out.vxx = s3 * Vxx;
    out.vyy = s3 * Vyy;
    out.vzz = s3 * Vzz;
    out.vxz = s3 * Vxz;
    out.vxy = s3 * Vxy;
    out.vyz = s3 * Vyz;
    return out;
  }
}
