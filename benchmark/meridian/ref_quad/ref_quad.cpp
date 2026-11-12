// ref_quad.cpp: extended-precision (IEEE binary128, 113-bit mantissa, ~1e-34) reference for the meridian benchmark.
//
// Purpose: say which engine (gsh or GrafLab) is closer to the exact value where the two differ by more than the
// 1e-12 criterion (dg and V_xz in benchmark/meridian). Nothing here is shared with gsh: the Legendre functions,
// their derivatives, the sums and the functional formulas are written again from the definitions below.
//
// What is evaluated (same definitions as gsh and GrafLab, on the same inputs):
//   * Fully normalized ALF (no Condon-Shortley phase, 4 pi normalisation), forward column recurrence in binary128.
//     The exponent range of binary128 (1e-4932) holds the sectoral underflow of double up to the pole region
//     where the terms that underflow are in the forbidden zone and are negligible (see README).
//   * dP/dphi from the identity  u dP_nm/dt = [(n+1) t P_nm - sqrt((2n+1)(n+m+1)(n-m+1)/(2n+3)) P_{n+1,m}] / u
//     (so d/dphi P = [..]/u), d2P/dphi2 from the Legendre ODE  P'' = tan(phi) P' - (n(n+1) - m^2/cos^2 phi) P.
//     At an exact pole the closed forms P_n0 = sqrt(2n+1) t^n and dP_n1/dphi = -t^(n+2) sqrt(2n+1) sqrt(n(n+1)/2) are
//     used (dg needs 1/cos(phi) and is not evaluated there; V_xx needs d2P and is not evaluated there).
//   * Raw sums S0 SR SP SL SRR SRP SPP at the longitude 0 (cos = 1, sin = 0), as in the project specification:
//       SR = sum q^n (n+1) A, SP = sum q^n A_phi, SL = sum q^n m (S cos - C sin) P, ...  with q = R/r.
//   * dg = |g| - gamma, g = (g_x, g_y, g_z) from dV/dr, dV/dphi, dV/dlambda plus the centrifugal terms, gamma = normal
//     gravity of the GRS80 spherical-harmonic series (even zonal terms up to degree 20, same formula as the project).
//   * V_xz = -(GM/r^3)(S_RP + S_P), V_zz = (GM/r^3) S_RR, V_xx = (GM/r^3)(S_PP - S_R): the last two are controls.
//
// Inputs that the engines also read in double are rounded to double first (coefficients as strtod, GM, R and the
// GRS80 constants), so that the reference measures the arithmetic of the engines and not their input representation.
// Positions: the geodetic latitude is converted to (r, phi) in binary128 from the exact decimal degree, or, with
// --engine-rad, from the double radian value lat_deg * (pi/180) that the engines use.
//
// Build (WSL, Ubuntu 24.04, gcc 13):  see build_and_run.sh.
// Usage: ref_quad <gfc> <out.csv> [--nmax N] [--lon-fixed] [--engine-rad] [--dlat-rad X] [--lat-list a,b,c]
#include <quadmath.h>
#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

typedef __float128 Q;

namespace
{
  std::vector<double> g_C, g_S;       // coefficient triangles, double (as the engines read them)
  std::vector<std::size_t> g_off;     // column offsets
  int g_nmax = 0, g_NM = 0;           // g_NM = nmax + 1 (the derivative of degree nmax needs P_{nmax+1})
  double g_GM = 0, g_R = 0;

  inline std::size_t idx(int n, int m) { return g_off[m] + static_cast<std::size_t>(n - m); }

  // Recurrence tables (independent of the point), binary128.
  std::vector<Q> g_a, g_b, g_df;      // a_nm, b_nm, dfac_nm
  std::vector<Q> g_sect;              // sectoral factors sqrt((2m+1)/(2m)), m >= 2
  std::vector<Q> g_c1;                // sqrt(2m+3)

  Q qpi() { return acosq(-1.0Q); }

  bool load_gfc(const char *path)
  {
    FILE *f = std::fopen(path, "r");
    if (!f)
      return false;
    static char line[1024];
    bool head = true;
    int maxdeg = 0;
    while (head && std::fgets(line, sizeof line, f))
    {
      char key[64], v[256];
      if (std::sscanf(line, "%63s %255s", key, v) == 2)
      {
        if (!std::strcmp(key, "earth_gravity_constant")) g_GM = std::strtod(v, nullptr);
        if (!std::strcmp(key, "radius")) g_R = std::strtod(v, nullptr);
        if (!std::strcmp(key, "max_degree")) maxdeg = std::atoi(v);
      }
      if (!std::strncmp(line, "end_of_head", 11))
        head = false;
    }
    if (g_nmax == 0)
      g_nmax = maxdeg;
    g_NM = g_nmax + 1;
    g_off.assign(g_NM + 2, 0);
    for (int m = 0; m <= g_NM; ++m)
      g_off[m + 1] = g_off[m] + static_cast<std::size_t>(g_NM - m + 1);
    const std::size_t tot = g_off[g_NM + 1];
    g_C.assign(tot, 0.0);
    g_S.assign(tot, 0.0);
    while (std::fgets(line, sizeof line, f))
    {
      if (std::strncmp(line, "gfc", 3))
        continue;
      int n, m;
      char cs[128], ss[128];
      if (std::sscanf(line + 3, "%d %d %127s %127s", &n, &m, cs, ss) != 4)
        continue;
      if (n > g_nmax)
        continue;
      for (char *p = cs; *p; ++p) if (*p == 'd' || *p == 'D') *p = 'e';
      for (char *p = ss; *p; ++p) if (*p == 'd' || *p == 'D') *p = 'e';
      g_C[idx(n, m)] = std::strtod(cs, nullptr);
      g_S[idx(n, m)] = std::strtod(ss, nullptr);
    }
    std::fclose(f);
    return true;
  }

  void build_tables()
  {
    const std::size_t tot = g_off[g_NM + 1];
    g_a.assign(tot, 0);
    g_b.assign(tot, 0);
    g_df.assign(tot, 0);
    g_sect.assign(g_NM + 2, 0);
    g_c1.assign(g_NM + 2, 0);
#pragma omp parallel for schedule(dynamic, 16)
    for (int m = 0; m <= g_NM; ++m)
    {
      g_c1[m] = sqrtq(static_cast<Q>(2 * m + 3));
      if (m >= 2)
        g_sect[m] = sqrtq(static_cast<Q>(2 * m + 1) / static_cast<Q>(2 * m));
      for (int n = m; n <= g_NM; ++n)
      {
        const long long N = n, M = m;
        const Q den = static_cast<Q>(N * N - M * M);
        if (n >= m + 2)
        {
          g_a[idx(n, m)] = sqrtq(static_cast<Q>(4 * N * N - 1) / den);
          g_b[idx(n, m)] = sqrtq(static_cast<Q>((2 * N + 1) * (N + M - 1) * (N - M - 1)) / (static_cast<Q>(2 * N - 3) * den));
        }
        if (n < g_NM)
          g_df[idx(n, m)] = sqrtq(static_cast<Q>((2 * N + 1) * (N + M + 1) * (N - M + 1)) / static_cast<Q>(2 * N + 3));
      }
    }
  }

  struct Out
  {
    Q r, phi, dg, g, gamma, vxz, vzz, vxx;
  };

  // Column of P_nm, m fixed, n = m..NM (u > 0).
  inline void column(int m, Q t, Q pmm, Q *P)
  {
    P[0] = pmm;
    if (m + 1 <= g_NM)
      P[1] = g_c1[m] * t * pmm;
    for (int n = m + 2; n <= g_NM; ++n)
    {
      const std::size_t i = idx(n, m);
      P[n - m] = g_a[i] * t * P[n - m - 1] - g_b[i] * P[n - m - 2];
    }
  }

  // Normal gravity of the GRS80 series (even zonals up to degree 20), at (r, phi) with t = sin, u = cos (u > 0).
  Q normal_gravity(Q r, Q t, Q u, Q &gx_out, Q &gz_out)
  {
    const Q a = 6378137.0Q, e2 = static_cast<Q>(0.00669438002290), GMe = static_cast<Q>(3.986005e14),
            J2 = static_cast<Q>(1.08263e-3), om = static_cast<Q>(7.292115e-5);
    // P_n0 and dP_n0/dphi for n <= 20 by the same machinery as the model (column m = 0).
    Q P[23];
    P[0] = 1;
    P[1] = sqrtq(3.0Q) * t;
    for (int n = 2; n <= 21; ++n)
    {
      const Q N = n;
      P[n] = sqrtq((4 * N * N - 1) / (N * N)) * t * P[n - 1] -
             ((N - 1) / N) * sqrtq((2 * N + 1) / (2 * N - 3)) * P[n - 2];
    }
    const Q q = a / r;
    const Q c20 = -J2 / sqrtq(5.0Q);
    Q sr = 0, sp = 0, qn = 1;
    for (int k = 0; k <= 10; ++k)
    {
      const int n = 2 * k;
      const Q N = n;
      Q c;
      if (k == 0) c = 1;
      else if (k == 1) c = c20;
      else
      {
        const Q K = k;
        const Q sign = (k % 2 == 0) ? 1 : -1;
        c = sign * 3 * powq(e2, K) / ((2 * K + 1) * (2 * K + 3) * sqrtq(4 * K + 1)) *
            (1 - K - 5 * sqrtq(5.0Q) * K * c20 / e2);
      }
      const Q dP = ((N + 1) * t * P[n] - sqrtq((2 * N + 1) * (N + 1) * (N + 1) / (2 * N + 3)) * P[n + 1]) / u;
      sr += qn * (N + 1) * c * P[n];
      sp += qn * c * dP;
      qn *= q * q;
    }
    const Q dVr = -GMe / (r * r) * sr, dVp = GMe / r * sp;
    gx_out = (dVp - om * om * r * r * u * t) / r;
    gz_out = dVr + om * om * r * u * u;
    return hypotq(gx_out, gz_out);
  }

  Out evaluate(double lat_deg, bool engine_rad, double dlat)
  {
    const Q a = 6378137.0Q, e2 = static_cast<Q>(0.00669438002290);
    const Q om = static_cast<Q>(7.292115e-5);
    const bool pole = std::fabs(lat_deg) == 90.0;
    Q t, u, r;
    if (pole)
    {
      t = lat_deg > 0 ? 1 : -1;
      u = 0;
      r = a * sqrtq(1 - e2);
    }
    else
    {
      Q phig;
      if (engine_rad)
        phig = static_cast<Q>(lat_deg * (3.14159265358979323846 / 180.0));
      else
        phig = static_cast<Q>(lat_deg) * qpi() / 180;
      phig += static_cast<Q>(dlat);
      const Q s = sinq(phig), c = cosq(phig);
      const Q N = a / sqrtq(1 - e2 * s * s);
      const Q x = N * c, z = N * (1 - e2) * s;
      const Q phi = atan2q(z, x);
      t = sinq(phi);
      u = cosq(phi);
      r = hypotq(x, z);
    }
    const Q GM = static_cast<Q>(g_GM), R = static_cast<Q>(g_R), q = R / r;
    std::vector<Q> qn(g_NM + 2);
    qn[0] = 1;
    for (int n = 1; n <= g_NM; ++n)
      qn[n] = qn[n - 1] * q;

    Q S0 = 0, SR = 0, SP = 0, SL = 0, SRR = 0, SRP = 0, SPP = 0;
    std::vector<Q> P(g_NM + 2);
    if (pole)
    {
      for (int n = 0; n <= g_nmax; ++n)
      {
        const Q N = n;
        // m = 0: P = sqrt(2n+1) t^n, dP/dphi = 0.
        const Q Pn0 = sqrtq(2 * N + 1) * powq(t, N);
        const Q C0 = static_cast<Q>(g_C[idx(n, 0)]);
        S0 += qn[n] * C0 * Pn0;
        SR += qn[n] * (N + 1) * C0 * Pn0;
        SRR += qn[n] * (N + 1) * (N + 2) * C0 * Pn0;
        // m = 1: dP/dphi = -t^(n+2) sqrt(2n+1) sqrt(n(n+1)/2); at lambda = 0 base = C_n1.
        if (n >= 1)
        {
          const Q dP = -powq(t, N + 2) * sqrtq(2 * N + 1) * sqrtq(N * (N + 1) / 2);
          const Q C1 = static_cast<Q>(g_C[idx(n, 1)]);
          SP += qn[n] * C1 * dP;
          SRP += qn[n] * (N + 1) * C1 * dP;
        }
      }
    }
    else
    {
      Q pmm = 1;
      for (int m = 0; m <= g_nmax; ++m)
      {
        if (m == 1) pmm = sqrtq(3.0Q) * u;
        else if (m >= 2) pmm = g_sect[m] * u * pmm;
        column(m, t, pmm, P.data());
        const Q M = m;
        for (int n = m; n <= g_nmax; ++n)
        {
          const std::size_t i = idx(n, m);
          const Q Pn = P[n - m], Pn1 = P[n - m + 1];
          const Q N = n;
          const Q dP = ((N + 1) * t * Pn - g_df[i] * Pn1) / u;
          const Q ddP = (t / u) * dP - (N * (N + 1) - M * M / (u * u)) * Pn;
          const Q C = static_cast<Q>(g_C[i]), S = static_cast<Q>(g_S[i]);
          const Q w = qn[n];
          const Q A = C * Pn, Ap = C * dP, App = C * ddP, Al = M * S * Pn;
          S0 += w * A;
          SR += w * (N + 1) * A;
          SP += w * Ap;
          SL += w * Al;
          SRR += w * (N + 1) * (N + 2) * A;
          SRP += w * (N + 1) * Ap;
          SPP += w * App;
        }
      }
    }
    Out o{};
    o.r = r;
    o.phi = atan2q(t, u);
    const Q k3 = GM / (r * r * r);
    o.vxz = -k3 * (SRP + SP);
    o.vzz = k3 * SRR;
    o.vxx = pole ? Q(NAN) : k3 * (SPP - SR);
    if (pole)
    {
      o.dg = Q(NAN);
      o.g = Q(NAN);
      o.gamma = Q(NAN);
    }
    else
    {
      const Q dVr = -GM / (r * r) * SR, dVp = GM / r * SP, dVl = GM / r * SL;
      const Q gx = (dVp - om * om * r * r * u * t) / r;
      const Q gy = -dVl / (r * u);
      const Q gz = dVr + om * om * r * u * u;
      o.g = sqrtq(gx * gx + gy * gy + gz * gz);
      Q ngx, ngz;
      o.gamma = normal_gravity(r, t, u, ngx, ngz);
      o.dg = o.g - o.gamma;
    }
    (void)S0;
    return o;
  }

  std::string fq(Q v)
  {
    char b[128];
    quadmath_snprintf(b, sizeof b, "%.33Qe", v);
    return b;
  }
}

int main(int argc, char **argv)
{
  if (argc < 3)
  {
    std::fprintf(stderr, "usage: ref_quad <gfc> <out.csv> [--nmax N] [--engine-rad] [--dlat-rad X] [--lat-list a,b,c]\n");
    return 2;
  }
  bool engine_rad = false;
  double dlat = 0;
  std::vector<double> lats;
  for (int i = 3; i < argc; ++i)
  {
    if (!std::strcmp(argv[i], "--nmax") && i + 1 < argc) g_nmax = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--engine-rad")) engine_rad = true;
    else if (!std::strcmp(argv[i], "--dlat-rad") && i + 1 < argc) dlat = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--lat-list") && i + 1 < argc)
    {
      char *s = argv[++i];
      for (char *tok = std::strtok(s, ","); tok; tok = std::strtok(nullptr, ","))
        lats.push_back(std::atof(tok));
    }
  }
  if (lats.empty())
    for (int k = -90; k <= 90; ++k)
      lats.push_back(k);
  if (!load_gfc(argv[1]))
  {
    std::fprintf(stderr, "cannot read %s\n", argv[1]);
    return 1;
  }
  std::fprintf(stderr, "ref_quad: nmax=%d GM=%.17g R=%.17g, %zu points, engine_rad=%d, dlat=%g\n", g_nmax, g_GM, g_R,
               lats.size(), engine_rad ? 1 : 0, dlat);
  build_tables();
  std::vector<Out> res(lats.size());
  double t0 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < static_cast<int>(lats.size()); ++i)
  {
    res[i] = evaluate(lats[i], engine_rad, dlat);
#pragma omp critical
    std::fprintf(stderr, "  lat %7.2f done (%.0f s)\n", lats[i], omp_get_wtime() - t0);
  }
  FILE *f = std::fopen(argv[2], "w");
  std::fprintf(f, "lat_deg,r_m,phi_sph_rad,dg_ms2,g_ms2,gamma_ms2,V_xz_s2,V_xx_s2,V_zz_s2\n");
  for (std::size_t i = 0; i < lats.size(); ++i)
  {
    const Out &o = res[i];
    std::fprintf(f, "%.17g,%s,%s,%s,%s,%s,%s,%s,%s\n", lats[i], fq(o.r).c_str(), fq(o.phi).c_str(), fq(o.dg).c_str(),
                 fq(o.g).c_str(), fq(o.gamma).c_str(), fq(o.vxz).c_str(), fq(o.vxx).c_str(), fq(o.vzz).c_str());
  }
  std::fclose(f);
  return 0;
}
