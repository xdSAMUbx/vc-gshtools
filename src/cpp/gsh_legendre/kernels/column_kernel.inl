// column_kernel.inl: MSFC recurrence fused with the synthesis, without storing the triangle P.
// Template over the lane type V (simd.hpp), B vectors per tile and the gradient order G (0: values;
// 1: values and first derivatives; 2: also second derivatives). Included once per ISA translation unit,
// after simd.hpp.
//
// Math for one point, with q = R/r, t = sin φ, u = cos φ:
//   y_n = q^n P̄_nm(φ), with the recurrence y_n = a (t q) y_{n-1} - b q^2 y_{n-2} (coefficients from legendre::Coefs).
//   Sectorial: y_mm = d_m (u q) y_{m-1,m-1}.
//   Mantissa and exponent per lane: y is kept near unity and the true value is y * R^-k with R = 2^500
//   (the MSFC idea of Fukushima 2012).
//   For each column m the kernel accumulates lumps a = Σ_n w_n y_n C_nm and b = Σ_n w_n y_n S_nm, where the
//   weight w_n depends on the field, and combines them with cos(mλ) and sin(mλ) to form the raw sums of shs.hpp.
//   Derivatives use the identities of Fukushima (2012, part II) and Bosch (2000) (with the 1/u form; the
//   1/u singularity is why the callers replace near-polar points by the oracle).
#include <cmath>
#include <cstring>
#include <vector>
#include "core.hpp"
#include "kernels.hpp"

namespace gsh::GSH_ISA_NS
{
  using kernels::PrepView;
  using kernels::TileOpts;

  // Recurrence steps between rescale and skip checks (one block of the column loop).
  inline constexpr int RESCALE_EVERY = 32;
  // Orders between exact recomputations of cos(mλ), sin(mλ) when rotate_trig is on; bounds the rotation drift.
  inline constexpr int TRIG_RESET = 64;
  // Mantissa threshold for exponent k = 1. The value is y * R^-1 = y * 2^-500, so |y| < 2^300 means
  // |value| < 2^-200, which is negligible in double precision next to the other terms.
  inline constexpr double NEGLIGIBLE_K1 = 0x1p300;

  // Number of (a, b) lump pairs per column for gradient order G: 1 (values), 3 (first derivatives), 6 (second).
  // Each pair is one vector per lane group, so a column produces 2P vectors per tile vector.
  template <int G>
  inline constexpr int LUMP_PAIRS = G == 0 ? 1 : (G == 1 ? 3 : 6);

  // Per-tile state: the per-lane geometry, constant along the m loop.
  //   t, u = sin φ, cos φ;  q = R/r;  tq = t q;  q2 = q²;  uq = u q: used by the recurrences.
  //   qu = q/u, tu = t/u, ru2 = 1/u²: only for gradient orders >= 1 (they contain 1/cos φ).
  //   c1, s1 = cos λ, sin λ: starting values of the trigonometric rotation; lon keeps λ for exact trig.
  //   polar_u_max: largest u in the tile, used by the polar truncation test.
  template <int B>
  struct TileState
  {
    V t[B], u[B], q[B], tq[B], q2[B], uq[B];
    V qu[B], tu[B], ru2[B];
    V c1[B], s1[B];
    alignas(64) double lon[B * WIDTH];
    double polar_u_max;
  };

  // Fills TileState for B * WIDTH points. Exact poles are set explicitly (t = ±1, u = 0), so that cos(π/2) does not
  // leave a rounding residue in u, which would feed the recurrences a spurious non-zero value.
  template <int B>
  GSH_INLINE void init_state(TileState<B> &st, const double *lat, const double *lon, const double *r, double R, int G)
  {
    alignas(64) double buf_t[WIDTH], buf_u[WIDTH], buf_q[WIDTH], buf_c[WIDTH], buf_s[WIDTH];
    double umax = 0.0;
    for (int v = 0; v < B; ++v)
    {
      for (int w = 0; w < WIDTH; ++w)
      {
        const int i = v * WIDTH + w;
        const double phi = lat[i];
        // Exact poles, so that cos(π/2) leaves no rounding residue.
        if (std::fabs(phi) == PI / 2)
        {
          buf_t[w] = phi > 0 ? 1.0 : -1.0;
          buf_u[w] = 0.0;
        }
        else
        {
          buf_t[w] = std::sin(phi);
          buf_u[w] = std::cos(phi);
        }
        umax = buf_u[w] > umax ? buf_u[w] : umax;
        buf_q[w] = R / r[i];
        buf_c[w] = std::cos(lon[i]);
        buf_s[w] = std::sin(lon[i]);
        st.lon[i] = lon[i];
      }
      st.t[v] = load(buf_t);
      st.u[v] = load(buf_u);
      st.q[v] = load(buf_q);
      st.c1[v] = load(buf_c);
      st.s1[v] = load(buf_s);
      st.tq[v] = mul(st.t[v], st.q[v]);
      st.q2[v] = mul(st.q[v], st.q[v]);
      st.uq[v] = mul(st.u[v], st.q[v]);
      if (G > 0)
      {
        const V one = set1(1.0);
        st.qu[v] = div(st.q[v], st.u[v]);
        st.tu[v] = div(st.t[v], st.u[v]);
        st.ru2[v] = div(one, mul(st.u[v], st.u[v]));
      }
    }
    st.polar_u_max = umax;
  }

  // Advances the sectorial term per lane: x <- d_m (u q) x. Lanes whose |x| falls below R^-1 are multiplied by R
  // and their exponent k is incremented (exact zeros are left alone).
  template <int B>
  GSH_INLINE void next_sectorial(const TileState<B> &st, double dm, V *x, V *k)
  {
    const V d = set1(dm), rinv = set1(R_INV), big = set1(R_BIG), zero = set1(0.0), one = set1(1.0);
    for (int v = 0; v < B; ++v)
    {
      x[v] = mul(mul(d, st.uq[v]), x[v]);
      M mk = mand(lt(abs(x[v]), rinv), ne(x[v], zero));
      while (any(mk))
      {
        x[v] = blend(mk, x[v], mul(x[v], big));
        k[v] = blend(mk, k[v], add(k[v], one));
        mk = mand(lt(abs(x[v]), rinv), ne(x[v], zero));
      }
    }
  }

  // Full column m. Inputs x0, k0: sectorial value and exponent per lane. Output: lump[2 * P * B] holds the
  // column lumps already on absolute scale (R^-k applied). Returns false when every lane has k >= 2, i.e. the
  // whole column is negligible and the caller may skip the combine step.
  template <int B, int G>
  GSH_INLINE bool column(const TileState<B> &st, const PrepView &pv, const TileOpts &opts, int m,
                         const V *x0, const V *k0, V *lump)
  {
    constexpr int P = LUMP_PAIRS<G>;
    const int nmax = pv.nmax;
    const std::int64_t off = col_offset(m, nmax);
    const double *ab = pv.ab + 2 * off;
    const double *cs = pv.cs + 2 * off;
    const double *ff = (G > 0) ? pv.f + off : nullptr;
    const int len = nmax - m;

    const V zero = set1(0.0), one = set1(1.0), two = set1(2.0);
    const V big = set1(R_BIG), rinv = set1(R_INV), negl = set1(NEGLIGIBLE_K1);
    const V m2 = set1(static_cast<double>(m) * m);

    V y1[B], y2[B], k[B];
    for (int i = 0; i < 2 * P * B; ++i)
      lump[i] = zero;

    // Term n = m (the sectorial value itself). For G >= 1, dy = -n (t/u) y is q^m dP̄_mm/dφ (Fukushima 2012, II).
    {
      const double N = m;
      const V C = set1(cs[0]), S = set1(cs[1]);
      const V w1 = set1(N + 1.0), w2 = set1(N + 2.0), nn = set1(N), nn1 = set1(N * (N + 1.0));
      for (int v = 0; v < B; ++v)
      {
        y2[v] = x0[v];
        k[v] = k0[v];
        const V y = x0[v];
        V *L = lump + 2 * P * v;
        const V yC = mul(C, y), yS = mul(S, y);
        L[0] = add(L[0], yC);
        L[1] = add(L[1], yS);
        if constexpr (G >= 1)
        {
          const V dy = sub(zero, mul(mul(nn, st.tu[v]), y));
          L[2] = fmadd(w1, yC, L[2]);
          L[3] = fmadd(w1, yS, L[3]);
          L[4] = fmadd(C, dy, L[4]);
          L[5] = fmadd(S, dy, L[5]);
          if constexpr (G == 2)
          {
            const V ddy = fmadd(st.tu[v], dy, mul(fmsub(m2, st.ru2[v], nn1), y));
            const V wpr = mul(w1, w2);
            L[6] = fmadd(wpr, yC, L[6]);
            L[7] = fmadd(wpr, yS, L[7]);
            L[8] = fmadd(mul(w1, C), dy, L[8]);
            L[9] = fmadd(mul(w1, S), dy, L[9]);
            L[10] = fmadd(C, ddy, L[10]);
            L[11] = fmadd(S, ddy, L[11]);
          }
        }
      }
    }

    if (len > 0)
    {
      // Term n = m+1 (b = 0 for this step, so only a is used). Its predecessor is the sectorial value y2.
      {
        const double N = m + 1.0;
        const V a = set1(ab[2]), C = set1(cs[2]), S = set1(cs[3]);
        for (int v = 0; v < B; ++v)
        {
          const V y = mul(mul(a, st.tq[v]), y2[v]);
          V *L = lump + 2 * P * v;
          const V yC = mul(C, y), yS = mul(S, y);
          L[0] = add(L[0], yC);
          L[1] = add(L[1], yS);
          if constexpr (G >= 1)
          {
            const V w1 = set1(N + 1.0), nn = set1(N);
            const V dy = fnmadd(mul(nn, st.tu[v]), y, mul(mul(set1(ff[1]), st.qu[v]), y2[v]));
            L[2] = fmadd(w1, yC, L[2]);
            L[3] = fmadd(w1, yS, L[3]);
            L[4] = fmadd(C, dy, L[4]);
            L[5] = fmadd(S, dy, L[5]);
            if constexpr (G == 2)
            {
              const V w2 = set1(N + 2.0), nn1 = set1(N * (N + 1.0));
              const V ddy = fmadd(st.tu[v], dy, mul(fmsub(m2, st.ru2[v], nn1), y));
              const V wpr = mul(w1, w2);
              L[6] = fmadd(wpr, yC, L[6]);
              L[7] = fmadd(wpr, yS, L[7]);
              L[8] = fmadd(mul(w1, C), dy, L[8]);
              L[9] = fmadd(mul(w1, S), dy, L[9]);
              L[10] = fmadd(C, ddy, L[10]);
              L[11] = fmadd(S, ddy, L[11]);
            }
          }
          y1[v] = y;
        }
      }

      int j = 2;
      while (j <= len)
      {
        const int jend = (j + RESCALE_EVERY - 1 < len) ? j + RESCALE_EVERY - 1 : len;

        // Can any lane contribute in this block? A lane contributes if k == 0, or if k == 1 and |y| is large.
        // When no lane can contribute, the block only advances the recurrence (no lump updates).
        bool accumulate = !opts.skip_forbidden;
        if (!accumulate)
          for (int v = 0; v < B; ++v)
          {
            const M m0 = lt(k[v], one);                                         // k == 0
            const M m1 = mand(lt(k[v], two), gt(abs(y1[v]), negl));             // k == 1 and |y| large
            const M m2c = mand(lt(k[v], two), gt(abs(y2[v]), negl));
            if (any(mor(m0, mor(m1, m2c))))
            {
              accumulate = true;
              break;
            }
          }

        if (!accumulate)
        {
          const double *pab = ab + 2 * j;
          for (int i = j; i <= jend; ++i, pab += 2)
          {
            const V a = set1(pab[0]), b = set1(pab[1]);
            for (int v = 0; v < B; ++v)
            {
              const V y = fnmadd(mul(b, st.q2[v]), y2[v], mul(mul(a, st.tq[v]), y1[v]));
              y2[v] = y1[v];
              y1[v] = y;
            }
          }
        }
        else if constexpr (G == 0)
        {
          const double *pab = ab + 2 * j, *pcs = cs + 2 * j;
          for (int i = j; i <= jend; ++i, pab += 2, pcs += 2)
          {
            const V a = set1(pab[0]), b = set1(pab[1]), C = set1(pcs[0]), S = set1(pcs[1]);
            for (int v = 0; v < B; ++v)
            {
              const V y = fnmadd(mul(b, st.q2[v]), y2[v], mul(mul(a, st.tq[v]), y1[v]));
              lump[2 * v] = fmadd(C, y, lump[2 * v]);
              lump[2 * v + 1] = fmadd(S, y, lump[2 * v + 1]);
              y2[v] = y1[v];
              y1[v] = y;
            }
          }
        }
        else
        {
          const double *pab = ab + 2 * j, *pcs = cs + 2 * j, *pf = ff + j;
          double N = static_cast<double>(m + j);
          for (int i = j; i <= jend; ++i, pab += 2, pcs += 2, ++pf, N += 1.0)
          {
            const V a = set1(pab[0]), b = set1(pab[1]), C = set1(pcs[0]), S = set1(pcs[1]), fv = set1(*pf);
            const V w1 = set1(N + 1.0), nn = set1(N);
            for (int v = 0; v < B; ++v)
            {
              const V y = fnmadd(mul(b, st.q2[v]), y2[v], mul(mul(a, st.tq[v]), y1[v]));
              const V dy = fnmadd(mul(nn, st.tu[v]), y, mul(mul(fv, st.qu[v]), y1[v]));
              V *L = lump + 2 * P * v;
              const V yC = mul(C, y), yS = mul(S, y);
              L[0] = add(L[0], yC);
              L[1] = add(L[1], yS);
              L[2] = fmadd(w1, yC, L[2]);
              L[3] = fmadd(w1, yS, L[3]);
              L[4] = fmadd(C, dy, L[4]);
              L[5] = fmadd(S, dy, L[5]);
              if constexpr (G == 2)
              {
                const V w2 = set1(N + 2.0), nn1 = set1(N * (N + 1.0));
                const V ddy = fmadd(st.tu[v], dy, mul(fmsub(m2, st.ru2[v], nn1), y));
                const V wpr = mul(w1, w2);
                L[6] = fmadd(wpr, yC, L[6]);
                L[7] = fmadd(wpr, yS, L[7]);
                L[8] = fmadd(mul(w1, C), dy, L[8]);
                L[9] = fmadd(mul(w1, S), dy, L[9]);
                L[10] = fmadd(C, ddy, L[10]);
                L[11] = fmadd(S, ddy, L[11]);
              }
              y2[v] = y1[v];
              y1[v] = y;
            }
          }
        }

        // Per-lane rescale: when |y1| > R_BIG and k > 0, divide the mantissas and the lump accumulators by R and
        // decrement k. The product y * R^-k is unchanged, so the result is exact up to rounding.
        for (int v = 0; v < B; ++v)
        {
          const M mr = mand(gt(abs(y1[v]), big), gt(k[v], zero));
          if (any(mr))
          {
            y1[v] = blend(mr, y1[v], mul(y1[v], rinv));
            y2[v] = blend(mr, y2[v], mul(y2[v], rinv));
            k[v] = blend(mr, k[v], sub(k[v], one));
            V *L = lump + 2 * P * v;
            for (int l = 0; l < 2 * P; ++l)
              L[l] = blend(mr, L[l], mul(L[l], rinv));
          }
        }
        j = jend + 1;
      }
    }

    // Apply R^-k: k == 0 -> 1; k == 1 -> R^-1; k >= 2 -> 0. Mantissas stay below about R = 2^500, so a term with
    // k >= 2 is below 2^-500 and is dropped.
    bool significant = false;
    for (int v = 0; v < B; ++v)
    {
      const M m0 = lt(k[v], one);
      const M m1 = lt(k[v], two);
      const V scale = blend(m0, blend(m1, zero, rinv), one);
      V *L = lump + 2 * P * v;
      if (any(m1))
        significant = true;
      for (int l = 0; l < 2 * P; ++l)
        L[l] = mul(L[l], scale);
    }
    return significant;
  }

  // Combines the lumps of one column with cos(mλ) and sin(mλ) and adds the result to the raw sums.
  // For example A_n = a cos(mλ) + b sin(mλ) and A^λ_n = m (b cos(mλ) - a sin(mλ)); G = 2 adds the m² terms.
  template <int B, int G>
  GSH_INLINE void combine(int m, const V *lump, const V *cm, const V *sm, V *f)
  {
    constexpr int P = LUMP_PAIRS<G>;
    constexpr int NF = kernels::fields_for_grad(G);
    const V mv = set1(static_cast<double>(m)), m2 = set1(static_cast<double>(m) * m);
    for (int v = 0; v < B; ++v)
    {
      const V *L = lump + 2 * P * v;
      V *F = f + NF * v;
      const V c = cm[v], s = sm[v];
      const V A0 = fmadd(L[0], c, mul(L[1], s));
      F[kernels::F_S0] = add(F[kernels::F_S0], A0);
      if constexpr (G >= 1)
      {
        F[kernels::F_SR] = add(F[kernels::F_SR], fmadd(L[2], c, mul(L[3], s)));
        F[kernels::F_SP] = add(F[kernels::F_SP], fmadd(L[4], c, mul(L[5], s)));
        F[kernels::F_SL] = fmadd(mv, fmsub(L[1], c, mul(L[0], s)), F[kernels::F_SL]);
        if constexpr (G == 2)
        {
          F[kernels::F_SRR] = add(F[kernels::F_SRR], fmadd(L[6], c, mul(L[7], s)));
          F[kernels::F_SRP] = add(F[kernels::F_SRP], fmadd(L[8], c, mul(L[9], s)));
          F[kernels::F_SRL] = fmadd(mv, fmsub(L[3], c, mul(L[2], s)), F[kernels::F_SRL]);
          F[kernels::F_SPP] = add(F[kernels::F_SPP], fmadd(L[10], c, mul(L[11], s)));
          F[kernels::F_SPL] = fmadd(mv, fmsub(L[5], c, mul(L[4], s)), F[kernels::F_SPL]);
          F[kernels::F_SLL] = fmadd(m2, A0, F[kernels::F_SLL]);
        }
      }
    }
  }

  // Exact cos(mλ) and sin(mλ) for every lane, computed from the stored longitude (no accumulated rotation error).
  template <int B>
  GSH_INLINE void exact_trig(const TileState<B> &st, int m, V *cm, V *sm)
  {
    alignas(64) double c[WIDTH], s[WIDTH];
    for (int v = 0; v < B; ++v)
    {
      for (int w = 0; w < WIDTH; ++w)
      {
        const double a = m * st.lon[v * WIDTH + w];
        c[w] = std::cos(a);
        s[w] = std::sin(a);
      }
      cm[v] = load(c);
      sm[v] = load(s);
    }
  }

  // True when column m cannot contribute to any lane: m - nmax * umax > threshold. The left side grows with m, so
  // once it is true it stays true and the caller may stop the m loop.
  inline bool polar_skip(double threshold, int m, int nmax, double umax)
  {
    return threshold >= 0.0 && static_cast<double>(m) - nmax * umax > threshold;
  }

  // Serial tile: one thread processes all the lanes of the tile. Sectorial values and trigonometry are carried
  // along the m loop, so the triangle is never stored. The caller parallelises over tiles.
  template <int B, int G>
  void tile_serial(const PrepView &pv, const TileOpts &opts, const double *lat, const double *lon, const double *r, double *out)
  {
    constexpr int P = LUMP_PAIRS<G>;
    constexpr int NF = kernels::fields_for_grad(G);
    TileState<B> st;
    init_state<B>(st, lat, lon, r, pv.R, G);

    alignas(64) V x[B], k[B], cm[B], sm[B], lump[2 * P * B], f[NF * B];
    const V zero = set1(0.0), one = set1(1.0);
    for (int v = 0; v < B; ++v)
    {
      x[v] = one;
      k[v] = zero;
      cm[v] = one;
      sm[v] = zero;
    }
    for (int i = 0; i < NF * B; ++i)
      f[i] = zero;

    for (int m = 0; m <= pv.nmax; ++m)
    {
      if (polar_skip(opts.polar_threshold, m, pv.nmax, st.polar_u_max))
        break; // the condition is monotone in m, so no later order contributes
      if (m > 0)
      {
        next_sectorial<B>(st, pv.d[m], x, k);
        // Rotation only for G = 0: in G = 1 and G = 2 the weights m and m² amplify the rotation error. Measurements
        // in this project found no gain there (the raw timings are not archived).
        if (G == 0 && opts.rotate_trig && m % TRIG_RESET != 0)
          for (int v = 0; v < B; ++v)
          {
            const V c2 = fnmadd(sm[v], st.s1[v], mul(cm[v], st.c1[v]));
            sm[v] = fmadd(sm[v], st.c1[v], mul(cm[v], st.s1[v]));
            cm[v] = c2;
          }
        else
          exact_trig<B>(st, m, cm, sm);
      }
      if (column<B, G>(st, pv, opts, m, x, k, lump))
        combine<B, G>(m, lump, cm, sm, f);
    }

    const int lanes = B * WIDTH;
    for (int v = 0; v < B; ++v)
      for (int nf = 0; nf < NF; ++nf)
        store(out + nf * lanes + v * WIDTH, f[NF * v + nf]);
  }

  // Tile with internal parallelism over m, for few points: the tile loop alone cannot keep all threads busy.
  // The sectorial values are precomputed serially. Each thread owns its accumulators, and the results are reduced
  // at the end. The loop order mirrors the CHarm code path (Bucha 2026).
  template <int B, int G>
  void tile_mparallel(const PrepView &pv, const TileOpts &opts, const double *lat, const double *lon, const double *r, double *out)
  {
    constexpr int P = LUMP_PAIRS<G>;
    constexpr int NF = kernels::fields_for_grad(G);
    TileState<B> st;
    init_state<B>(st, lat, lon, r, pv.R, G);

    const int nmax = pv.nmax;
    AlignedVec<double> sx((nmax + 1) * B * WIDTH), sk((nmax + 1) * B * WIDTH);
    {
      V x[B], k[B];
      for (int v = 0; v < B; ++v)
      {
        x[v] = set1(1.0);
        k[v] = set1(0.0);
      }
      for (int m = 0; m <= nmax; ++m)
      {
        if (m > 0)
          next_sectorial<B>(st, pv.d[m], x, k);
        for (int v = 0; v < B; ++v)
        {
          store(sx.data() + (m * B + v) * WIDTH, x[v]);
          store(sk.data() + (m * B + v) * WIDTH, k[v]);
        }
      }
    }

    const int lanes = B * WIDTH;
    for (int i = 0; i < NF * lanes; ++i)
      out[i] = 0.0;

#pragma omp parallel
    {
      alignas(64) V x[B], k[B], cm[B], sm[B], lump[2 * P * B], f[NF * B];
      for (int i = 0; i < NF * B; ++i)
        f[i] = set1(0.0);
#pragma omp for schedule(dynamic, 4)
      for (int m = 0; m <= nmax; ++m)
      {
        if (polar_skip(opts.polar_threshold, m, nmax, st.polar_u_max))
          continue;
        for (int v = 0; v < B; ++v)
        {
          x[v] = load(sx.data() + (m * B + v) * WIDTH);
          k[v] = load(sk.data() + (m * B + v) * WIDTH);
        }
        exact_trig<B>(st, m, cm, sm);
        if (column<B, G>(st, pv, opts, m, x, k, lump))
          combine<B, G>(m, lump, cm, sm, f);
      }
// Reduce this thread's accumulators into out. The critical section serialises the short final sum.
#pragma omp critical(gsh_mpar_reduce)
      {
        alignas(64) double buf[WIDTH];
        for (int v = 0; v < B; ++v)
          for (int nf = 0; nf < NF; ++nf)
          {
            store(buf, f[NF * v + nf]);
            for (int w = 0; w < WIDTH; ++w)
              out[nf * lanes + v * WIDTH + w] += buf[w];
          }
      }
    }
  }

  // ---- Grids: latitude lumps, with the even and odd (n - m) terms kept apart. ----
  // Layout lump[((m * 2 + parity) * 2P + l) * lanes + lane], already on absolute scale (R^-k applied).
  // Adds the term of degree N to the lumps of one parity class. y = q^N P̄_NM and ym1 = q^(N-1) P̄_(N-1)M are the
  // current and previous values; the derivatives follow Fukushima (2012, II) with fnm = f(N, m) = e_nm:
  //   dy = fnm (q/u) ym1 - N (t/u) y,    ddy = (t/u) dy + (m²/u² - N(N+1)) y.
  template <int G>
  GSH_INLINE void add_term(V *L, V C, V S, V y, V ym1, double N, double fnm, V qu, V tu, V ru2, V m2)
  {
    const V yC = mul(C, y), yS = mul(S, y);
    L[0] = add(L[0], yC);
    L[1] = add(L[1], yS);
    if constexpr (G >= 1)
    {
      const V w1 = set1(N + 1.0);
      const V dy = fnmadd(mul(set1(N), tu), y, mul(mul(set1(fnm), qu), ym1));
      L[2] = fmadd(w1, yC, L[2]);
      L[3] = fmadd(w1, yS, L[3]);
      L[4] = fmadd(C, dy, L[4]);
      L[5] = fmadd(S, dy, L[5]);
      if constexpr (G == 2)
      {
        const V w2 = set1(N + 2.0), nn1 = set1(N * (N + 1.0));
        const V ddy = fmadd(tu, dy, mul(fmsub(m2, ru2, nn1), y));
        const V wpr = mul(w1, w2);
        L[6] = fmadd(wpr, yC, L[6]);
        L[7] = fmadd(wpr, yS, L[7]);
        L[8] = fmadd(mul(w1, C), dy, L[8]);
        L[9] = fmadd(mul(w1, S), dy, L[9]);
        L[10] = fmadd(C, ddy, L[10]);
        L[11] = fmadd(S, ddy, L[11]);
      }
    }
  }

  // Grid lumps for WIDTH latitudes at once (one vector, B = 1). The longitude is not involved: the lumps are
  // combined with cos(mλ) and sin(mλ) later, in grid.cpp. Orders beyond the polar truncation are written as zeros
  // and the loop stops. Layout: see the section comment above.
  template <int G>
  void grid_lumps(const PrepView &pv, const TileOpts &opts, int, int, const double *lat, const double *r, double *out)
  {
    constexpr int B = 1;
    constexpr int P = LUMP_PAIRS<G>;
    const int nmax = pv.nmax, lanes = WIDTH;
    const int per_m = 2 * 2 * P * lanes;
    const double lon0[WIDTH] = {};
    TileState<B> st;
    init_state<B>(st, lat, lon0, r, pv.R, G);

    const V zero = set1(0.0), one = set1(1.0), two = set1(2.0);
    const V big = set1(R_BIG), rinv = set1(R_INV), negl = set1(NEGLIGIBLE_K1);
    V x[B] = {one}, k0[B] = {zero};

    for (int m = 0; m <= nmax; ++m)
    {
      double *dst = out + static_cast<std::int64_t>(m) * per_m;
      if (polar_skip(opts.polar_threshold, m, nmax, st.polar_u_max))
      {
        std::memset(dst, 0, sizeof(double) * static_cast<std::size_t>(nmax - m + 1) * per_m);
        break;
      }
      if (m > 0)
        next_sectorial<B>(st, pv.d[m], x, k0);

      const std::int64_t off = col_offset(m, nmax);
      const double *ab = pv.ab + 2 * off, *cs = pv.cs + 2 * off;
      const double *ff = G > 0 ? pv.f + off : nullptr;
      const int len = nmax - m;
      const V m2 = set1(static_cast<double>(m) * m);

      alignas(64) V L[2][2 * P];
      for (int p = 0; p < 2; ++p)
        for (int l = 0; l < 2 * P; ++l)
          L[p][l] = zero;
      V y2 = x[0], k = k0[0], y1 = zero;

      add_term<G>(L[0], set1(cs[0]), set1(cs[1]), y2, zero, m, 0.0, st.qu[0], st.tu[0], st.ru2[0], m2);
      if (len > 0)
      {
        y1 = mul(mul(set1(ab[2]), st.tq[0]), y2);
        add_term<G>(L[1], set1(cs[2]), set1(cs[3]), y1, y2, m + 1.0, G > 0 ? ff[1] : 0.0, st.qu[0], st.tu[0], st.ru2[0], m2);
      }
      int j = 2;
      while (j <= len)
      {
        const int jend = (j + RESCALE_EVERY - 1 < len) ? j + RESCALE_EVERY - 1 : len;
        bool accumulate = !opts.skip_forbidden;
        if (!accumulate)
          accumulate = any(mor(lt(k, one), mand(lt(k, two), mor(gt(abs(y1), negl), gt(abs(y2), negl)))));
        if (!accumulate)
          for (int i = j; i <= jend; ++i)
          {
            const V y = fnmadd(mul(set1(ab[2 * i + 1]), st.q2[0]), y2, mul(mul(set1(ab[2 * i]), st.tq[0]), y1));
            y2 = y1;
            y1 = y;
          }
        else
          for (int i = j; i <= jend; ++i)
          {
            const V y = fnmadd(mul(set1(ab[2 * i + 1]), st.q2[0]), y2, mul(mul(set1(ab[2 * i]), st.tq[0]), y1));
            add_term<G>(L[i & 1], set1(cs[2 * i]), set1(cs[2 * i + 1]), y, y1, static_cast<double>(m + i), G > 0 ? ff[i] : 0.0, st.qu[0], st.tu[0], st.ru2[0], m2);
            y2 = y1;
            y1 = y;
          }
        const M mr = mand(gt(abs(y1), big), gt(k, zero));
        if (any(mr))
        {
          y1 = blend(mr, y1, mul(y1, rinv));
          y2 = blend(mr, y2, mul(y2, rinv));
          k = blend(mr, k, sub(k, one));
          for (int p = 0; p < 2; ++p)
            for (int l = 0; l < 2 * P; ++l)
              L[p][l] = blend(mr, L[p][l], mul(L[p][l], rinv));
        }
        j = jend + 1;
      }
      const V scale = blend(lt(k, one), blend(lt(k, two), zero, rinv), one);
      for (int p = 0; p < 2; ++p)
        for (int l = 0; l < 2 * P; ++l)
          store(dst + (p * 2 * P + l) * lanes, mul(L[p][l], scale));
    }
  }

  // Dispatch helpers: map the vector-count index bi (0..4 for B = 1, 2, 4, 8, 16) to the tile function of gradient G.
  template <int G>
  constexpr kernels::TileFn serial_for_b(int bi)
  {
    switch (bi)
    {
    case 0:
      return &tile_serial<1, G>;
    case 1:
      return &tile_serial<2, G>;
    case 2:
      return &tile_serial<4, G>;
    case 3:
      return &tile_serial<8, G>;
    default:
      return &tile_serial<16, G>;
    }
  }
  template <int G>
  constexpr kernels::TileFn mpar_for_b(int bi)
  {
    switch (bi)
    {
    case 0:
      return &tile_mparallel<1, G>;
    case 1:
      return &tile_mparallel<2, G>;
    case 2:
      return &tile_mparallel<4, G>;
    case 3:
      return &tile_mparallel<8, G>;
    default:
      return &tile_mparallel<16, G>;
    }
  }

  // Fills the function-pointer table of this ISA: tile[G][bi], tile_mpar[G][bi] and grid[G].
  inline kernels::Table make_table(const char *name)
  {
    kernels::Table t{};
    t.name = name;
    t.width = WIDTH;
    for (int bi = 0; bi < kernels::NUM_B; ++bi)
    {
      t.tile[0][bi] = serial_for_b<0>(bi);
      t.tile[1][bi] = serial_for_b<1>(bi);
      t.tile[2][bi] = serial_for_b<2>(bi);
      t.tile_mpar[0][bi] = mpar_for_b<0>(bi);
      t.tile_mpar[1][bi] = mpar_for_b<1>(bi);
      t.tile_mpar[2][bi] = mpar_for_b<2>(bi);
    }
    t.grid[0] = &grid_lumps<0>;
    t.grid[1] = &grid_lumps<1>;
    t.grid[2] = &grid_lumps<2>;
    return t;
  }
}
