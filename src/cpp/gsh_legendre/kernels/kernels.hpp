// kernels.hpp: internal interface between the synthesis schedulers (shs.cpp, grid.cpp) and the per-ISA
// column kernels (column_kernel.inl, compiled once for scalar, AVX2 and AVX-512).
//
// The kernels never see the scheduler's data layout directly: they receive a PrepView (coefficients, recurrence
// tables) and a TileOpts (optimisation switches), plus per-tile point arrays. The ISA is chosen at run time by
// shs.cpp (select()) from the Table returned by table_scalar(), table_avx2() or table_avx512().
#pragma once
#include <cstdint>

namespace gsh::kernels
{
  // Non-owning view of the prepared data for one synthesis call. All arrays are packed by column m, see
  // legendre.hpp for the triangle layout (idx_tri / col_offset).
  struct PrepView
  {
    int nmax;
    const double *ab; // [a, b] per (n, m): 2 * tri_size entries, interleaved as ab[2i], ab[2i+1]
    const double *cs; // [C, S] per (n, m): 2 * tri_size entries (already truncated to nmin..nmax in prepare())
    const double *f;  // f(n,m) for dP/dφ; only read when grad >= 1
    const double *d;  // sectorial scale factors d[m], nmax + 1 entries
    double R;         // model reference radius (m); the radial factor is q = R / r
  };

  // Per-call options, forwarded from SynthOptions (shs.hpp / synth_types.hpp).
  struct TileOpts
  {
    // Polar truncation threshold μ. Negative disables it; otherwise columns with m - nmax*u > μ are skipped
    // (Jekeli-type truncation, Jekeli, Lee and Kwon 2007).
    double polar_threshold;
    // When true, a block of recurrence steps is not accumulated where every lane's value is negligible
    // (|P| < 2^-200 in all lanes), which saves the coefficient work for those terms.
    bool skip_forbidden;
    // When true, cos(mλ) and sin(mλ) are advanced by angle-addition rotation (restarted periodically) instead
    // of being recomputed with std::cos / std::sin for every order.
    bool rotate_trig;
  };

  // Raw per-point sums, in this fixed order. Each field is a sum over (n, m) of coefficient x Legendre term, with no
  // GM/r^k factor (the functionals in functionals.cpp apply it). Definitions: Bucha and Janák (2019)
  // (A_n, A^φ_n, A^λ_n, ...). q_n = (R/r)^n is included in each sum.
  // grad 0 computes F_S0 only, grad 1 the first four, grad 2 all ten.
  enum Field
  {
    F_S0 = 0, // Σ q_n A_n                      (potential V or T)
    F_SR,     // Σ q_n (n+1) A_n               (radial derivative)
    F_SP,     // Σ q_n A^φ_n                    (latitude derivative)
    F_SL,     // Σ q_n A^λ_n                    (longitude derivative, includes the factor m)
    F_SRR,    // Σ q_n (n+1)(n+2) A_n          (second radial derivative)
    F_SRP,    // Σ q_n (n+1) A^φ_n
    F_SRL,    // Σ q_n (n+1) A^λ_n
    F_SPP,    // Σ q_n A^φφ_n                   (second latitude derivative)
    F_SPL,    // Σ q_n A^φλ_n
    F_SLL,    // Σ q_n A^λλ_n                   (second longitude derivative, includes m²)
    F_COUNT
  };
  // Number of fields produced for a gradient order: 0 -> 1 (value), 1 -> 4 (value + first derivatives), 2 -> 10.
  inline constexpr int fields_for_grad(int grad) { return grad == 0 ? 1 : (grad == 1 ? 4 : 10); }

  // Tile function: processes WIDTH * B points (lanes) at once. lat, lon (radians), r (metres) hold one value per
  // lane. Results go to out[field * lanes + lane] in the Field order above.
  using TileFn = void (*)(const PrepView &, const TileOpts &, const double *lat, const double *lon, const double *r, double *out);

  // Grid kernel: computes the latitude "lumps" of one tile of WIDTH latitudes (longitude is applied later, in grid.cpp).
  // Parameters 3 and 4 (the m range of the original signature) are unused by the current implementation.
  // The lumps are split into even and odd (n - m) terms, because the mirror latitude -φ flips their sign differently.
  // Layout: lump[((m * 2 + parity) * 2P + l) * WIDTH + lane], with P = grid_groups(grad) and l = a/b of each group.
  // Groups (grid_lumps in column_kernel.inl):
  //   0: Σ y C, Σ y S   1: (n+1) y   2: dy   3: (n+1)(n+2) y   4: (n+1) dy   5: ddy
  // where y = q^n P̄_nm(φ), dy = q^n dP̄/dφ and ddy = q^n d²P̄/dφ², and C, S are the coefficients of (n, m).
  using GridColumnFn = void (*)(const PrepView &, const TileOpts &, int, int, const double *lat, const double *r, double *lump);

  // Number of supported vector counts per tile, B = 1, 2, 4, 8, 16; table index bi = log2(B).
  inline constexpr int NUM_B = 5;
  // Function-pointer table for one ISA. Indexed by [grad][bi] with bi = log2(B), see make_table() in column_kernel.inl.
  struct Table
  {
    const char *name;           // "scalar", "avx2" or "avx512"; returned by active_isa_name()
    int width;                  // lanes per SIMD vector: 1, 4 or 8
    TileFn tile[3][NUM_B];      // serial over m; the caller parallelises over tiles
    TileFn tile_mpar[3][NUM_B]; // parallel over m inside the tile; used when there are few points
    GridColumnFn grid[3];       // one vector per call (B = 1)
  };

  // One table per ISA. Each is defined in its own translation unit (kernels_scalar.cpp, kernels_avx2.cpp,
  // kernels_avx512.cpp), compiled with the matching instruction set.
  const Table &table_scalar();
  const Table &table_avx2();
  const Table &table_avx512();

  // Number of lump groups for a gradient order: 1 (grad 0), 3 (grad 1) or 6 (grad 2). Each group has an a and a b part.
  inline constexpr int grid_groups(int grad) { return grad == 0 ? 1 : (grad == 1 ? 3 : 6); }
}
