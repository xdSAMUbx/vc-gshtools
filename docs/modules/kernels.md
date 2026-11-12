# kernels: CPU SIMD kernel

> **Note on the examples.** The snippets below are illustrative fragments of internal headers (`namespace gsh::kernels`, ISA-specific templates). They are not self-contained translation units and are not compile-checked; the public-API examples in the other module docs are.


Internal API in `src/cpp/gsh_legendre/kernels/`. The public headers do not include it.

| File | Contents |
|---|---|
| `kernels.hpp` | Interface with the scheduler: views, fields and function table |
| `simd.hpp` | Lane type `V` for each ISA |
| `column_kernel.inl` | Algorithm, templated on `B` and `G`; the lane type `V` comes from `simd.hpp` |
| `kernels_{scalar,avx2,avx512}.cpp` | One compilation unit per ISA |

Each `.cpp` defines `GSH_ISA_*`, includes `simd.hpp` and `column_kernel.inl`, and exports `table_<isa>()`. The binary contains the three versions and `shs.cpp` selects one at run time.

## Algorithm

The recurrence is fused with the synthesis, without storing P:

```text
y_n = q^n P_nm
y_n = a (t q) y_{n-1} − b q² y_{n-2}
lumps: a = Σ w y C,  b = Σ w y S     per column m
sum:   field += a cos(mλ) + b sin(mλ)
```

**Exponent per lane.** The real value is `y · R^-k`, with `R = 2^500` (`R_BIG`, `R_INV` in `core.hpp`). Every `RESCALE_EVERY` (32) recurrence steps, a block boundary checks whether `|y| > R` and, if so, rescales mantissa and lumps and reduces k. At the end of the column: `k = 0` changes nothing, `k = 1` applies `R^-1`, and `k ≥ 2` gives 0.

**Derivatives per column.** As in CHarm, from `y` and `y_{n-1}`:

```text
dy  = f q/u · y_{n-1} − n t/u · y
ddy = t/u · dy + (m²/u² − n(n+1)) y
```

They are singular when `u → 0`, which is why near-polar points are routed to the full-triangle oracle (see [Polar handling](#polar-handling)). The error at `u = 1e-2` and `u = 1e-3` is not verified in this repository: the figures 1e-13 and 6e-11 appear only in `context.md`, not in a test or a code comment.

## Interface (`kernels.hpp`)

| Element | Meaning |
|---|---|
| `PrepView` | Pointers to `ab`, `cs`, `f`, `d`, plus `nmax` and `R` |
| `TileOpts` | `polar_threshold` (< 0 disables it), `skip_forbidden`, `rotate_trig` |
| `Field` | Order of the fields: `F_S0 … F_SLL` |
| `fields_for_grad(g)` | 1, 4 or 10 fields |
| `grid_groups(g)` | 1, 3 or 6 lump groups on grids |
| `TileFn` | Computes one tile. Output: `out[field * lanes + lane]` |
| `GridColumnFn` | Lumps for one batch of latitudes |
| `Table` | `name`, `width`, `tile[grad][bi]`, `tile_mpar[grad][bi]`, `grid[grad]`, with `bi = log2(vectors per tile)` |
| `table_scalar()`, `table_avx2()`, `table_avx512()` | Table for each ISA |

**Example** (`PrepView`)

```cpp
// Internal: built inside shs.cpp; normally reached through gsh::synth_points (src/cpp/gsh_expansion/shs.hpp).
#include "kernels.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = Model::zeros(3.986004415e14, 6378136.3, 0); // GM [m^3 s^-2], R [m], nmax
model.C[idx_tri(0, 0, 0)] = 1.0;                             // only C00 = 1
Prepared prep = prepare(model);
kernels::PrepView pv{prep.nmax, prep.rec->ab.data(), prep.cs.data(),
                     prep.rec->f.data(), prep.rec->d.data(), prep.R};
// pv.nmax == 0, pv.R == 6378136.3, pv.cs[0] == 1.0 (C00), pv.cs[1] == 0.0 (S00)
```

**Example** (`TileOpts`)

```cpp
// Internal: built inside shs.cpp from SynthOptions; normally reached through gsh::synth_points.
kernels::TileOpts exact{-1.0, false, false};   // no polar cut, no block skipping, cos/sin computed per order
kernels::TileOpts defaults{100.0, true, true}; // SynthOptions defaults at nmax = 2190: max(100, 0.01 * 2190) = 100
// A negative polar_threshold disables the polar cut, so polar_skip() below never returns true.
```

**Example** (`Field`)

```cpp
// Internal: the field order of the raw sums; normally reached through gsh::synth_points (RawFields::f is indexed by it).
static_assert(kernels::F_S0 == 0 && kernels::F_SL == 3 && kernels::F_SLL == 9 && kernels::F_COUNT == 10);
// With grad = 1 synth_points returns 4 vectors (F_S0..F_SL); with grad = 2 all 10 are filled.
// double dlambda = raw.f[kernels::F_SL][i];   // sum of q^n A^lambda_n for point i
```

**Example** (`fields_for_grad`)

```cpp
// Internal: sizes the output of synth_points; normally reached through gsh::synth_points.
static_assert(kernels::fields_for_grad(0) == 1);
static_assert(kernels::fields_for_grad(1) == 4);
static_assert(kernels::fields_for_grad(2) == 10);
```

**Example** (`grid_groups`)

```cpp
// Internal: sizes the grid lump buffer in grid.cpp; normally reached through the grid path of gsh::synth_points / gsh::compute_functionals.
static_assert(kernels::grid_groups(0) == 1 && kernels::grid_groups(1) == 3 && kernels::grid_groups(2) == 6);
// Lump buffer length: (nmax + 1) * 2 * 2 * grid_groups(grad) * lanes doubles.
// For nmax = 0, grad = 0 and the scalar ISA (lanes = 1): 4 doubles.
```

**Example** (`TileFn`)

```cpp
// Internal: TileFn values come from Table::tile; normally reached through gsh::synth_points (shs.cpp run_tile).
#include "kernels.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = Model::zeros(3.986004415e14, 6378136.3, 0);
model.C[idx_tri(0, 0, 0)] = 1.0;
Prepared prep = prepare(model);
kernels::PrepView pv{prep.nmax, prep.rec->ab.data(), prep.cs.data(), prep.rec->f.data(), prep.rec->d.data(), prep.R};
kernels::TileOpts opts{-1.0, true, true};

const kernels::Table &t = kernels::table_scalar();  // t.width == 1: one lane per point
kernels::TileFn fn = t.tile[0][0];                  // grad 0, bi = 0 -> one vector per tile
double lat[1] = {deg2rad(45.0)}, lon[1] = {deg2rad(30.0)}, r[1] = {6371000.0};
double out[1];                                      // out[field * lanes + lane]: one field, one lane
fn(pv, opts, lat, lon, r, out);
// out[0] == 1.0: with only C00 the value sum is C00 * q^0 * P00 = 1 at every point.
```

**Example** (`GridColumnFn`)

```cpp
// Internal: GridColumnFn values come from Table::grid; normally reached through the grid path (grid.cpp) of gsh::synth_points / gsh::compute_functionals.
#include "kernels.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = Model::zeros(3.986004415e14, 6378136.3, 0);
model.C[idx_tri(0, 0, 0)] = 1.0;
Prepared prep = prepare(model);
kernels::PrepView pv{prep.nmax, prep.rec->ab.data(), prep.cs.data(), prep.rec->f.data(), prep.rec->d.data(), prep.R};
kernels::TileOpts opts{-1.0, true, true};

kernels::GridColumnFn gfn = kernels::table_scalar().grid[0];  // grad 0, one vector of latitudes
double lat[1] = {deg2rad(45.0)}, r[1] = {6371000.0};
double lump[4];                                               // (nmax+1) * 2 * 2 * P * lanes = 4 doubles
gfn(pv, opts, 0, 0, lat, r, lump);                            // the two int arguments are unused
// lump == {1, 0, 0, 0}: m = 0, even parity: sum y C = 1 and sum y S = 0; odd parity is empty.
```

**Example** (`Table`)

```cpp
// Internal: each Table is built once per ISA by make_table(); normally reached through gsh::synth_points (shs.cpp select()).
const kernels::Table &s = kernels::table_scalar();
// s.name == "scalar", s.width == 1
// s.tile[grad][bi]: bi = log2(B), B = vectors per tile, so lanes = B * width. With SynthOptions::tile_vectors = 0 (automatic),
// shs.cpp uses B = 4 (bi = 2) for grad 0, B = 2 (bi = 1) for grad 1 and B = 1 (bi = 0) for grad 2.
// s.tile_mpar[grad][bi] has the same indexing and is used when there are few tiles; s.grid[grad] serves the grids (B = 1 always).
```

**Example** (`table_scalar()`, `table_avx2()`, `table_avx512()`)

```cpp
// Internal: one table per ISA; normally reached through gsh::synth_points with SynthOptions::isa (shs.cpp select()).
#include <string>
#include "kernels.hpp"
#include "shs.hpp"
using namespace gsh;
const kernels::Table &scalar = kernels::table_scalar();    // always available, width 1
std::string isa = active_isa_name(Isa::Avx2);              // "avx2" (width 4), or throws if the CPU lacks AVX2+FMA
// const kernels::Table &avx = kernels::table_avx2();      // only after the check above: on a CPU without AVX2 the kernels fault
// const kernels::Table &x512 = kernels::table_avx512();   // width 8; needs AVX-512F/DQ/VL and OS ZMM state
```

## Lane (`simd.hpp`)

`V` and the mask `M` have the same interface in the three ISAs. `WIDTH` is 8, 4 or 1.

| Function | Operation |
|---|---|
| `set1`, `load`, `store` | Fill, load and store without requiring alignment |
| `add`, `sub`, `mul`, `div`, `abs` | Per-lane arithmetic |
| `fmadd(a,b,c)` | `a*b + c` |
| `fnmadd(a,b,c)` | `c − a*b` |
| `fmsub(a,b,c)` | `a*b − c` |
| `gt`, `lt`, `ne` | Comparisons that return a mask |
| `mand`, `mor`, `any` | Operations on masks |
| `blend(m,f,t)` | `m ? t : f` |

Not callable on their own: the lane names exist only in a translation unit compiled for one ISA (one `GSH_ISA_*` macro defined, matching compile flags). The examples below use the scalar ISA.

**Example** (lane setup)

```cpp
// Internal: needs a TU that defines one GSH_ISA_* macro; normally reached through gsh::synth_points (the kernels call these).
#define GSH_ISA_SCALAR   // or GSH_ISA_AVX2 (gsh::isa_avx2, compile with -mavx2 -mfma) or GSH_ISA_AVX512 (gsh::isa_avx512)
#include "simd.hpp"
using namespace gsh::isa_scalar;
static_assert(WIDTH == 1);  // 4 for AVX2 and 8 for AVX-512
```

**Example** (`set1`, `load`, `store`)

```cpp
// Internal: lane primitive, needs the lane setup above; normally reached through gsh::synth_points.
double buf[WIDTH];                 // WIDTH = 1 here; 4 or 8 doubles in the SIMD ISAs
store(buf, set1(2.5));             // buf[0] == 2.5: every lane gets the broadcast value
V b = load(buf);                   // unaligned load: b.v == 2.5 in the scalar ISA
```

**Example** (`add`, `sub`, `mul`, `div`, `abs`)

```cpp
// Internal: lane primitive, needs the lane setup above; normally reached through gsh::synth_points.
V x = set1(-3.0), y = set1(2.0);
V s = add(x, y);   // -1.0
V d = sub(x, y);   // -5.0
V p = mul(x, y);   // -6.0
V q = div(x, y);   // -1.5
V ax = abs(x);     //  3.0
```

**Example** (`fmadd`, `fnmadd`, `fmsub`)

```cpp
// Internal: lane primitive, needs the lane setup above; normally reached through gsh::synth_points.
V a = set1(2.0), b = set1(3.0), c = set1(1.0);
V f1 = fmadd(a, b, c);   // a*b + c =  7.0
V f2 = fnmadd(a, b, c);  // c - a*b = -5.0
V f3 = fmsub(a, b, c);   // a*b - c =  5.0
```

**Example** (`gt`, `lt`, `ne`)

```cpp
// Internal: lane primitive returning a mask M, needs the lane setup above; normally reached through gsh::synth_points.
V a = set1(2.0), b = set1(3.0);
M g = gt(a, b);   // false (2 > 3)
M l = lt(a, b);   // true  (2 < 3)
M n = ne(a, a);   // false (a is equal to itself)
```

**Example** (`mand`, `mor`, `any`)

```cpp
// Internal: mask logic, needs the lane setup above; normally reached through gsh::synth_points.
M t = lt(set1(1.0), set1(2.0));   // true
M f = gt(set1(1.0), set1(2.0));   // false
bool both = any(mand(t, f));      // false
bool either = any(mor(t, f));     // true
```

**Example** (`blend`)

```cpp
// Internal: lane primitive, needs the lane setup above; normally reached through gsh::synth_points.
V lo = set1(1.0), hi = set1(9.0);
M pick_hi = gt(hi, lo);           // true in every lane
V r = blend(pick_hi, lo, hi);     // blend(m, f, t) = m ? t : f, so the result is hi (9.0)
```

## Functions of `column_kernel.inl`

### Constants
`RESCALE_EVERY = 32`, `TRIG_RESET = 64` (exact reset of the rotation), `NEGLIGIBLE_K1 = 2^300` (with k = 1, a `|y|` smaller than this implies a value below 2^-200). `LUMP_PAIRS<G>` is 1, 3 or 6.

Not callable on their own: every name in this section lives in `gsh::isa_<isa>` and exists only in a translation unit that defines one `GSH_ISA_*` macro and includes `column_kernel.inl`. From ordinary code, reach the same instantiations through `Table::tile`, `Table::tile_mpar` and `Table::grid` (see the `Table` example above).

**Example** (constants, and the setup shared by the examples below)

```cpp
// Internal: needs this kernel TU prelude; normally reached through gsh::synth_points / gsh::compute_functionals (src/cpp/gsh_expansion/shs.hpp).
#define GSH_ISA_SCALAR   // same TU structure as kernels_scalar.cpp
#include "simd.hpp"
#include "column_kernel.inl"
#include "shs.hpp"
using namespace gsh::isa_scalar;
using gsh::Model; using gsh::Prepared; using gsh::prepare; using gsh::idx_tri; using gsh::deg2rad;

static_assert(RESCALE_EVERY == 32 && TRIG_RESET == 64);
static_assert(NEGLIGIBLE_K1 == 0x1p300);   // k = 1 means value = y * 2^-500, so |y| < 2^300 gives |value| < 2^-200
static_assert(LUMP_PAIRS<0> == 1 && LUMP_PAIRS<1> == 3 && LUMP_PAIRS<2> == 6);

// Setup used by the examples below: degree 0, only C00 = 1.
Model model = Model::zeros(3.986004415e14, 6378136.3, 0);
model.C[idx_tri(0, 0, 0)] = 1.0;
Prepared prep = prepare(model);
PrepView pv{prep.nmax, prep.rec->ab.data(), prep.cs.data(), prep.rec->f.data(), prep.rec->d.data(), prep.R};
TileOpts opts{-1.0, true, true};
```

### `TileState<B>` and `init_state<B>(st, lat, lon, r, R, G)`
Precomputes per lane `t, u, q, tq, q2, uq` and, with derivatives, `q/u, t/u, 1/u²`. Also `cos λ`, `sin λ` and `polar_u_max`, the largest `u` in the tile. Handles exact poles without the residual of `cos(π/2)`.

**Example** (`TileState`, `init_state`)

```cpp
// Internal: needs the kernel TU prelude and setup from the constants example; normally reached through gsh::synth_points.
TileState<1> st;
double lat[1] = {deg2rad(45.0)}, lon[1] = {deg2rad(30.0)}, r[1] = {6371000.0};
init_state<1>(st, lat, lon, r, pv.R, 0);   // B = 1 vector, G = 0 (the q/u, t/u, 1/u^2 fields are filled only for G > 0)
// st.t[0].v ~ 0.70711 (sin 45), st.u[0].v ~ 0.70711, st.q[0].v = R/r ~ 1.00112, st.polar_u_max ~ 0.70711.
// A pole (lat = 90 deg) gives t = 1 and u = 0 exactly.
```

### `next_sectorial<B>(st, dm, x, k)`
Advances `x ← d_m · u q · x`. If `|x| < R^-1`, it multiplies by `R` and increments k for that lane.

**Example** (`next_sectorial`)

```cpp
// Internal: needs the kernel TU prelude and setup from the constants example; normally reached through gsh::synth_points.
TileState<1> st;
double lat[1] = {deg2rad(45.0)}, lon[1] = {0.0}, r[1] = {6371000.0};
init_state<1>(st, lat, lon, r, pv.R, 0);
V x[1] = {set1(0x1p-501)}, k[1] = {set1(0.0)};
next_sectorial<1>(st, 1.0, x, k);   // x <- d_m (u q) x with d_m = 1 and u q ~ 0.70790
// 0.7079 * 2^-501 is below R^-1 = 2^-500, so it is multiplied by R once:
// x[0].v ~ 0.354 and k[0].v == 1 (the value is x * R^-1).
```

### `column<B, G>(st, pv, opts, m, x0, k0, lump)`
Traverses the whole column m and leaves the lumps on an absolute scale. Returns `false` when every lane ends with k ≥ 2 (the column is negligible for all lanes), in which case the caller skips `combine`.

- Terms `n = m` and `n = m+1` are handled separately, because `b = 0`.
- Blocks of `RESCALE_EVERY` terms. If no lane can contribute, it only runs the recurrence and avoids reading C and S. This is skipping the forbidden zone.
- A G0-only loop, and one loop shared by G1 and G2 (G2 adds the second-derivative lumps under `if constexpr`).

**Example** (`column`)

```cpp
// Internal: needs the kernel TU prelude and setup from the constants example; normally reached through gsh::synth_points.
TileState<1> st;
double lat[1] = {deg2rad(45.0)}, lon[1] = {0.0}, r[1] = {6371000.0};
init_state<1>(st, lat, lon, r, pv.R, 0);
V x[1] = {set1(1.0)}, k[1] = {set1(0.0)};   // sectorial value and exponent for m = 0
V lump[2];                                  // 2 * P * B = 2 entries for G = 0
bool sig = column<1, 0>(st, pv, opts, 0, x, k, lump);
// lump[0].v == 1.0 (C00 * y_0), lump[1].v == 0.0 (S00), sig == true (the lane has k = 0).
```

### `combine<B, G>(m, lump, cm, sm, f)`
Adds the contribution of one column to the fields with `cos/sin(mλ)`. `SL`, `SRL`, `SPL` and `SLL` come from the same lumps with a factor m or m².

**Example** (`combine`)

```cpp
// Internal: needs the kernel TU prelude and setup from the constants example; normally reached through gsh::synth_points.
V lump[2] = {set1(2.0), set1(3.0)};          // a = sum y C and b = sum y S for column m = 1
V cm[1] = {set1(1.0)}, sm[1] = {set1(0.0)};  // cos(1 * lambda), sin(1 * lambda) with lambda = 0
V f[1] = {set1(0.0)};                        // F_S0 accumulator (one field for G = 0)
combine<1, 0>(1, lump, cm, sm, f);           // F_S0 += a cos(lambda) + b sin(lambda)
// f[0].v == 2.0. At lambda = 90 deg (cm = 0, sm = 1) the same call would add b = 3 instead.
```

### `exact_trig<B>(st, m, cm, sm)`
Direct `cos(mλ)` and `sin(mλ)`, lane by lane.

**Example** (`exact_trig`)

```cpp
// Internal: needs the kernel TU prelude and setup from the constants example; normally reached through gsh::synth_points.
TileState<1> st;
double lat[1] = {deg2rad(45.0)}, lon[1] = {deg2rad(30.0)}, r[1] = {6371000.0};
init_state<1>(st, lat, lon, r, pv.R, 0);
V cm[1], sm[1];
exact_trig<1>(st, 2, cm, sm);   // cos(2 * 30 deg) and sin(2 * 30 deg), lane by lane
// cm[0].v ~ 0.5 and sm[0].v ~ 0.86603.
```

### `polar_skip(threshold, m, nmax, umax)`
True if `threshold >= 0` and `m − nmax·umax > threshold`, where `umax` is `TileState::polar_u_max` (the largest cos φ in the tile), so the test is per tile. The threshold is passed in by the scheduler; this function does not compute it (see [Polar handling](#polar-handling)). The condition is monotonic in m, so the loop can be cut off.

**Example** (`polar_skip`)

```cpp
// Internal: inline helper in gsh::isa_*, needs the kernel TU prelude; normally reached through gsh::synth_points.
bool a = polar_skip(100.0, 121, 2190, 0.01);  // 121 - 2190 * 0.01 = 99.1 -> false, column kept
bool b = polar_skip(100.0, 122, 2190, 0.01);  // 122 - 21.9 = 100.1 -> true, the m loop stops here
bool c = polar_skip(-1.0, 122, 2190, 0.01);   // negative threshold disables the cut -> false
```

### `tile_serial<B, G>(pv, opts, lat, lon, r, out)`
Full tile on one thread. Sectorial terms and trigonometry advance on the fly, without being stored. Rotation is only used in G0: in G1 and G2 the weights m and m² amplify its error. The claim that rotation gives no gain in G1 and G2 is not backed by an archived measurement in this repository (see the TODO in `column_kernel.inl`).

**Example** (`tile_serial`)

```cpp
// Internal: tile_serial<B, G> is named only inside the kernel TU prelude; normally reached through gsh::synth_points via Table::tile.
double lat[1] = {deg2rad(45.0)}, lon[1] = {deg2rad(30.0)}, r[1] = {6371000.0}, out[1];
tile_serial<1, 0>(pv, opts, lat, lon, r, out);                     // direct call, B = 1, G = 0
gsh::kernels::table_scalar().tile[0][0](pv, opts, lat, lon, r, out); // same instantiation through the table
// out[0] == 1.0 (only C00 is nonzero).
```

### `tile_mparallel<B, G>(pv, opts, lat, lon, r, out)`
For few points. It first stores all sectorial terms and then distributes the orders m among threads, with a reduction in a critical section.

**Example** (`tile_mparallel`)

```cpp
// Internal: tile_mparallel<B, G> is named only inside the kernel TU prelude; normally reached through gsh::synth_points via Table::tile_mpar.
double lat[1] = {deg2rad(45.0)}, lon[1] = {deg2rad(30.0)}, r[1] = {6371000.0}, out[1];
tile_mparallel<1, 0>(pv, opts, lat, lon, r, out);                          // OpenMP over m; serial if built without OpenMP
gsh::kernels::table_scalar().tile_mpar[0][0](pv, opts, lat, lon, r, out);  // same kernel through the table
// out[0] == 1.0. With more orders in the sum, the reduction order can change the last bits, not the result.
```

### `add_term<G>(L, C, S, y, ym1, N, fnm, qu, tu, ru2, m2)`
Adds one term (n, m) to the lumps of one parity. Used only for grids.

**Example** (`add_term`)

```cpp
// Internal: add_term<G> is named only inside the kernel TU prelude; normally reached through the grid path of gsh::synth_points (grid.cpp).
V L[2] = {set1(0.0), set1(0.0)};                  // one parity class, G = 0: the a and b lumps
V zero = set1(0.0), one = set1(1.0);
add_term<0>(L, set1(1.0), set1(0.0), one, zero,   // C, S, y, ym1
            0.0, 0.0, zero, zero, zero, zero);    // N, fnm, then qu, tu, ru2, m2 (unused for G = 0)
// L[0].v == 1.0 (C * y) and L[1].v == 0.0 (S * y).
```

### `grid_lumps<G>(pv, opts, _, _, lat, r, out)`
Lumps for `WIDTH` latitudes, separating even and odd terms in `(n − m)` for the ±φ symmetry. Layout:

```text
out[((m * 2 + parity) * 2P + l) * lanes + lane]
```

Each group is an (a, b) pair of lumps, so `l` runs over 2P entries. With `C` and `S` the coefficients of (n, m):

- group 0: a = Σ y C, b = Σ y S
- group 1: Σ (n+1) y C, Σ (n+1) y S
- group 2: Σ dy C, Σ dy S
- group 3: Σ (n+1)(n+2) y C, Σ (n+1)(n+2) y S
- group 4: Σ (n+1) dy C, Σ (n+1) dy S
- group 5: Σ ddy C, Σ ddy S

G0 uses group 0 only, G1 groups 0 to 2, G2 groups 0 to 5.

**Example** (`grid_lumps`)

```cpp
// Internal: grid_lumps<G> is named only inside the kernel TU prelude; normally reached through the grid path (grid.cpp) via Table::grid.
double lat[1] = {deg2rad(45.0)}, r[1] = {6371000.0}, lump[4];   // (nmax+1) * 2 parities * 2 * P * lanes = 4 doubles
grid_lumps<0>(pv, opts, 0, 0, lat, r, lump);                    // the two int parameters are unused
// lump == {1, 0, 0, 0}: m = 0, even parity: sum y C = 1 and sum y S = 0; odd parity is empty.
```

### `serial_for_b<G>(bi)`, `mpar_for_b<G>(bi)`, `make_table(name)`
Build the table of function pointers for each combination of B and G.

**Example** (`serial_for_b`, `mpar_for_b`, `make_table`)

```cpp
// Internal: named only inside the kernel TU prelude; normally reached through gsh::synth_points, which picks the entry from the vectors-per-tile count.
gsh::kernels::TileFn f4 = serial_for_b<1>(2);   // bi = 2 -> tile_serial<4, 1>: 4 vectors, value and first derivatives
gsh::kernels::TileFn p8 = mpar_for_b<2>(3);     // bi = 3 -> tile_mparallel<8, 2>: 8 vectors, second derivatives
gsh::kernels::Table t = make_table("scalar");   // same layout as table_scalar(), built again
// t.tile[1][2] is the same instantiation as f4, and t.grid[0] is grid_lumps<0>.
```

## Polar handling

Three different thresholds are involved. They are easy to confuse:

- **Truncation in m (`polar_threshold`).** `synth_points` (shs.cpp) and `synth_grid` (grid.cpp) compute `μ = max(100, 0.01 · nmax)` when `SynthOptions::polar_optimization` is true (default), and `-1` (disabled) otherwise. The value goes into `TileOpts::polar_threshold`. The kernels apply it through `polar_skip()` (column_kernel.inl) in `tile_serial`, `tile_mparallel` and `grid_lumps`: once `m − nmax · umax > μ`, the remaining orders are dropped (`grid_lumps` writes zeros for them). `umax` is the largest cos φ in the tile, so the test is per tile, not per point. The formula for μ has a TODO in shs.cpp: its justification is not archived in this repository.
- **Oracle for near-polar points (`pole_cos_threshold`).** `FunctionalRequest::pole_cos_threshold` (functionals.hpp, default 5e-3) is tested in `compute_functionals` (functionals.cpp). Points with cos φ below it are recomputed with `reference_point()` after the synthesis, which overwrites their raw sums. `synth_grid` (grid.cpp) uses the same value, hard-coded as 5e-3 and only for grad > 0: those rows go to the full-triangle oracle. The 5e-3 value has a TODO; no accuracy measurement is archived in this repository.
- **Not a kernel threshold.** The kernels never see `pole_cos_threshold`; they only see the μ of the first bullet.

## Putting it together

The public path ends in these kernels. `synth_points` (shs.cpp) calls `select(opts.isa)` to get a `Table`, takes `tile[grad][bi]` (or `tile_mpar[grad][bi]` when there are few tiles) and runs it on each tile of `lanes = vectors * width` points. `compute_functionals` (functionals.cpp) calls `synth_points` for all points and then replaces the points flagged by `pole_cos_threshold` with `reference_point()`; see [Polar handling](#polar-handling). The grid path (grid.cpp) reads `Table::grid` instead.

```cpp
// Public path: user code calls these; the kernels above run inside synth_points.
#include <iostream>
#include "functionals.hpp"   // includes shs.hpp, model.hpp and synth_types.hpp
using namespace gsh;

Model model = load_gfc("model.gfc", 360);                  // ICGEM file, truncated to degree 360
Prepared prep = prepare(model);                            // recurrence tables and packed C, S
Points pts = Points::from_spherical({45.0, -12.5}, {30.0, 100.0}, {6371000.0, 6378137.0});  // deg, deg, m

SynthOptions opts;
opts.isa = Isa::Auto;                                      // or Isa::Avx2, Isa::Avx512, Isa::Scalar
std::cout << active_isa_name(opts.isa) << "\n";           // "avx512", "avx2" or "scalar", as the CPU allows

// grad 1 uses 2 vectors per tile (lanes = 2 * width). With 2 points and several threads, synth_points
// takes the tile_mpar[1][1] path (tile_mparallel<2, 1>) instead of tile[1][1].
RawFields raw = synth_points(prep, pts, 1, opts);          // raw.f holds F_S0, F_SR, F_SP, F_SL
std::cout << raw.f[kernels::F_S0][0] << "\n";              // raw sum at point 0 (dimensionless)

FunctionalRequest req;
req.outputs = {Functional::V, Functional::T};
FunctionalResult res = compute_functionals(model, pts, req);
double V0 = res.values.at(Functional::V)[0];               // potential at point 0 (m^2 s^-2)
```
