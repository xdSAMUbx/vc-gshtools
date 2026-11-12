# grid: synthesis on grids

Header `src/cpp/gsh_expansion/grid.hpp`, implementation `src/cpp/gsh_expansion/grid.cpp`.

## Examples: setup

All examples assume:

```cpp
#include "grid.hpp"   // also brings in shs.hpp, model.hpp, synth_types.hpp and legendre.hpp
using namespace gsh;
```

Add `src/cpp/gsh_expansion`, `src/cpp/gsh_utils` and `src/cpp/gsh_legendre` to the include path. Examples that need a model use `Model::zeros` with only one coefficient set, so the expected values follow by hand from the formulas. Real models come from `load_gfc` (see "Putting it together").

## `struct Grid`

| Field | Meaning |
|---|---|
| `lat` | Latitudes in radians |
| `r` | One radius per latitude |
| `lon0`, `dlon` | First longitude and step, in radians |
| `nlon` | Number of longitudes |

`nlat()` returns the number of latitudes.

**Example**

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

Grid g;                            // all fields are public
g.lat = {0.5, 0.4};                // radians, north to south
g.r = {6378136.3, 6378136.3};      // metres, one per row
g.lon0 = 0.0;                      // first longitude, radians
g.dlon = PI / 180.0;               // 1 degree step, in radians
g.nlon = 3;                        // columns at 0, 1 and 2 degrees
std::printf("%lld rows\n", (long long)g.nlat());     // 2 rows
std::printf("full circle: %d\n", g.full_circle());   // 0: 3 degrees is not 2 pi
```

### Global grids

Same as CHarm, with `L = nmax + 1`, latitudes from north to south, and `λ_j = j Δλ` starting at 0.

| Function | Latitudes × longitudes | Nodes |
|---|---|---|
| `gauss_legendre(nmax, r)` | L × 2L | Legendre roots, with exact symmetry |
| `driscoll_healy1(nmax, r)` | 2L × 2L | Equally spaced; includes the north pole |
| `driscoll_healy2(nmax, r)` | 2L × 4L | Same, with twice as many longitudes |

**Example** (`gauss_legendre`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// nmax = 4 gives L = 5 latitudes (roots of P_5) and 2L = 10 longitudes.
const Grid g = Grid::gauss_legendre(4, 6378136.3);
std::printf("%lld x %lld\n", (long long)g.nlat(), (long long)g.nlon);      // 5 x 10
for (std::int64_t i = 0; i < g.nlat(); ++i)
  std::printf("lat %8.4f deg\n", rad2deg(g.lat[i]));  // 64.98, 32.58, 0, -32.58, -64.98 (approx.)
std::printf("dlon = %.1f deg\n", rad2deg(g.dlon));                          // 36.0
```

**Example** (`driscoll_healy1`, `driscoll_healy2`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// nmax = 4: L = 5, so 2L = 10 latitudes at 90, 72, ..., -72 deg (the row at 0 is exact).
const Grid dh1 = Grid::driscoll_healy1(4, 6378136.3);   // 10 x 10, dlon = 36 deg
const Grid dh2 = Grid::driscoll_healy2(4, 6378136.3);   // 10 x 20, dlon = 18 deg
std::printf("DH1: %lld x %lld\n", (long long)dh1.nlat(), (long long)dh1.nlon);   // 10 x 10
std::printf("DH2: %lld x %lld\n", (long long)dh2.nlat(), (long long)dh2.nlon);   // 10 x 20
std::printf("first row %.1f deg, row 5 %.1f deg\n", rad2deg(dh1.lat[0]), rad2deg(dh1.lat[5]));  // 90.0, 0.0
```

### `regular(lat_north, lat_south, nlat, lon_west, lon_east, nlon, r)`
Regular grid in degrees, with the end points included. It can be regional.

**Example**

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// Regional: 7 latitudes from 45 N to 30 N (step -2.5 deg), 11 longitudes from 250 to 260 deg E (step 1 deg).
const Grid reg = Grid::regular(45.0, 30.0, 7, 250.0, 260.0, 11, 6378136.3 + 5000.0);
std::printf("%lld x %lld\n", (long long)reg.nlat(), (long long)reg.nlon);            // 7 x 11
std::printf("lat[1] = %.1f deg\n", rad2deg(reg.lat[1]));                              // 42.5
std::printf("lon0 = %.1f, dlon = %.1f deg\n", rad2deg(reg.lon0), rad2deg(reg.dlon));  // 250.0, 1.0

// Global: longitudes 0 to 359 deg with 360 columns give dlon = 1 deg, so nlon * dlon = 2 pi (FFT path).
// Using 0 to 360 with 360 columns would give dlon = 360/359 deg and the direct path.
const Grid global = Grid::regular(89.0, -89.0, 179, 0.0, 359.0, 360, 6378136.3);
std::printf("global full circle: %s\n", global.full_circle() ? "yes" : "no");      // yes
```

### `bool full_circle() const`
True if `nlon·dlon = 2π`, with a tolerance of 1e-12. Decides whether an FFT is used.

**Example**

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

const Grid gl = Grid::gauss_legendre(40, 6378136.3);                          // nlon = 82
const Grid reg = Grid::regular(45.0, 30.0, 7, 250.0, 260.0, 11, 6378136.3);   // 11 deg of longitude
std::printf("gauss:    %s\n", gl.full_circle() ? "true" : "false");   // true: 82 * 2 pi / 82 = 2 pi
std::printf("regional: %s\n", reg.full_circle() ? "true" : "false");  // false: 11 deg != 2 pi
```

## `RawFields synth_grid(prep, grid, grad, opts = {}, timing = nullptr)`

Computes the [raw sums](../README.md#sumas-crudas) on the grid. Each field has `nlat·nlon` values, stored by rows: `f[k][i*nlon + j]`.

**CPU steps:**

1. Pairs latitudes `±φ` with the same radius. One task serves both rows.
2. With `grad > 0`, separates the rows with `cos φ < 5e-3`, which are computed with the full triangle.
3. Computes lumps in batches of `WIDTH` latitudes (the SIMD width of the ISA in use: 1, 4 or 8) with `grid_lumps`, one thread per batch.
4. Converts the lumps into rows:
   - **Full circle:** FFTW `c2r` FFT. The plan is shared and each thread runs it with its own buffers. If `m` exceeds `nlon/2`, the frequency is folded, so few longitudes are supported.
   - **Regional:** direct summation by rotation along the parallel.

**Errors:** `invalid_argument` if `grad` is not 0, 1 or 2, or if `r` does not have one value per latitude; `runtime_error` if FFTW cannot create the plan.

**Example**

```cpp
#include <cstdio>
#include <stdexcept>
#include "grid.hpp"
using namespace gsh;

Model model = Model::zeros(3.986004415e14, 6378136.3, 4);  // GM [m^3/s^2], R [m], nmax 4
model.C[idx_tri(0, 0, 4)] = 1.0;                            // only C00 = 1 (fully normalized)
const auto prep = prepare(model);                           // nmin = 0, nmax = 4
const Grid g = Grid::gauss_legendre(4, model.R * 1.001);    // 5 x 10 nodes
const RawFields raw = synth_grid(prep, g, 0);
std::printf("fields %zu, values per field %zu\n", raw.f.size(), raw.f[0].size());  // 1, 50
// C00 is the only term: no latitude or longitude dependence, so every entry of a row is equal.
// Entry (i, j) is raw.f[0][i * g.nlon + j]; row 0 is the northernmost latitude.
const double v = raw.f[0][0 * g.nlon + 3];

SynthOptions opts;
opts.threads = 4;                        // OpenMP threads (0 = maximum)
const RawFields r1 = synth_grid(prep, g, 1, opts);   // grad 1: 4 fields
const RawFields r2 = synth_grid(prep, g, 2, opts);   // grad 2: 10 fields

try { synth_grid(prep, g, 3); }          // grad outside 0..2
catch (const std::invalid_argument &e) { std::printf("%s\n", e.what()); }
// "synth_grid: grad must be 0, 1 or 2"
Grid bad = g;
bad.r.pop_back();                        // one radius fewer than latitudes
try { synth_grid(prep, bad, 0); }
catch (const std::invalid_argument &e) { std::printf("%s\n", e.what()); }
// "synth_grid: one radius per latitude is required"
```

### `GridTiming`
`lumps` and `fft` are thread-seconds (summed over threads). `total` is wall-clock time.

**Example**

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

const Model model = Model::zeros(3.986004415e14, 6378136.3, 4);
const auto prep = prepare(model);
const Grid g = Grid::gauss_legendre(4, 6378136.3);
GridTiming tm;                                        // members default to 0
const RawFields raw = synth_grid(prep, g, 1, SynthOptions{}, &tm);
std::printf("lumps %.6f  fft %.6f  total %.6f s\n", tm.lumps, tm.fft, tm.total);
// lumps and fft are summed over threads, total is wall clock time. With more than one thread,
// lumps + fft can be larger than total.
```

### `Points grid_to_points(grid)`
Expands the grid into scattered points. Used for validation; very slow at high degrees.

**Example**

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

const Grid reg = Grid::regular(45.0, 30.0, 7, 250.0, 260.0, 11, 6378136.3 + 5000.0);
const Points pts = grid_to_points(reg);
std::printf("%lld points\n", (long long)pts.size());                    // 77 = 7 * 11
const std::int64_t k = 2 * 11 + 4;                                       // row 2 (40 N), column 4 (254 E)
std::printf("lat %.3f  lon %.3f  r %.1f\n", rad2deg(pts.lat[k]), rad2deg(pts.lon[k]), pts.r[k]);
// lat 40.000  lon 254.000  r 6383136.3
```

## Internal functions (`grid.cpp`)

| Function | What it does |
|---|---|
| `auto_table(isa)` | Kernel table for the active ISA |
| `gauss_nodes(L, lat)` | Newton iteration on `P_L` starting from `cos(π(i+0.75)/(L+0.5))`, with exact symmetrisation |
| `uniform(lat, nlon, r)` | Creates a global grid with constant radius |
| `build_tasks(g)` | Sorts by \|φ\| and creates `LatTask{lat, r, row, row_mirror}` |
| `oracle_lumps(...)` | Lumps for one latitude with MSFC and `dphi`, without singularity |
| `emit` (lambda) | Converts lumps into one row per field, with the signs of the mirror row |

These functions are in an anonymous namespace in `grid.cpp` (`emit` is a local lambda), so user code cannot call them. The examples below show the behaviour each one produces, through the public API. `emit` is covered in "Symmetry ±φ".

**Example** (`auto_table`, through `active_isa_name`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// auto_table(isa) selects the kernel table; the public query for that choice is active_isa_name().
std::printf("auto:   %s\n", active_isa_name(Isa::Auto).c_str());     // "avx512", "avx2" or "scalar"
std::printf("forced: %s\n", active_isa_name(Isa::Scalar).c_str());   // "scalar"
```

**Example** (`gauss_nodes`, through `gauss_legendre`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// gauss_nodes(L, lat) gives the roots of P_L. For nmax = 2, L = 3 and the roots are 0 and +/- sqrt(3/5).
const Grid g = Grid::gauss_legendre(2, 6378136.3);
std::printf("%.4f %.4f %.4f deg\n", rad2deg(g.lat[0]), rad2deg(g.lat[1]), rad2deg(g.lat[2]));  // 50.7685 0.0000 -50.7685
std::printf("exact symmetry: %d\n", g.lat[0] == -g.lat[2] && g.lat[1] == 0.0);                  // 1
```

**Example** (`uniform`, through `gauss_legendre` and `driscoll_healy1`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// uniform(lat, nlon, r) sets one radius per row, dlon = 2 pi / nlon and lon0 = 0.
const Grid g = Grid::gauss_legendre(3, 6378136.3);   // L = 4, nlon = 8
std::printf("same radius: %d\n", g.r[0] == g.r[3]);                      // 1
std::printf("dlon = %.1f deg, lon0 = %.1f deg\n", rad2deg(g.dlon), rad2deg(g.lon0));  // 45.0, 0.0
```

**Example** (`build_tasks`, through `synth_grid`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

const Model model = Model::zeros(3.986004415e14, 6378136.3, 4);
const auto prep = prepare(model);
Grid g = Grid::gauss_legendre(4, 6378136.3);   // rows 0 and 4 are +/-64.98 deg, same radius: one mirror task
const RawFields a = synth_grid(prep, g, 0);
g.r[4] += 1000.0;                              // radii differ: the pair is no longer mirrored
const RawFields b = synth_grid(prep, g, 0);    // both rows are computed on their own
std::printf("row 0 difference: %g\n", a.f[0][0] - b.f[0][0]);   // ~0, up to rounding
```

**Example** (`oracle_lumps`, through `synth_grid` with `grad > 0`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// Rows with cos(lat) < 5e-3 use oracle_lumps() when grad > 0. driscoll_healy1 row 0 is the north pole.
const Model model = Model::zeros(3.986004415e14, 6378136.3, 4);   // C00 only
const auto prep = prepare(model);
const Grid dh = Grid::driscoll_healy1(4, 6378136.3);
const RawFields raw = synth_grid(prep, dh, 2);
// Compare entry (row 0, column 3) with the point oracle at the same place (longitude 3 * dlon).
const ReferencePoint ref = reference_point(model, legendre::make_coefs(4), PI / 2, 3 * dh.dlon, 6378136.3);
std::printf("field 0: grid %.6e, oracle %.6e\n", raw.f[0][3], ref.raw[0]);   // agree to about 1e-10 relative
```

**Example** (`field_map`, through the number of fields)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// field_map is internal. Its effect is visible in the field count: grad 0 -> 1, grad 1 -> 4, grad 2 -> 10.
const Model model = Model::zeros(3.986004415e14, 6378136.3, 4);
const auto prep = prepare(model);
const Grid g = Grid::gauss_legendre(4, 6378136.3);
for (int grad = 0; grad <= 2; ++grad)
  std::printf("grad %d: %zu fields\n", grad, synth_grid(prep, g, grad).f.size());   // 1, 4, 10
```

## Symmetry ±φ

`P_nm(−t) = (−1)^(n−m) P_nm(t)`. Hence the row at `−φ` reuses the even and odd lumps, changing signs. The groups with dP reverse the parity (`group_is_dp`). `field_map` gives, for each field, its group and factor (1, m or m²).

**Example** (parity signs in `emit`, `group_is_dp` and `field_map`)

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

// Only the degree-1 zonal term C10 is set. P_10 is odd in sin(phi).
Model model = Model::zeros(3.986004415e14, 6378136.3, 4);
model.C[idx_tri(1, 0, 4)] = 1e-5;
const auto prep = prepare(model);
const Grid g = Grid::gauss_legendre(4, 6378136.3);   // rows 0 and 4 are +/-64.98 deg, same radius
const RawFields raw = synth_grid(prep, g, 2);
const std::int64_t north = 0, south = 4 * g.nlon;
// Field 0 (group 0, plain P): the mirror row is minus the north row.
std::printf("f0 north %+.6e  south %+.6e\n", raw.f[0][north + 3], raw.f[0][south + 3]);
// Field 2 (group 2, dP/dphi, listed in group_is_dp): the mirror row equals the north row.
std::printf("f2 north %+.6e  south %+.6e\n", raw.f[2][north + 3], raw.f[2][south + 3]);
```

## Putting it together

A typical workflow: load a model, prepare it once, synthesise on a global Gauss-Legendre grid, then scale the result. The model path is a placeholder for a local ICGEM `.gfc` file.

```cpp
#include <cstdio>
#include "grid.hpp"
using namespace gsh;

int main()
{
  // 1. Load a fully normalized model (placeholder path; needs a local .gfc file).
  const Model model = load_gfc("EGM2008.gfc", 360);

  // 2. Prepare the coefficients once: degrees 2 to 360, reused for every grid below.
  const auto prep = prepare(model, 2, 360);

  // 3. Gauss-Legendre grid for nmax = 360: L = 361 latitudes x 2L = 722 longitudes, radius R.
  const Grid g = Grid::gauss_legendre(prep.nmax, model.R);

  // 4. Raw sums with grad 1 (4 fields), timed. threads = 8 overrides the default (0 = OpenMP maximum).
  GridTiming tm;
  SynthOptions opts;
  opts.threads = 8;
  const RawFields raw = synth_grid(prep, g, 1, opts, &tm);

  // 5. Scale field 0 to the potential in m^2/s^2: S0 * GM / R, as tests/py/test_gsh_tools.py does.
  const double potential_north_row_0 = raw.f[0][0] * model.GM / model.R;

  std::printf("rows %lld, cols %lld, fields %zu\n", (long long)g.nlat(), (long long)g.nlon, raw.f.size());
  std::printf("potential at row 0, column 0: %.6e m^2/s^2\n", potential_north_row_0);
  std::printf("time: %.3f s wall (lumps %.3f, fft %.3f thread-s)\n", tm.total, tm.lumps, tm.fft);
  return 0;
}
```
