# shs: synthesis at scattered points

Header `src/cpp/gsh_expansion/shs.hpp`. Implementation in `src/cpp/gsh_expansion/shs.cpp` (scheduler) and `src/cpp/gsh_expansion/reference.cpp` (oracle). The CPU kernels are described in [kernels](kernels.md).

**Examples.** Each example is self-contained: it lists the headers it needs and uses `using namespace gsh;`. The headers are in `src/cpp/gsh_utils`, `gsh_legendre` and `gsh_expansion`, which the `gsh` target exports as public include directories, so `#include "shs.hpp"` works when you link `gsh`. `kernels.hpp` is private to the library, so the field order of `RawFields::f` is written out in the comments instead. Model paths assume that the working directory is the repository root. The numeric values in the comments are approximate unless stated otherwise.

## Points

### `struct Points`
`lat`, `lon` (radians) and `r` (metres), in parallel arrays. `size()` returns the number of points.

**Example**

```cpp
#include "synth_types.hpp"
using namespace gsh;

Points p;
p.lat = {0.785398163397448};    // 45 deg, in radians
p.lon = {-1.2915436466866565};  // -74 deg, in radians
p.r   = {6378237.0};            // metres (about 100 m above the equatorial radius)
// p.size() == 1
```

### `static Points from_spherical(lat_deg, lon_deg, r)`
Spherical latitude and longitude in degrees, radius in metres. Throws `invalid_argument` if the sizes differ.

**Example**

```cpp
#include <stdexcept>
#include <vector>
#include "synth_types.hpp"
using namespace gsh;

std::vector<double> lat_deg = {45.0, -12.5};
std::vector<double> lon_deg = {-74.0, 77.0};
std::vector<double> r_m     = {6378237.0, 6381137.0};  // about 100 m and 3000 m above the equatorial radius
Points pts = Points::from_spherical(lat_deg, lon_deg, r_m);
// pts.size() == 2; pts.lat[0] == 0.7853981... (45 deg in radians); pts.r[1] == 6381137.0

try {
  Points::from_spherical({0.0, 1.0}, {0.0}, {6.4e6, 6.4e6});  // sizes 2, 1, 2 differ
} catch (const std::invalid_argument &) {
  // thrown as documented
}
```

### `static Points from_geodetic(lat_deg, lon_deg, h, a, e2)`
Geodetic coordinates in degrees and ellipsoidal height in metres, converted with `geodetic_to_spherical`. Throws `invalid_argument` if the sizes differ.

**Example**

```cpp
#include "ellipsoid.hpp"
#include "synth_types.hpp"
using namespace gsh;

Ellipsoid ell = Ellipsoid::GRS80();  // a = 6378137 m, e2 = 0.00669438002290
Points pts = Points::from_geodetic({40.0}, {-3.7}, {650.0}, ell.a, ell.e2);
// pts.lat[0] is the GEOCENTRIC latitude, about 39.81 deg (not 40.0), in radians
// pts.r[0] is the geocentric radius, about 6.37e6 m; the height h = 650 m is folded into it
```

## Options

### `struct SynthOptions`

| Field | Default | Effect |
|---|---|---|
| `threads` | 0 | Threads. 0 uses the OpenMP maximum |
| `tile_vectors` | 0 | SIMD vectors per tile: 1, 2, 4, 8 or 16. 0 selects 4 (G0), 2 (G1) or 1 (G2), according to measurements |
| `isa` | `Auto` | `Auto`, `Scalar`, `Avx2`, `Avx512`. Forcing an ISA that is not available throws `runtime_error` |
| `polar_optimization` | true | Skips m columns that do not contribute at high latitudes |
| `skip_forbidden` | true | Does not accumulate where \|P\| < 2^-200 in all lanes |
| `rotate_trig` | true | `cos/sin(mλ)` by rotation. Applied only in G0 |
| `sort_latitude` | true | Groups points of similar \|φ\| in each tile |
| `parallel_orders` | -1 | -1 automatic, 0 per tile, 1 per order m |

**Example**

```cpp
#include "shs.hpp"
using namespace gsh;

SynthOptions opts;            // defaults: all OpenMP threads, automatic tiles, Isa::Auto
opts.threads = 8;             // 0 = OpenMP maximum
opts.tile_vectors = 4;        // 1, 2, 4, 8 or 16; other values throw invalid_argument when synth_points runs
opts.isa = Isa::Avx2;         // a CPU without AVX2/FMA makes synth_points throw runtime_error
opts.sort_latitude = true;    // group points of similar |latitude| inside each tile
```

## Preparation

### `struct Recurrence` and `make_recurrence(int nmax)`
Recurrence coefficients for one degree: `ab` (pairs a, b per (n, m)), `f` (for dP/dφ) and `d` (sectorials). Shared with `shared_ptr` between models with the same nmax. `compute_functionals` reuses it between the model and ΔC̄.

**Example**

```cpp
#include <memory>
#include "shs.hpp"
using namespace gsh;

std::shared_ptr<const Recurrence> rec = make_recurrence(360);
// rec->nmax == 360; rec->d has 361 entries (d[0] = 1, d[1] = sqrt(3))
// rec->ab holds an (a, b) pair per packed (n, m) index; rec->f holds the dP/dphi factors
// Build it once and pass it to every prepare() call that uses nmax = 360.
```

### `Prepared prepare(model, nmin = 0, nmax = -1, rec = nullptr)`
Prepares the model for the kernels.

- Packs `[C, S]` per (n, m) into a single array.
- Sets to 0 the degrees `n < nmin`.
- `nmax < 0` uses the model's degree.
- With `rec = nullptr` it creates the recurrence.

**Errors:** `invalid_argument` if `nmax` exceeds the model's degree, `nmin` is out of range, or `rec` is for another degree.

**Example**

```cpp
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);  // keeps degrees <= 360
Prepared full = prepare(model, 0, 360);            // builds its own recurrence tables
Prepared tail = prepare(model, 2, 360, full.rec);  // shares them; degrees 0 and 1 are zeroed
// full.nmin == 0, full.nmax == 360, full.GM == model.GM; full.cs holds [C, S] per (n, m)
// prepare(model, 0, 400) throws invalid_argument (400 > 360); prepare(model, -1) throws (nmin < 0)
```

## Synthesis

### `struct RawFields`
`grad` (0, 1 or 2) and `f`, one array per field, each with one value per input point. Field order is that of the internal `kernels::Field` enum (a private header, not reachable from user code; use the indices): 0 `S0`, 1 `SR`, 2 `SP`, 3 `SL`, 4 `SRR`, 5 `SRP`, 6 `SRL`, 7 `SPP`, 8 `SPL`, 9 `SLL`. Grad 0 fills only index 0, grad 1 fills 0 to 3.

**Example**

```cpp
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
Prepared prep = prepare(model, 0, 360);
Points pts = Points::from_spherical({45.0}, {-74.0}, {6378237.0});
RawFields rf = synth_points(prep, pts, 2);
// rf.grad == 2; rf.f.size() == 10; rf.f[k][i] is field k at point i
// rf.f[0][0] (S0) is close to 1 near the surface, because the n = 0 term of the potential is 1
```

### `RawFields synth_points(prep, pts, grad, opts = {})`

Computes the [raw sums](../README.md#sumas-crudas) at each point, without `GM/r^k`.

- **`grad`:** 0, 1 or 2, which give 1, 4 or 10 fields in `RawFields::f`. Any other value throws `invalid_argument`.
- **Result:** each field has `pts.size()` values, in input order.

**CPU steps:**

1. Chooses the kernel table according to the ISA.
2. Sorts the points by |φ| if `sort_latitude` is active.
3. Forms tiles of `vectors × WIDTH` points. The last tile is padded by repeating its last point.
4. Parallelises:
   - **Per tile (`schedule(dynamic,1)`):** with many points. The result does not depend on the number of threads.
   - **Per order m:** if there are fewer than `2·threads` tiles. This is the CHarm method; the summation order varies, so results are not bit-for-bit identical between runs. It calls `omp_set_num_threads`, which changes the global OpenMP state.
5. Writes each field to the original position of the point.

**Poles:** the result with grad ≥ 1 is not reliable if `cos φ` is very small, because dP is singular per column. `compute_functionals` corrects it with the oracle; if you call `synth_points` directly, you must handle this yourself.

**Example**

```cpp
#include <cstdint>
#include <cstdio>
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
Prepared prep = prepare(model, 0, 360);
Points pts = Points::from_spherical({45.0, -12.5}, {-74.0, 77.0}, {6378237.0, 6381137.0});

RawFields rf = synth_points(prep, pts, 0);   // grad 0: rf.f holds one array
SynthOptions opts;
opts.threads = 4;                            // optional tuning
RawFields rf2 = synth_points(prep, pts, 0, opts);

for (std::int64_t i = 0; i < pts.size(); ++i) {
  double V = prep.GM / pts.r[i] * rf.f[0][i];  // V in m^2 s^-2, about 6.2e7 near the surface
  std::printf("point %lld: S0 = %.9f, V = %.3f m2/s2\n", static_cast<long long>(i), rf.f[0][i], V);
}

RawFields g1 = synth_points(prep, pts, 1);   // rf.f has 4 arrays: S0, SR, SP, SL
double dV_dphi = prep.GM / pts.r[0] * g1.f[2][0];  // dV/dphi in m^2 s^-2 per radian
// With grad >= 1, keep |lat| away from 90 deg here; use reference_point() or compute_functionals() near the poles.
```

### `std::string active_isa_name(Isa isa = Isa::Auto)`
Name of the ISA that would be used: `"avx512"`, `"avx2"` or `"scalar"`.

**Example**

```cpp
#include <iostream>
#include <stdexcept>
#include "shs.hpp"
using namespace gsh;

std::cout << active_isa_name() << "\n";             // "avx512", "avx2" or "scalar": the best one this CPU has
std::cout << active_isa_name(Isa::Scalar) << "\n";  // always "scalar"
try {
  std::cout << active_isa_name(Isa::Avx512) << "\n";  // prints only if the CPU supports AVX-512
} catch (const std::runtime_error &e) {
  std::cout << "not available: " << e.what() << "\n";
}
```

## Oracle

### `ReferencePoint reference_point(model, coefs, lat, lon, r, nmin = 0, nmax = -1)`

Evaluates one point sequentially with the full triangle, using MSFC and `dphi`, and the literal formulas of the specification §2 and §4.2 (internal study notes, not distributed). It is slow. Used as a reference in tests and for near-polar points.

- **`coefs`:** from `legendre::make_coefs`, with degree ≥ `nmax`.
- **Returns:** `raw[10]` (raw sums) and `vxx … vzz` (local tensor in s⁻², already multiplied by `GM/r³`).
- **Poles:** no singularity.

**Example**

```cpp
#include <cstdio>
#include "legendre.hpp"
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
legendre::Coefs coefs = legendre::make_coefs(360);  // coefs.nmax >= nmax

// 45 N, 74 W, 100 m above the equatorial radius: radians and metres
double lat = deg2rad(45.0), lon = deg2rad(-74.0), r = 6378237.0;
ReferencePoint ref = reference_point(model, coefs, lat, lon, r, 0, 360);
double V = model.GM / r * ref.raw[0];  // V in m^2 s^-2; raw[0] is dimensionless
std::printf("V = %.3f m2/s2, vzz = %.6e s^-2\n", V, ref.vzz);

// At the exact pole the oracle has no singularity, so this call is valid
ReferencePoint np = reference_point(model, coefs, PI / 2, 0.0, 6356752.0, 0, 360);
// Oracle cost is O(nmax^2) per point: fine for a few points, too slow for a grid
```

### `struct ReferencePoint`
`raw[10]` holds the raw sums in `kernels::Field` order (the same order as `RawFields::f`). `vxx`, `vxy`, `vxz`, `vyy`, `vyz`, `vzz` is the LNOF tensor in s⁻², already multiplied by GM/r³.

**Example**

```cpp
#include <cstdio>
#include "legendre.hpp"
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
legendre::Coefs coefs = legendre::make_coefs(360);
ReferencePoint ref = reference_point(model, coefs, deg2rad(45.0), deg2rad(-74.0), 6378237.0, 0, 360);

// raw[0] = S0, raw[1] = SR (same order as RawFields::f at the same point)
std::printf("raw[0] = %.12f  raw[1] = %.12f\n", ref.raw[0], ref.raw[1]);
std::printf("vxx = %.6e  vzz = %.6e s^-2\n", ref.vxx, ref.vzz);  // already scaled by GM/r^3
```

## Internal functions

| Function | File | What it does |
|---|---|---|
| `cpu_has_avx512()`, `cpu_has_avx2()` | `shs.cpp` | Detect the ISA with CPUID and check that the operating system saves the registers (XGETBV) |
| `select(isa)` | `shs.cpp` | Returns the kernel table. `Auto` picks the best available one |
| `b_index(vectors)` | `shs.cpp` | Converts 1, 2, 4, 8, 16 into a table index |
| `sq(v)` | `reference.cpp` | `√v`, with negative radicand → 0 (§4.2) |
| `kron(a, b)` | `reference.cpp` | Kronecker delta |

These helpers live in anonymous namespaces inside `shs.cpp` and `reference.cpp`, so user code cannot call them and no example is given. Use `active_isa_name()` to see which ISA is selected.

## Putting it together

The example below loads a model, prepares it, synthesises the raw sums at three points, converts them to V and compares with the oracle. Expected result: V of about 6.2e7 m²/s² at each point, and differences at round-off level relative to V.

```cpp
#include <cmath>
#include <cstdint>
#include <cstdio>
#include "legendre.hpp"
#include "model.hpp"
#include "shs.hpp"
using namespace gsh;

int main()
{
  const int nmax = 360;
  Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", nmax);  // 1. load model
  Prepared prep = prepare(model, 0, nmax);                                  // 2. prepare
  legendre::Coefs coefs = legendre::make_coefs(nmax);                       // needed by the oracle

  // 3. Points: geocentric latitude and longitude in degrees, radius in metres
  Points pts = Points::from_spherical({45.0, -12.5, 0.0}, {-74.0, 77.0, 0.0},
                                      {6378237.0, 6381137.0, 6378137.0});
  RawFields rf = synth_points(prep, pts, 0);                                // 4. raw sums (dimensionless)

  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    double V = prep.GM / pts.r[i] * rf.f[0][i];                             // 5. V = GM/r * S0 (m^2 s^-2)
    ReferencePoint ref = reference_point(model, coefs, pts.lat[i], pts.lon[i], pts.r[i], 0, nmax);  // 6. oracle
    double Vref = model.GM / pts.r[i] * ref.raw[0];
    std::printf("point %lld: V = %.6f  Vref = %.6f  |diff| = %.3e m2/s2\n",
                static_cast<long long>(i), V, Vref, std::fabs(V - Vref));
  }
  return 0;
}
```
