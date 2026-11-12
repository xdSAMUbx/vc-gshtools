# functionals: geopotential functionals

Header `src/cpp/gsh_expansion/functionals.hpp`, implementation `src/cpp/gsh_expansion/functionals.cpp`. Formulas: internal study notes (not distributed).

**Examples setup.** The examples assume `#include "functionals.hpp"` and `using namespace gsh;`. That header includes `ellipsoid.hpp`, `height.hpp` and `shs.hpp`, so `Model`, `load_gfc`, `Points`, `Ellipsoid`, `HeightGrid`, `load_gravsoft_grid` and `height_from_grid` are visible too. Add `#include <cstdio>` for `std::printf` and `#include <stdexcept>` for `std::invalid_argument`. Link against the `gsh` target, which exports the `gsh_utils`, `gsh_legendre` and `gsh_expansion` include directories. `Points` stores geocentric latitude and longitude in radians; `Points::from_geodetic` takes degrees and converts to geocentric. Values are SI unless a comment says otherwise. The data paths are placeholders; point them at your copy of the files.

## `FunctionalResult compute_functionals(model, pts, req)`

Computes the functionals requested in `req.outputs` with the minimum number of synthesis passes.

### `FunctionalRequest`

| Field | Default | Meaning |
|---|---|---|
| `outputs` | | List of `Functional` |
| `nmin`, `nmax` | 0, -1 | Degrees. `nmax = -1` uses the model's degree |
| `ellipsoid` | GRS80 | Normal field |
| `height` | null | `HeightSource`. **Required** for N and ζ |
| `density` | 2670 | Density ρ in kg/m³, for N and ζ |
| `G` | 6.67259e-11 | Gravitational constant, for N and ζ |
| `pole_cos_threshold` | 5e-3 | Below this `cos φ` the point is evaluated with the oracle |
| `synth` | `{}` | [`SynthOptions`](shs.md#struct-synthoptions) |

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);                       // ICGEM file, degrees up to 360
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});    // about 1 km above the equator

FunctionalRequest rq;
rq.outputs = {Functional::g, Functional::g_sa};  // only these are evaluated
rq.nmin = 2;                                     // lowest degree (default 0)
rq.nmax = 360;                                   // highest degree (default -1 = model.nmax)
rq.ellipsoid = Ellipsoid::GRS80();               // default, written out here
rq.pole_cos_threshold = 5e-3;                    // default: cos(lat) below this uses the oracle
rq.synth.threads = 4;                            // OpenMP threads (default 0 = all)

const auto res = compute_functionals(model, pts, rq);
std::printf("g = %.6f m/s^2\n", res.values.at(Functional::g)[0]);  // about 9.8 m/s^2
```

### `FunctionalResult`

| Field | Meaning |
|---|---|
| `values` | Map `Functional` → values, one per point |
| `points_without_height` | Points with NaN in N and ζ |
| `polar_points` | Points resolved with the oracle |
| `timing_s` | Seconds per stage: `model_synthesis`, `disturbing_synthesis`, `height`, `ellipsoid_synthesis` |

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0, 89.8}, {-70.0, 0.0}, {6378137.0, 6378137.0});
FunctionalRequest rq;
rq.outputs = {Functional::g};
const FunctionalResult res = compute_functionals(model, pts, rq);

const AlignedVec<double> &g = res.values.at(Functional::g);  // one entry per input point
std::printf("n = %zu, g[0] = %.6f m/s^2\n", g.size(), g[0]);
std::printf("polar_points = %lld\n", static_cast<long long>(res.polar_points));                    // 1 (lat 89.8 deg)
std::printf("points_without_height = %lld\n", static_cast<long long>(res.points_without_height));  // 0, no N or zeta
for (const auto &kv : res.timing_s)  // seconds per stage; only "model_synthesis" runs for g
  std::printf("%-22s %.4f s\n", kv.first.c_str(), kv.second);
```

### Errors

- `invalid_argument` if `nmax` exceeds the model's degree.
- `invalid_argument` if a functional containing `1/cos φ` is requested at an exact pole: `g_y`, `g`, `dg`, `V_rlambda`, `V_philambda`, `V_lambdalambda`, `T_rlambda`, `T_philambda`, `T_lambdalambda`, `eta`, `Theta`.
- `invalid_argument` if N or ζ is requested without `height`.

### Example

```cpp
const Ellipsoid ell = Ellipsoid::GRS80();
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc");
const Points pts = Points::from_geodetic({38.2}, {-105.6}, {5376.0}, ell.a, ell.e2);

FunctionalRequest rq;
rq.outputs = {Functional::T, Functional::Dg_sa, Functional::N, Functional::zeta};
rq.height = height_from_grid(load_gravsoft_grid("data/cache/tablas/Colorado_dem/Colorado_dem.grd"), 5.0 / 60.0);

const auto res = compute_functionals(model, pts, rq);
const double zeta_m = res.values.at(Functional::zeta)[0];
const double Dg_mGal = res.values.at(Functional::Dg_sa)[0] * 1e5;  // the API returns m/s²
```

**Example (errors)**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6378137.0});
FunctionalRequest rq;
rq.outputs = {Functional::V};

rq.nmax = 361;  // larger than model.nmax = 360
try { compute_functionals(model, pts, rq); }
catch (const std::invalid_argument &e) { std::printf("%s\n", e.what()); }
// compute_functionals: nmax is larger than the model degree

rq.nmax = -1;
rq.outputs = {Functional::N};  // N needs rq.height, which is not set
try { compute_functionals(model, pts, rq); }
catch (const std::invalid_argument &e) { std::printf("%s\n", e.what()); }
// N and zeta require a height source (HeightSource)
```

```cpp
// Exact north pole (geocentric latitude pi/2): eta contains 1/cos(lat) and is rejected.
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
Points pole;
pole.lat.push_back(PI / 2);
pole.lon.push_back(0.0);
pole.r.push_back(6356752.314);
FunctionalRequest rq;
rq.outputs = {Functional::eta};
try { compute_functionals(model, pole, rq); }
catch (const std::invalid_argument &e) { std::printf("%s\n", e.what()); }
// The functional eta is not defined at an exact pole
```

## Available functionals

`S0`, `SR`, etc. are [raw sums](../README.md#sumas-crudas). "Model" uses C̄, S̄; "disturbing" uses ΔC̄, ΔS̄. γ is the normal gravity at `(r, φ)`. The `Functional` enum uses these names.

### Potential and gravity (model)

| Name | Formula | Unit |
|---|---|---|
| `V` | `GM/r · S0` | m²/s² |
| `dV_r` | `−GM/r² · SR` | m/s² |
| `dV_phi` | `GM/r · SP` | m²/s² |
| `dV_lambda` | `GM/r · SL` | m²/s² |
| `W` | `V + ½ω²r²cos²φ` | m²/s² |
| `g_x` | `(dV_phi − ω²r²cosφ sinφ)/r` | m/s² |
| `g_y` | `−dV_lambda/(r cosφ)` | m/s² |
| `g_z` | `dV_r + ω²r cos²φ` | m/s² |
| `g` | `√(g_x² + g_y² + g_z²)` | m/s² |
| `g_sa` | `\|g_z\|` | m/s² |
| `W_rr` | `V_rr + ω²cos²φ` | s⁻² |

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::V, Functional::dV_r, Functional::dV_phi, Functional::dV_lambda,
              Functional::W, Functional::g_x, Functional::g_y, Functional::g_z,
              Functional::g, Functional::g_sa, Functional::W_rr};
const auto res = compute_functionals(model, pts, rq);
for (Functional f : rq.outputs)
  std::printf("%-10s %.9g\n", functional_name(f), res.values.at(f)[0]);
// V, dV_phi, dV_lambda, W: m^2/s^2. dV_r, g_x, g_y, g_z, g, g_sa: m/s^2 (g about 9.8).
// W_rr: s^-2; the centrifugal part omega^2 cos^2(45 deg) = 2.66e-9 s^-2 is included.
```

### Tensor (model)

With `k = GM/r³`. The `T_*` equivalents of the disturbing potential use the same formulas on ΔC̄, ΔS̄.

| Name | Formula |
|---|---|
| `V_rr` | `k · SRR` |
| `V_rphi` | `−k · SRP` |
| `V_rlambda` | `−k/cosφ · SRL` |
| `V_phiphi` | `k · SPP` |
| `V_philambda` | `k/cosφ · SPL` |
| `V_lambdalambda` | `−k/cos²φ · SLL` |

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::V_rr, Functional::V_rphi, Functional::V_rlambda,
              Functional::V_phiphi, Functional::V_philambda, Functional::V_lambdalambda};
const auto res = compute_functionals(model, pts, rq);
for (Functional f : rq.outputs)  // second derivatives in s^-2 (1 E = 1e-9 s^-2)
  std::printf("%-16s %+.4e s^-2\n", functional_name(f), res.values.at(f)[0]);
// V_rr is of order 1e-6 s^-2 (thousands of Eotvos).
```

Tensor in the local frame, with `tn = tanφ`. Units are s⁻²:

| Name | Formula |
|---|---|
| `V_xx` | `k (SPP − SR)` |
| `V_xy` | `−k (SPL + tn SL)/cosφ` |
| `V_xz` | `−k (SRP + SP)` |
| `V_yy` | `k (−SLL/cos²φ − SR − tn SP)` |
| `V_yz` | `k (SRL + SL)/cosφ` |
| `V_zz` | `k · SRR` |

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::V_xx, Functional::V_xy, Functional::V_xz,
              Functional::V_yy, Functional::V_yz, Functional::V_zz, Functional::V_rr};
const auto res = compute_functionals(model, pts, rq);
const auto v = [&](Functional f) { return res.values.at(f)[0]; };
std::printf("V_zz - V_rr = %.3e s^-2\n", v(Functional::V_zz) - v(Functional::V_rr));  // V_zz is the same raw sum: ~0
std::printf("trace       = %.3e s^-2\n", v(Functional::V_xx) + v(Functional::V_yy) + v(Functional::V_zz));  // Laplace: ~0
```

### Disturbing potential and anomalies

| Name | Formula | Unit | What it is |
|---|---|---|---|
| `T` | `GM/r · S0` (disturbing) | m²/s² | Disturbing potential |
| `dg_sa` | `GM/r² · SR` | m/s² | Disturbance, spherical approximation: `−∂T/∂r` |
| `Dg_sa` | `GM/r² · (SR − 2 S0)` | m/s² | Gravity anomaly, spherical approximation |
| `dg` | `g − γ` | m/s² | Exact disturbance. Uses the full model |

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::T, Functional::dg_sa, Functional::Dg_sa, Functional::dg};
const auto res = compute_functionals(model, pts, rq);
std::printf("T      = %.4f m^2/s^2\n", res.values.at(Functional::T)[0]);  // of order 1e3
std::printf("dg_sa  = %.4f mGal\n", res.values.at(Functional::dg_sa)[0] * 1e5);  // m/s^2 -> mGal
std::printf("Dg_sa  = %.4f mGal\n", res.values.at(Functional::Dg_sa)[0] * 1e5);
std::printf("dg     = %.4f mGal\n", res.values.at(Functional::dg)[0] * 1e5);
// Dg_sa = dg_sa - 2T/r, so the two differ by about 2T/r (roughly 30 mGal for T of 1e3 m^2/s^2 at r = 6.38e6 m).
```

### Disturbing tensor

`T_rr`, `T_rphi`, `T_rlambda`, `T_phiphi`, `T_philambda`, `T_lambdalambda`, `T_xx`, `T_xy`, `T_xz`, `T_yy`, `T_yz`, `T_zz`.

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::T_rr, Functional::T_rphi, Functional::T_rlambda, Functional::T_phiphi,
              Functional::T_philambda, Functional::T_lambdalambda, Functional::T_xx, Functional::T_xy,
              Functional::T_xz, Functional::T_yy, Functional::T_yz, Functional::T_zz};
const auto res = compute_functionals(model, pts, rq);
for (Functional f : rq.outputs)  // s^-2; multiply by 1e9 for Eotvos
  std::printf("%-16s %+.4e s^-2 (%+.3f E)\n", functional_name(f), res.values.at(f)[0], res.values.at(f)[0] * 1e9);
// T_zz is the same raw sum as T_rr, so the two agree.
```

### Deflections and heights

| Name | Formula | Unit | Where it is evaluated |
|---|---|---|---|
| `xi` | `−GM/(r²γ) · SP` | rad | At the point |
| `eta` | `−GM/(r²γ cosφ) · SL` | rad | At the point |
| `Theta` | `√(ξ² + η²)` | rad | At the point |
| `zeta_ell` | `T/γ` | m | At the point: **at flight height** if the point is in flight |
| `N` | `(T − 2πGρH²)/γ` | m | On the ellipsoid |
| `zeta` | `T/γ − δg_sa (H + N)/γ` | m | On the ellipsoid |

For N and ζ, T, δg_sa and γ are evaluated at `r_ell(φ)`, at the same spherical latitude as the point and at its longitude. H is looked up at the geodetic latitude of that foot point.

**To get T/γ on the ellipsoid without a DEM,** create `Points` with `r = ell.radius(lat)` and request `zeta_ell`. This is N without the topography term. `zeta_on_ellipsoid` in `examples/cpp/ex_aerogravimetria.cpp` does this.

**Example (deflections)**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::xi, Functional::eta, Functional::Theta};
const auto res = compute_functionals(model, pts, rq);
const double xi = res.values.at(Functional::xi)[0], eta = res.values.at(Functional::eta)[0];
std::printf("xi    = %+.3f arcsec\n", rad2deg(xi) * 3600.0);   // values are stored in radians
std::printf("eta   = %+.3f arcsec\n", rad2deg(eta) * 3600.0);
std::printf("Theta = %.3f arcsec\n", rad2deg(res.values.at(Functional::Theta)[0]) * 3600.0);  // sqrt(xi^2 + eta^2)
```

**Example (N and zeta, with a DEM)**

```cpp
const Ellipsoid ell = Ellipsoid::GRS80();
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_geodetic({38.2}, {-105.6}, {5376.0}, ell.a, ell.e2);
FunctionalRequest rq;
rq.outputs = {Functional::N, Functional::zeta};
rq.height = height_from_grid(load_gravsoft_grid("data/cache/tablas/Colorado_dem/Colorado_dem.grd"), 5.0 / 60.0);
const auto res = compute_functionals(model, pts, rq);
std::printf("N    = %.3f m\n", res.values.at(Functional::N)[0]);
std::printf("zeta = %.3f m\n", res.values.at(Functional::zeta)[0]);
std::printf("without height: %lld\n", static_cast<long long>(res.points_without_height));  // 0 if the DEM covers the point
std::printf("r_ell(0) = %.1f m, r_ell(pi/2) = %.1f m\n", ell.radius(0.0), ell.radius(PI / 2));  // 6378137.0, about 6356752.3
```

**Example (zeta_ell on the ellipsoid, no DEM)**

```cpp
const Ellipsoid ell = Ellipsoid::GRS80();
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_geodetic({38.2}, {-105.6}, {5376.0}, ell.a, ell.e2);
Points foot = pts;  // same geocentric latitudes, radius moved to the ellipsoid
for (std::int64_t i = 0; i < foot.size(); ++i)
  foot.r[i] = ell.radius(foot.lat[i]);
FunctionalRequest rq2;
rq2.outputs = {Functional::zeta_ell};
const double T_over_gamma = compute_functionals(model, foot, rq2).values.at(Functional::zeta_ell)[0];  // metres
```

## What is computed and when

`need_of(f)` indicates, for each functional, the gradient order required on the model and on ΔC̄, and whether it needs γ or a pass on the ellipsoid. `compute_functionals` takes the maximum over all requested functionals, so requesting `g` and `V` costs a single G1 pass. N and ζ always add a G1 pass on the ellipsoid.

**Near-polar points** (`cos φ < pole_cos_threshold`): the sums with derivatives and the local tensor are replaced by `reference_point`, one call per point. This is slow if there are many such points.

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({89.8, 89.9999}, {0.0, 0.0}, {6378137.0, 6356752.314});
FunctionalRequest rq;
rq.outputs = {Functional::V_xx, Functional::V_zz, Functional::g_x};
const auto res = compute_functionals(model, pts, rq);
// cos(89.8 deg) = 3.5e-3 and cos(89.9999 deg) = 1.7e-6, both below 5e-3: both go through reference_point().
std::printf("polar points: %lld\n", static_cast<long long>(res.polar_points));  // 2
```

**Actual cost.** `dg_sa`, `Dg_sa`, `T_rr` and `V_rr` use G1 or G2 kernels even though a G0 with weights would be enough. The plan foresees that kernel and it is not implemented.

## Other functions

| Function | What it does |
|---|---|
| `functional_name(f)` | Name of the functional as a string |
| `all_functionals()` | List of all 45 |
| `local_tensor(S, k3, lat)` (internal) | Local tensor from the raw sums |
| `need_of(f)` (internal) | Passes required by a functional |

**`functional_name(f)` and `all_functionals()`**

**Example**

```cpp
std::printf("%s\n", functional_name(Functional::Dg_sa));  // "Dg_sa"
std::printf("%s\n", functional_name(Functional::zeta));   // "zeta"
const std::vector<Functional> all = all_functionals();
std::printf("count = %zu\n", all.size());                  // 45
for (Functional f : all)                                   // enum order, starting with "V"
  std::printf("%s ", functional_name(f));
std::printf("\n");
```

**`need_of(f)` (internal)**

`need_of` is in an anonymous namespace in `functionals.cpp` and cannot be called from user code. Its effect is visible through `timing_s`, which lists only the passes that ran.

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::T, Functional::dg_sa};  // disturbing model only, gradient order 1
FunctionalResult res = compute_functionals(model, pts, rq);
std::printf("model_synthesis: %d\n", static_cast<int>(res.timing_s.count("model_synthesis")));          // 0
std::printf("disturbing_synthesis: %d\n", static_cast<int>(res.timing_s.count("disturbing_synthesis")));  // 1
rq.outputs = {Functional::V, Functional::g};      // one gradient-1 model pass serves both
res = compute_functionals(model, pts, rq);
std::printf("model_synthesis: %d\n", static_cast<int>(res.timing_s.count("model_synthesis")));          // 1
```

**`local_tensor(S, k3, lat)` (internal)**

`local_tensor` is defined in `functionals.cpp` (anonymous namespace) and is not declared in the header. The public path to it is the `V_xx`...`V_zz` (and `T_*`) outputs, as below.

**Example**

```cpp
const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
const Points pts = Points::from_spherical({45.0}, {-70.0}, {6379137.0});
FunctionalRequest rq;
rq.outputs = {Functional::V_xx, Functional::V_yy, Functional::V_zz};
const auto res = compute_functionals(model, pts, rq);
const double trace = res.values.at(Functional::V_xx)[0] + res.values.at(Functional::V_yy)[0] + res.values.at(Functional::V_zz)[0];
std::printf("trace = %.3e s^-2\n", trace);  // Laplace equation: roundoff-sized
```

## Putting it together

A typical gravity workflow: load the model and the observation points, request the functionals once, and convert the outputs to gravimetric units.

```cpp
#include <cstdio>
#include "functionals.hpp"
using namespace gsh;

int main()
{
  const Ellipsoid ell = Ellipsoid::GRS80();
  const Model model = load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);  // ICGEM .gfc, fully normalized
  // Two observation points: geodetic degrees and ellipsoidal height in metres.
  const Points pts = Points::from_geodetic({38.2, 39.0}, {-105.6, -105.0}, {5376.0, 2200.0}, ell.a, ell.e2);

  FunctionalRequest rq;
  rq.outputs = {Functional::g, Functional::dg, Functional::Dg_sa, Functional::xi, Functional::eta,
                Functional::N, Functional::zeta};
  rq.height = height_from_grid(load_gravsoft_grid("data/cache/tablas/Colorado_dem/Colorado_dem.grd"), 5.0 / 60.0);  // 5 arcmin blocks
  const FunctionalResult res = compute_functionals(model, pts, rq);

  for (std::int64_t i = 0; i < pts.size(); ++i)
  {
    std::printf("point %lld\n", static_cast<long long>(i));
    std::printf("  g     = %.4f mGal\n", res.values.at(Functional::g)[i] * 1e5);   // m/s^2 -> mGal
    std::printf("  dg    = %.4f mGal\n", res.values.at(Functional::dg)[i] * 1e5);
    std::printf("  Dg_sa = %.4f mGal\n", res.values.at(Functional::Dg_sa)[i] * 1e5);
    std::printf("  xi    = %+.3f arcsec\n", rad2deg(res.values.at(Functional::xi)[i]) * 3600.0);
    std::printf("  eta   = %+.3f arcsec\n", rad2deg(res.values.at(Functional::eta)[i]) * 3600.0);
    std::printf("  N     = %.3f m\n", res.values.at(Functional::N)[i]);
    std::printf("  zeta  = %.3f m\n", res.values.at(Functional::zeta)[i]);
  }
  std::printf("points without height: %lld\n", static_cast<long long>(res.points_without_height));
}
```
