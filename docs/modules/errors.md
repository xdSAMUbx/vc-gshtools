# errors: commission and omission errors

Header `src/cpp/gsh_expansion/errors.hpp`, implementation in `errors.cpp`. Omission model: Tscherning and Rapp (1974); functional weights: Bucha and Janák (2019). Include `errors.hpp` and link against `gsh`.

All functionals are linear in the coefficients, so the errors use the spherical approximation of the references (no ellipsoidal corrections).

| `ErrorFunctional` | Degree weight w_n | Unit |
|---|---|---|
| `T` | GM/r (R/r)^n | m²/s² |
| `ZetaEll` (T/γ) | GM/(r γ) (R/r)^n | m |
| `GravityAnomaly` (spherical Δg, as `Dg_sa`) | GM/r² (n−1) (R/r)^n | m/s² |
| `GravityDisturbance` (radial) | GM/r² (n+1) (R/r)^n | m/s² |

## Commission error

The `.gfc` format has only σC and σS, so the covariance of the coefficients is taken as **diagonal** (correlations are ignored). `load_gfc_sigma(path, nmax)` returns a `Model` whose `C` and `S` hold σC and σS.

`σ_F² = Σ_n w_n² Σ_m P̄_nm² (σC² cos² mλ + σS² sin² mλ)`

- `commission_error(sigma, pts, f, opts)`: error at every point. One Legendre triangle per point (`O(nmax²)`), parallel over points. `ErrorOptions` has `nmin`, `nmax`, `ellipsoid` (for γ) and `threads`.
- `commission_per_degree(sigma, f, r, gamma, nmin, nmax)`: area-mean contribution of each degree.
- `commission_global(sigma, f, r, gamma, nmin, nmax)`: quadrature sum of the previous one. For `T`: `CE(T) = (GM/a) sqrt(Σ (a/r)^(2n+2) Σ_m (σC²+σS²))`, and `CE(N) = CE(T)/γ ≈ R sqrt(Σ Σ (σC²+σS²))`.

The result is an absolute error only if the model's sigmas are calibrated (EGM2008 is); otherwise it measures relative geographic precision.

```cpp
#include "errors.hpp"
using namespace gsh;

Model sigma = load_gfc_sigma("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
Points pts = Points::from_geodetic({38.2}, {-105.6}, {5376.0}, 6378137.0, 0.00669437999014);
std::vector<double> ce = commission_error(sigma, pts, ErrorFunctional::ZetaEll);  // metres, one per point
double global = commission_global(sigma, ErrorFunctional::ZetaEll, 6378136.3, 9.81, 2, 360);  // ~0.078 m
```

## Omission error

Global (area-mean) signal of the degrees above `nmax`, from degree variances `σ_n² = Σ_m (C²+S²)`.

- `degree_variances(model)`: from the coefficients of a model.
- `tscherning_rapp_variances(GM, R, n_first, n_last, p)`: from the Tscherning–Rapp model `c_n = A (n−1)/((n−2)(n+B)) (R_B²/R_E²)^(n+2)` in mGal². Defaults `A = 425.28`, `B = 24`, `R_E = 6371 km` and `R_B = 6369.78 km`, so that `s = R_B²/R_E² = 0.999617` (Tscherning and Rapp 1974, OSU Report 208, Table 7). `R_B = 6370 km` would give `s = 0.999686` and is not the published value; at nmax = 2190 it changes the omission error of N by 2.5 mm (25.017 → 22.507 mm).
  - `B` choice: `B = 24` (default) is the original Tscherning–Rapp (1974) form. `B = 4` is the variant used with least-squares collocation (GRAVSOFT) when a high-degree model has been removed. Effect on the omission error of the undulation N (tail to 10800, γ = 9.81): at nmax = 2190, 22.507 mm (B = 24) vs 22.587 mm (B = 4), a difference of 0.080 mm; at nmax = 360, 227.4 mm vs 231.7 mm, a difference of 4.3 mm. So `B` only matters for low truncation degrees. Set `TscherningRapp::B` to switch. The B = 4 variant is the one used with least-squares collocation (GRAVSOFT).
- `degree_variances_with_tail(model, n_end, p)`: model variances up to `model.nmax`, Tscherning–Rapp beyond.
- `omission_error(sigma2, nmax, f, GM, R, r, gamma)`: `sqrt(Σ_{n>nmax} w_n² σ_n²)` for any `ErrorFunctional`.
- `omission_tscherning_rapp_undulation(nmax, n_tail, gamma, p)`: `OE(N) = (R_E/γ) sqrt(Σ c_n/(n−1)²)`, with `n_tail ≈ 10800`.
- `omission_kaula_undulation(R, nmax)`: Kaula rule, `R·1e-5/(nmax+1)`.

```cpp
double oe_tr = omission_tscherning_rapp_undulation(360, 10800, 9.81);   // 0.232 m
double oe_kaula = omission_kaula_undulation(6371000.0, 360);             // 0.177 m
```

## Numerical error helpers

- `relative_error(calc, ref)`: `|calc − ref| / |ref|`.
- `digit_loss_bound(M, theta)`: `log10(M) + log10(10 + 2|θ cot θ|)`, upper bound of the digits lost by the recurrence (Wittwer et al. 2008, Eq. 3, attributed there to Smith, Olver and Lozier 1981).

## Validation

`tests/test_errors.cpp`: a single-coefficient hand formula; a Monte Carlo (3000 random coefficient sets synthesised with `compute_functionals`) whose sample standard deviation matches the commission error within 6 % for T, ζ and Δg; Parseval checks for the global commission and omission errors; and consistency between the generic and direct Tscherning–Rapp routes. Values for EGM2008 at n = 360: commission 0.078 m, omission 0.232 m (Tscherning–Rapp), 0.177 m (Kaula).

## Error report

Header `src/cpp/gsh_expansion/error_report.hpp` (implementation `error_report.cpp`). `make_error_report` gathers, for one run, the numerical error, the omission error and the commission error of the undulation N:

```cpp
#include "error_report.hpp"
Model sigma = load_gfc_sigma("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
ReportOptions opt;                       // optional: colatitude, n_tail, coefficient model (load_gfc), TR parameters
ErrorReport rep = make_error_report(sigma, 360, 6378136.3, 9.81, calc_ref, &pts, opt);
set_relative_error(rep, 3e-14);          // or pass a precomputed max relative error instead of calc_ref
std::string text = format_error_report(rep);
std::string json = error_report_json(rep);
```

- Numerical: `E_rel` from `(calc, ref)` pairs, whether `E_rel <= 1e-12`, the digits `-log10(E_rel)` (capped at 16), and `digit_loss_bound(nmax, colatitude)`.
- Omission: Kaula and Tscherning–Rapp for degrees above `nmax` up to `n_tail`, on the model sphere `R`. If `opt.coefficients` (a `load_gfc` model) has degree above `nmax`, the Tscherning–Rapp value also uses its own degree variances.
- Commission: global `CE(N)` (`commission_global`, ZetaEll, degrees 2 to min(`nmax`, sigma degree)), and point min/mean/max when `pts` is given. Point values use the per-point normal gravity; the global value uses `gamma`.
