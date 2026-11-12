# ellipsoid: normal field

Header `src/cpp/gsh_utils/ellipsoid.hpp`, implementation `src/cpp/gsh_utils/ellipsoid.cpp`. Specification §6 (internal study notes, not distributed).

> **Examples.** Each example assumes `src/cpp` is on the include path (`#include "gsh_utils/ellipsoid.hpp"`) and uses the `gsh` namespace. Numbers marked "~" are hand-derived approximations; the others are exact values taken from the source.

## `struct Ellipsoid`

| Field | GRS80 |
|---|---|
| `a` | 6378137.0 m |
| `e2` | 0.00669438002290 |
| `GM` | 3.986005e14 m³/s² |
| `J2` | 1.08263e-3 |
| `omega` | 7.292115e-5 rad/s |

`Ellipsoid::GRS80()` returns these values. The fields can be changed to use another ellipsoid.

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include <cstdio>

int main() {
  gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  std::printf("a = %.1f m, e2 = %.14f, GM = %.6e\n", e.a, e.e2, e.GM);
  // Output: a = 6378137.0 m, e2 = 0.00669438002290, GM = 3.986005e+14
  e.GM = 3.986004418e14;  // use a different GM (e.g. the WGS84 value)
}
```

### `double zonal(int k) const`
Coefficient `C̄^Ell_{2k,0}`, with k in 0..10 (PDF 52-53). Throws `out_of_range` outside that interval.

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include <cstdio>
#include <stdexcept>

int main() {
  const gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  std::printf("C_0,0 = %.1f\n", e.zonal(0));  // 1.0 by definition
  std::printf("C_2,0 = %.6e\n", e.zonal(1));  // ~ -4.8416e-04 = -J2/sqrt(5)
  try { e.zonal(11); } catch (const std::out_of_range &) { std::puts("k = 11 rejected"); }
}
```

Do not set `e2 = 0` and then call `zonal(k)` for k >= 2: the formula divides by `e2`, so the result is NaN.

### `double radius(double lat) const`
Radius of the ellipsoid at spherical latitude φ:

```text
r_ell = a √((1 − e²) / (1 − e² cos²φ))
```

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include <cstdio>

int main() {
  const gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  std::printf("r(0)    = %.1f m\n", e.radius(0.0));              // a = 6378137.0 m
  std::printf("r(pi/4) = %.1f m\n", e.radius(gsh::PI / 4));      // ~ 6367418 m
  std::printf("r(pi/2) = %.1f m\n", e.radius(gsh::PI / 2));      // b = a sqrt(1 - e2) ~ 6356752 m
}
```

### `double geodetic_of_spherical(double lat) const`
Geodetic latitude of the point **on the ellipsoid** with spherical latitude φ: `atan2(sinφ, (1−e²) cosφ)`.

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include <cstdio>

int main() {
  const gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  const double phi = e.geodetic_of_spherical(gsh::PI / 4);  // input and output in radians
  std::printf("geodetic = %.6f rad = %.4f deg\n", phi, gsh::rad2deg(phi));
  // Output: ~ 0.78876 rad = ~45.19 deg (the geodetic latitude is larger than the spherical one)
  std::printf("pole: %.6f rad\n", e.geodetic_of_spherical(gsh::PI / 2));  // exactly pi/2
}
```

### `Normal normal_gravity(double lat, double r) const`
Normal gravity from a spherical-harmonic expansion (γ_SH), with zonal terms up to degree 20 and the centrifugal term.

- **Input:** spherical latitude in radians, radius in metres.
- **Returns:** `Normal{g_x, g_z, gamma}` in m/s². `g_y` is 0 by symmetry.
- **Poles:** no singularity; uses `dP_n0/dφ = √(n(n+1)/2) P_n1`.

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include <cstdio>

int main() {
  const gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  // Equator: lat = 0, r = a.
  const gsh::Ellipsoid::Normal eq = e.normal_gravity(0.0, e.a);
  // North pole: lat = pi/2, r = radius at the pole (b).
  const gsh::Ellipsoid::Normal pole = e.normal_gravity(gsh::PI / 2, e.radius(gsh::PI / 2));
  std::printf("equator: g_x = %.6f  g_z = %.6f  gamma = %.6f m/s^2\n", eq.g_x, eq.g_z, eq.gamma);
  std::printf("pole:    gamma = %.6f m/s^2\n", pole.gamma);
  // Expected: equator gamma ~ 9.7803 (close to gamma_e), g_x = 0 by symmetry;
  //           pole gamma ~ 9.8322 (close to gamma_p).
  // Sign note: g_z is negative (~ -9.780 at the equator) because dVr = -GM/r^2 * (positive sum).
  // The header calls g_z the "up" component, so check the sign convention before combining it.
}
```

### `double somigliana(lat_geodetic, gamma_e = 9.7803267715, gamma_p = 9.8321863685) const`
Somigliana closed-form normal gravity on the ellipsoid, for a **geodetic** latitude in radians. The defaults are the GRS80 equatorial and polar values; they are not derived from the ellipsoid members, so pass the right values if you change `a`, `GM`, `J2` or `omega`. Used only to validate `normal_gravity`.

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include <cstdio>

int main() {
  const gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  std::printf("%.10f\n", e.somigliana(0.0));     // 9.7803267715 (gamma_e, the default)
  std::printf("%.10f\n", e.somigliana(gsh::PI / 2));  // 9.8321863685 (gamma_p)
  std::printf("%.4f\n", e.somigliana(gsh::PI / 4));   // ~ 9.8062 m/s^2 at 45 deg geodetic
  // Custom constants are passed explicitly:
  std::printf("%.10f\n", e.somigliana(0.0, 9.78, 9.83));  // 9.78 (equatorial value replaced)
}
```

## `Model disturbing_model(const Model &model, const Ellipsoid &ell)`

Difference coefficients for T (PDF 50-51), scaled to the GM and R of the model. Only the even zonal terms (n = 2k, m = 0, k = 0..10 and 2k ≤ nmax) are changed:

```text
ΔC̄_nm = C̄_nm − C̄^Ell_nm (GM^Ell/GM) (a/R)^n
ΔS̄_nm = S̄_nm
```

With `nmin = 0` the degree-zero term is included, `ΔC̄₀₀ = (GM − GM^Ell)/GM`. With EGM2008 and GRS80 this equals −0.937 m in N. See [results](../resultados.md).

**Example**

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include "gsh_utils/model.hpp"
#include <cstdio>

int main() {
  // Stand-in model with GM and R of the file, nmax = 20, and only C00 = 1.
  // For real data use gsh::load_gfc("path/to/model.gfc", 20) (needs the file on disk).
  gsh::Model m = gsh::Model::zeros(3.986004418e14, 6378137.0, 20);
  m.C[gsh::idx_tri(0, 0, 20)] = 1.0;
  const gsh::Ellipsoid e = gsh::Ellipsoid::GRS80();
  const gsh::Model d = gsh::disturbing_model(m, e);
  std::printf("dC00 = %.3e\n", d.c(0, 0));  // 1 - GM_ell/GM ~ -1.46e-07
  std::printf("dC20 = %.4e\n", d.c(2, 0));  // 0 - C20_ell * GM_ell/GM * (a/R)^2 ~ +4.8416e-04
  std::printf("dC21 = %.4e\n", d.c(2, 1));  // S-type and m > 0 terms are untouched: 0
}
```

All other coefficients, and every S coefficient, are copied unchanged. The ellipsoid has no odd-degree or m > 0 terms, so for those the formula reduces to the copy.

## Putting it together

A typical workflow: take the GRS80 ellipsoid, compute the disturbing model, and check the normal gravity at a surface point against the closed form.

```cpp
#include "gsh_utils/ellipsoid.hpp"
#include "gsh_utils/model.hpp"
#include <cstdio>

int main() {
  // 1. Reference ellipsoid and a geopotential model (stand-in here; use load_gfc for real data).
  const gsh::Ellipsoid ell = gsh::Ellipsoid::GRS80();
  gsh::Model full = gsh::Model::zeros(ell.GM, ell.a, 20);
  full.C[gsh::idx_tri(0, 0, 20)] = 1.0;

  // 2. Disturbing coefficients: the ellipsoidal zonal part removed (here all ~0 since GM and a match).
  const gsh::Model dist = gsh::disturbing_model(full, ell);
  std::printf("dC00 = %.3e\n", dist.c(0, 0));  // 0 exactly, because GM = GM_ell and a = R

  // 3. Normal gravity at spherical latitude 45 deg, on the surface (r = radius(lat)).
  const double lat = gsh::PI / 4;
  const double r = ell.radius(lat);                 // ~ 6367418 m
  const gsh::Ellipsoid::Normal n = ell.normal_gravity(lat, r);

  // 4. Compare with the closed form at the matching geodetic latitude.
  const double phi = ell.geodetic_of_spherical(lat);  // ~ 0.78876 rad
  std::printf("r = %.1f m, geodetic = %.4f deg\n", r, gsh::rad2deg(phi));
  std::printf("gamma_SH = %.6f  gamma_Somigliana = %.6f m/s^2\n", n.gamma, ell.somigliana(phi));
  // Expected: both ~ 9.8062 m/s^2; the difference is the validation residual.
}
```
