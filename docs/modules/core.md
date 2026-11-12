# core

Header `src/cpp/gsh_utils/core.hpp`, implementation `src/cpp/gsh_utils/core.cpp`. Namespace `gsh`.

Examples in this section use the include root `src/cpp`, so the header is included as `"gsh_utils/core.hpp"`. They are not part of the build; numeric results are derived by hand from the formulas above, not from a program run.

## Constants

| Name | Value | Use |
|---|---|---|
| `PI` | π | |
| `RBITS` | 500 | Bits of the MSFC rescaling |
| `R_BIG` | 2^500 | Exponent base: value = `x · R^-k` |
| `R_INV` | 2^-500 | Inverse of `R_BIG` |

**Example**

```cpp
#include <cassert>
#include <cmath>
#include "gsh_utils/core.hpp"

int main() {
  double pi = gsh::PI;                                  // 3.141592653589793
  assert(gsh::RBITS == 500);
  assert(gsh::R_BIG * gsh::R_INV == 1.0);               // exact: both are powers of two
  assert(gsh::R_BIG == std::ldexp(1.0, gsh::RBITS));    // R = 2^500
  double x = 3.0;
  double scaled = x * gsh::R_BIG;                       // magnitude 3 * 2^500, still finite
  double back = scaled * gsh::R_INV;                    // exactly 3.0, no rounding
  return (back == x && pi > 3.14) ? 0 : 1;
}
```

## Functions

### `double deg2rad(double v)`, `double rad2deg(double v)`
Conversion between degrees and radians.

**Example**

```cpp
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  double r180 = gsh::deg2rad(180.0);          // PI, about 3.141592653590 rad
  double r45  = gsh::deg2rad(45.0);           // PI/4, about 0.785398163397 rad
  double d90  = gsh::rad2deg(gsh::PI / 2.0);  // 90.0 deg
  std::printf("%.12f %.12f %.1f\n", r180, r45, d90);
  // Output: 3.141592653590 0.785398163397 90.0
}
```

### `std::int64_t tri_size(int nmax)`
Number of elements in the triangle: `(nmax+1)(nmax+2)/2`. For nmax = 2190 it is 2 401 336.

**Example**

```cpp
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  // Number of (n, m) pairs with 0 <= m <= n <= nmax
  std::printf("%lld\n", (long long)gsh::tri_size(0));     // 1
  std::printf("%lld\n", (long long)gsh::tri_size(2));     // 6  = 3*4/2
  std::printf("%lld\n", (long long)gsh::tri_size(2190));  // 2401336
  // Storage for a full coefficient set at nmax = 2190 in doubles: 2401336 * 8 bytes, about 19.2 MB
}
```

### `std::int64_t col_offset(int m, int nmax)`
Position of the first element of column m (the element `n = m`): `m(nmax+1) − m(m−1)/2`.

**Example**

```cpp
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  const int nmax = 2;
  // Column 0 holds (0,0), (1,0), (2,0); column 1 holds (1,1), (2,1); column 2 holds (2,2).
  for (int m = 0; m <= nmax; ++m)
    std::printf("col_offset(%d, %d) = %lld\n", m, nmax, (long long)gsh::col_offset(m, nmax));
  // Output: col_offset(0, 2) = 0, col_offset(1, 2) = 3, col_offset(2, 2) = 5
  // The last column has one element, so the total is 5 + 1 = 6 = tri_size(2).
}
```

### `std::int64_t idx_tri(int n, int m, int nmax)`
Index of `(n, m)`: `col_offset(m, nmax) + (n − m)`. Columns are contiguous because the recurrence and the kernels traverse them in order. Requires `0 ≤ m ≤ n ≤ nmax` and does not check it.

**Example**

```cpp
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  const int nmax = 2;
  for (int m = 0; m <= nmax; ++m)
    for (int n = m; n <= nmax; ++n)
      std::printf("(n=%d, m=%d) -> %lld\n", n, m, (long long)gsh::idx_tri(n, m, nmax));
  // Output: (0,0)->0, (1,0)->1, (2,0)->2, (1,1)->3, (2,1)->4, (2,2)->5
  // The indices cover 0 .. tri_size(nmax) - 1 with no gaps.

  // Last element for nmax = 2190 is tri_size(2190) - 1
  std::printf("%lld\n", (long long)gsh::idx_tri(2190, 2190, 2190));  // 2401335
}
```

### `Spherical geodetic_to_spherical(lat_deg, lon_deg, h, a, e2)`
Converts geodetic coordinates to spherical coordinates on the ellipsoid `(a, e2)`.

- **Input:** latitude and longitude in degrees, ellipsoidal height in metres.
- **Returns:** `Spherical{lat, lon, r}`, in radians and metres.
- **Formula:** `N = a/√(1−e² sin²φ)`, `p = (N+h) cosφ`, `z = (N(1−e²)+h) sinφ`; then `lat = atan2(z, p)` and `r = hypot(p, z)`.

**Example**

```cpp
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  // WGS84 ellipsoid
  const double a  = 6378137.0;            // semi-major axis [m]
  const double f  = 1.0 / 298.257223563;  // flattening
  const double e2 = f * (2.0 - f);        // first eccentricity squared

  // Equator, on the ellipsoid: lat = 0, lon = 0, r = a
  gsh::Spherical s1 = gsh::geodetic_to_spherical(0.0, 0.0, 0.0, a, e2);
  // s1.lat = 0 rad, s1.lon = 0 rad, s1.r = 6378137.0 m

  // Equator, 1000 m above the ellipsoid, lon = 90 deg: r = a + h
  gsh::Spherical s2 = gsh::geodetic_to_spherical(0.0, 90.0, 1000.0, a, e2);
  // s2.lon = PI/2 rad (1.5707963), s2.r = 6379137.0 m

  // North pole: lat = PI/2 rad, r = semi-minor axis b = a(1-f), about 6356752.31 m
  gsh::Spherical s3 = gsh::geodetic_to_spherical(90.0, 0.0, 0.0, a, e2);
  std::printf("%.7f %.2f\n", s3.lat, s3.r);  // Output: 1.5707963 6356752.31
  (void)s1;
  (void)s2;
}
```

The geocentric latitude differs from the geodetic one away from the equator and the poles. For example, at 45 deg geodetic latitude the geocentric latitude is about 44.8 deg.

## Types

### `Spherical`
Fields `lat` and `lon` (radians) and `r` (metres).

**Example**

```cpp
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  gsh::Spherical s{gsh::deg2rad(40.0), gsh::deg2rad(-3.0), 6370000.0};
  // s.lat = 0.698132 rad (40 deg), s.lon = -0.052360 rad (-3 deg), s.r = 6370000.0 m
  std::printf("lat=%.6f lon=%.6f r=%.1f\n", s.lat, s.lon, s.r);
  // Output: lat=0.698132 lon=-0.052360 r=6370000.0
  gsh::Spherical copy = s;  // trivially copyable plain struct
  (void)copy;
}
```

### `AlignedAllocator<T>` and `AlignedVec<T>`
`AlignedVec<T>` is a `std::vector` aligned to 64 bytes (cache line and ZMM register). It is used for all numerical arrays.

**Example** (`AlignedAllocator<T>`)

```cpp
#include <cstdint>
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  gsh::AlignedAllocator<double> alloc;
  double *p = alloc.allocate(1000);  // 8000 bytes, returned on a 64-byte boundary
  std::printf("aligned: %d\n", reinterpret_cast<std::uintptr_t>(p) % 64 == 0);  // aligned: 1
  alloc.deallocate(p, 1000);
  // Throws std::bad_alloc if the underlying aligned allocation fails.
}
```

**Example** (`AlignedVec<T>`)

```cpp
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  gsh::AlignedVec<double> v(1000, 0.0);  // behaves like std::vector<double>
  v[0] = 1.5;
  std::printf("size=%zu aligned=%d\n", v.size(),
              reinterpret_cast<std::uintptr_t>(v.data()) % 64 == 0);
  // Output: size=1000 aligned=1
  gsh::AlignedVec<double> copy = v;      // copy keeps the aligned allocator
  std::printf("copy[0]=%.1f\n", copy[0]);  // copy[0]=1.5
}
```

## Putting it together

A typical workflow: place a station on the ellipsoid, allocate a coefficient set for degree 2190, and address one coefficient by its triangular index.

```cpp
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include "gsh_utils/core.hpp"

int main() {
  const int nmax = 2190;
  const double a = 6378137.0, f = 1.0 / 298.257223563, e2 = f * (2.0 - f);

  // 1. Coefficient storage for all (n, m) pairs up to nmax, zero-filled (about 19.2 MB).
  gsh::AlignedVec<double> cnm(static_cast<std::size_t>(gsh::tri_size(nmax)), 0.0);

  // 2. Station at 40.4 N, 3.7 W, 650 m ellipsoidal height. Inputs in degrees and metres.
  gsh::Spherical sp = gsh::geodetic_to_spherical(40.4, -3.7, 650.0, a, e2);
  std::printf("lat=%.6f rad  lon=%.6f rad  r=%.1f m\n", sp.lat, sp.lon, sp.r);
  // r is close to 6.37e6 m; lat is about 0.7018 rad (geocentric, slightly below the 40.4 deg geodetic value).

  // 3. Address the normalized coefficient C_20 (n = 2, m = 0). Index = 2.
  const std::int64_t k20 = gsh::idx_tri(2, 0, nmax);
  cnm[static_cast<std::size_t>(k20)] = -4.84165e-4;  // illustrative value, typical of a C20 term

  // 4. Walk the first column (m = 0) in order; its indices are 0, 1, 2, ..., nmax.
  for (int n = 0; n <= 2; ++n)
    std::printf("C(%d,0) at index %lld\n", n, (long long)gsh::idx_tri(n, 0, nmax));
  // Output: C(0,0) at index 0, C(1,0) at index 1, C(2,0) at index 2
  return 0;
}
```
