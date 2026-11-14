# Legendre functions

Associated Legendre functions are one of the solutions of the differential equation:

$$
\sin(\theta)g''(\theta) + \cos(\theta)g'(\theta) + \left[n(n+1)\sin(\theta)-\frac{m^2}{\sin(\theta)}\right] g(\theta) = 0
$$

Where:

$$
g(\theta) = P_{nm}(\cos\theta)
$$

Spherical harmonics are a linear combination with the general form:

$$
Y_n(\theta,\lambda) = \sum_{m=0}^{n}\left[a_{nm}P_{nm}(\cos\theta)\cos(m\lambda) + b_{nm}P_{nm}(\cos\theta)\sin(m\lambda)\right]
$$

Where $a_{nm}$ and $b_{nm}$ are arbitrary constants. This is the general expression for the surface spherical harmonics $Y_n(\theta,\lambda)$.

These relations are important for computing the solutions of Laplace's equation $\Delta V = 0$, which are:

- **Interior:**

$$
V_i(r,\theta,\lambda) = \sum_{n=0}^{\infty} r^n \sum_{m=0}^{n}\left[a_{nm}P_{nm}(\cos\theta)\cos(m\lambda) + b_{nm}P_{nm}(\cos\theta)\sin(m\lambda)\right]
$$

- **Exterior:**

$$
V_e(r,\theta,\lambda) = \sum_{n=0}^{\infty} r^{-(n+1)} \sum_{m=0}^{n}\left[a_{nm}P_{nm}(\cos\theta)\cos(m\lambda) + b_{nm}P_{nm}(\cos\theta)\sin(m\lambda)\right]
$$

Once the Legendre function $P_{nm}\left(\cos\theta\right)$ is introduced (the indices _n_ and _m_ are known as the degree and order of the function, respectively) as a solution of the differential equation above, it is convenient to use the substitution:

$$
t = \cos\left(\theta\right)
$$

Solving the differential equation then gives that the Legendre function $\bar{g}\left(t\right) = P_{nm}\left(t\right)$ is mathematically defined as:

$$
P_{nm}\left(t\right) = \frac{1}{2^nn!}\left(1-t^2\right)^{\frac{m}{2}}\frac{d^{n+m}}{dt^{n+m}}\left(t^2-1\right)^n
$$

Evaluating this function has been a substantial problem in numerical computation and remains hard to implement directly. Researchers have nevertheless proposed ways to compute it in fully normalised form (Full Normalized Associated Legendre Function) with specific computational methods, presented below.

## Computational methods

### Standard Forward Column method (SFC)

The standard forward column method is one of the most common ways to compute the normalised associated Legendre function. This method is used by the Legendre module of the library and is computed up to degree 1800 for intermediate degrees (0° - 45°). The method uses the following initial values for the iterative calculation:

$$
\bar{P}_{0,0} = 1 \quad; \quad \bar{P}_{1,0} = \sqrt{3} \cdot t \quad; \quad \bar{P}_{1,1} = \sqrt{3} \cdot u \quad; \quad t = \cos\left(\theta\right) \quad; \quad u = \sin\left(\theta\right)
$$

#### Sectorials

The sectorials serve as seed values for the recurrences presented below, and are computed from the initial values shown above:

$$
\bar{P}_{m,m}\left(\theta\right) = u\sqrt{\frac{2m+1}{2m}} \bar{P}_{m-1,m-1}\left(\theta\right), \quad m \geq 2
$$

The library module defines:

$$
d_m = \sqrt{\frac{2m+1}{2m}}
$$

with $d_1 = \sqrt{3}$, because $\bar{P}_{1,1} = \sqrt{3}\,u$.

#### First term of each column

Once the sectorial is known, the term $n = m+1$ is computed without the $b$ term:

$$
\bar{P}_{m+1,m} = e_{m+1,m} \cdot t \cdot \bar{P}_{m,m}\left(\theta\right)
$$

Where the library defines $e_{m+1,m}$ as:

$$
e_m = \sqrt{2m+3}
$$

#### Tesserals

The tesserals are computed from the second position on, so they are computed once $n \geq m+2$:

$$
\bar{P}_{n,m}\left(\theta\right) = a_{n,m} \cdot t \cdot \bar{P}_{n-1,m}\left(\theta\right) - b_{n,m} \cdot \bar{P}_{n-2,m}\left(\theta\right)
$$

Where:

$$
a_{n,m} = \sqrt{\frac{(2n - 1)(2n+1)}{(n-m)(n+m)}}
$$

$$
b_{n,m} = \sqrt{\frac{(2n+1)(n+m-1)(n-m-1)}{(n-m)(n+m)(2n-3)}}
$$

#### Derivatives

The library differentiates with respect to the **spherical latitude** $\varphi$, with $t = \sin\varphi$ and $u = \cos\varphi$. The kernels use the column form:

$$
\frac{d\bar{P}_{n,m}}{d\varphi} = \frac{1}{u}\left(f_{n,m} \cdot \bar{P}_{n-1,m} - n \cdot t \cdot \bar{P}_{n,m}\right)
$$

$$
\frac{d^2\bar{P}_{n,m}}{d\varphi^2} = \frac{t}{u}\frac{d\bar{P}_{n,m}}{d\varphi} + \left(\frac{m^2}{u^2} - n(n+1)\right)\bar{P}_{n,m}
$$

Where:

$$
f_{n,m} = \sqrt{\frac{(n^2-m^2)(2n+1)}{(2n-1)}}
$$

These forms are singular when $u \to 0$. The function `dphi` instead uses the non-singular form, which combines functions of the same degree:

$$
\frac{d\bar{P}_{n,0}}{d\varphi} = \sqrt{\frac{n(n+1)}{2}}\,\bar{P}_{n,1}
$$

$$
\frac{d\bar{P}_{n,m}}{d\varphi} = \frac{1}{2}\left[-\sqrt{(n+m)(n-m+1)}\sqrt{1+\delta_{m,1}}\,\bar{P}_{n,m-1} + \sqrt{(n-m)(n+m+1)}\,\bar{P}_{n,m+1}\right]
$$

### Modified method (MSFC)

For high degrees the sectorial falls below the range of `double`. MSFC stores each value as a mantissa and an exponent, $\bar{P} = x \cdot R^{-k}$ with $R = 2^{500}$. When $|x|$ exceeds $R$ in the column, it is divided by $R$ and $k$ is reduced. With $k \geq 3$ the value is returned as 0. It shares one exponent per column, unlike extended-range arithmetic (ERA), which carries one per number; ERA is in `era`.

## Storage of results

The triangles are stored in a contiguous vector, with contiguous columns $m$:

$$
idx(n,m,nmax) = m\left(nmax + 1\right) - \frac{m(m-1)}{2} + (n - m)
$$

$$
size(nmax) = \frac{(nmax+1)(nmax+2)}{2}
$$

The coefficients $a$, $b$ and $f$ use the same index, without offsets, so that $a_{m,m} = b_{m,m} = 0$ and $b_{m+1,m} = 0$.

The helpers below are in `core.hpp` (namespace `gsh`, `constexpr`, no bounds checks).

### `tri_size(int nmax)`

Number of entries of the packed triangle for degree `nmax`.

**Example**

```cpp
#include <cstdio>
#include <vector>
#include "core.hpp"
int main()
{
  static_assert(gsh::tri_size(2) == 6);                  // (n,m): 00 10 11 20 21 22
  static_assert(gsh::tri_size(2190) == 2401336);         // (2191 * 2192) / 2
  std::vector<double> P(gsh::tri_size(2190));            // about 19.2 MB of doubles
  std::printf("%zu\n", P.size());                        // 2401336
}
```

### `idx_tri(int n, int m, int nmax)`

Position of `(n, m)` in the packed triangle: `col_offset(m, nmax) + (n - m)`. Valid for `0 <= m <= n <= nmax`; nothing checks this.

**Example**

```cpp
#include <cstdio>
#include "core.hpp"
int main()
{
  static_assert(gsh::idx_tri(0, 0, 2) == 0);
  static_assert(gsh::idx_tri(2, 0, 2) == 2);
  static_assert(gsh::idx_tri(1, 1, 2) == 3);
  static_assert(gsh::idx_tri(2, 2, 2) == 5);
  std::printf("%lld\n", static_cast<long long>(gsh::idx_tri(3, 1, 3))); // 6: column 1 starts at 4, n = 3 is its third entry
}
```

### `col_offset(int m, int nmax)`

Index of the first element of column `m`, that is `(n, m) = (m, m)`. Its value is `m(nmax + 1) - m(m - 1)/2`.

**Example**

```cpp
#include <cstdio>
#include "core.hpp"
int main()
{
  static_assert(gsh::col_offset(0, 2) == 0);
  static_assert(gsh::col_offset(1, 2) == 3);
  static_assert(gsh::col_offset(2, 2) == gsh::idx_tri(2, 2, 2));
  std::printf("%lld\n", static_cast<long long>(gsh::col_offset(1000, 2190))); // 1691500 = 1000*2191 - 1000*999/2
}
```

## Functions

Header `src/cpp/gsh_legendre/legendre.hpp`, implementation `src/cpp/gsh_legendre/legendre.cpp`. Namespace `gsh::legendre`. Arguments: `t = sinφ = cosθ`, `u = cosφ = sinθ`. At the exact poles it receives `u = 0` and `t = ±1`, not `cos(π/2)`.

**Includes and namespace.** Include `"legendre.hpp"` (folder `src/cpp/gsh_legendre/`). That header includes `"core.hpp"`, which lives in `src/cpp/gsh_utils/`, so both folders must be on the include path. The Legendre functions are in `gsh::legendre`; `tri_size`, `idx_tri`, `col_offset`, `deg2rad` and `PI` are in `gsh`. Each example includes the headers it uses and qualifies every name. Numbers in comments are the expected results for the stated inputs, derived by hand from the recurrences unless the comment says otherwise.

### `Coefs make_coefs(int nmax)`
Computes `a`, `b`, `f` and `d`. Throws `invalid_argument` if `nmax < 0`.

**Example**

```cpp
#include <cstdio>
#include <stdexcept>
#include "legendre.hpp"
int main()
{
  const auto c = gsh::legendre::make_coefs(2100);   // O(nmax^2): build once per degree, reuse for every latitude
  std::printf("%d %zu\n", c.nmax, c.d.size());      // 2100 2101
  try { gsh::legendre::make_coefs(-1); }
  catch (const std::invalid_argument &e) { std::printf("%s\n", e.what()); } // "make_coefs: nmax negativo" (source message is in Spanish)
}
```

### `Coefs` (struct)

`Coefs` holds `nmax`, the packed triangles `a`, `b` and `f` (each `tri_size(nmax)` long, indexed like `P`) and `d` (`nmax + 1` entries, one per order `m`). It is built by `make_coefs`.

**Example**

```cpp
#include <cstdio>
#include "legendre.hpp"
int main()
{
  const auto c = gsh::legendre::make_coefs(2);
  std::printf("%.7f %.7f\n", c.d[1], c.a[gsh::idx_tri(1, 0, 2)]);                       // 1.7320508 1.7320508 (both sqrt(3))
  std::printf("%.7f %.7f\n", c.a[gsh::idx_tri(2, 0, 2)], c.b[gsh::idx_tri(2, 0, 2)]);   // 1.9364917 1.1180340 (sqrt(15/4), sqrt(5/4))
  std::printf("%.1f %.1f\n", c.b[gsh::idx_tri(1, 0, 2)], c.a[gsh::idx_tri(0, 0, 2)]);   // 0.0 0.0 (b needs n >= m+2, a(m,m) = 0)
}
```

### `void sfc(double t, double u, const Coefs &c, std::span<double> P)`
SFC method. Fast; valid up to nmax ≈ 1800.

**Example**

```cpp
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>
#include "legendre.hpp"
int main()
{
  const auto c = gsh::legendre::make_coefs(2);
  std::vector<double> P(gsh::tri_size(2)), small(5);
  const double lat = gsh::deg2rad(30.0);                  // geocentric latitude
  gsh::legendre::sfc(std::sin(lat), std::cos(lat), c, P);
  std::printf("%.7f %.7f %.7f\n", P[gsh::idx_tri(0, 0, 2)], P[gsh::idx_tri(1, 0, 2)], P[gsh::idx_tri(1, 1, 2)]); // 1 0.8660254 1.5
  std::printf("%.7f %.7f %.7f\n", P[gsh::idx_tri(2, 0, 2)], P[gsh::idx_tri(2, 1, 2)], P[gsh::idx_tri(2, 2, 2)]); // -0.2795085 1.6770510 1.4523688
  try { gsh::legendre::sfc(0.5, 0.8, c, small); }        // 5 < tri_size(2) = 6
  catch (const std::invalid_argument &) { std::printf("buffer too small\n"); }
}
```

### `void era(double t, double u, const Coefs &c, std::span<double> P)`
Extended-range arithmetic (Smith-Olver-Lozier; equivalent to `fnALFs = 3` in GrafLab). Each number is a mantissa in [0.5, 1) plus its own integer exponent, renormalised after every operation and aligned before subtraction. No scale factors, and valid at any latitude and degree. It is the precision reference for `msfc`; for ALFs alone it costs ~1.2× MSFC (nmax 2190: 17 ms versus 12 ms), but it is not used in the full synthesis. Same error as MSFC against a 60-digit oracle (mpmath): ≤ 1e-12 relative at mid latitudes and up to nmax 10800; up to ~1e-10 at the least favourable point 0.0001° from the pole (zeros of P included), which is the limit of the recurrence in double precision, not of the range.

**Example**

```cpp
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const auto c = gsh::legendre::make_coefs(2);
  std::vector<double> P(gsh::tri_size(2));
  gsh::legendre::era(0.0, 1.0, c, P);                     // equator: t = 0, u = 1
  std::printf("%.7f %.7f %.7f\n", P[gsh::idx_tri(1, 1, 2)], P[gsh::idx_tri(2, 0, 2)], P[gsh::idx_tri(2, 2, 2)]); // 1.7320508 -1.1180340 1.9364917
  std::printf("%.1f\n", P[gsh::idx_tri(2, 1, 2)]);      // 0.0: P(2,1) carries a factor t
}
```

### `void msfc(double t, double u, const Coefs &c, std::span<double> P)`
MSFC method. Valid at ultra-high degree; tested at nmax = 2700 with u = 0.34.

**Example**

```cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const int nmax = 2700;
  const auto c = gsh::legendre::make_coefs(nmax);
  std::vector<double> P(gsh::tri_size(nmax));
  const double u = 0.34, t = std::sqrt(1.0 - u * u);      // the case used in the doc's test
  gsh::legendre::msfc(t, u, c, P);
  const int n = 2000;                                      // above the SFC limit of 1800
  double s = 0.0;
  for (int m = 0; m <= n; ++m) { const double v = P[gsh::idx_tri(n, m, nmax)]; s += v * v; }
  std::printf("sum = %.12f (expected %d)\n", s, 2 * n + 1); // close to 4001; not hand-derived, the sum rule is the check
}
```

### `void compute(double t, double u, const Coefs &c, std::span<double> P)`
Uses `sfc` if `nmax ≤ 1800`, and `msfc` otherwise.

**Example**

```cpp
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const auto small = gsh::legendre::make_coefs(2);        // compute() takes the sfc branch
  const auto big = gsh::legendre::make_coefs(2000);       // compute() takes the msfc branch
  std::vector<double> P(gsh::tri_size(2)), Q(gsh::tri_size(2000));
  gsh::legendre::compute(0.0, 1.0, small, P);             // equator
  gsh::legendre::compute(0.0, 1.0, big, Q);
  std::printf("%.7f %.7f\n", P[gsh::idx_tri(1, 1, 2)], Q[gsh::idx_tri(1, 1, 2000)]); // 1.7320508 1.7320508
}
```

### `void dphi(int nmax, std::span<const double> P, std::span<double> dP)`
Computes $d\bar{P}/d\varphi$ from the triangle, without singularity at the poles. Applied twice it gives the second derivative. Throws `invalid_argument` if the buffers are too small.

**Example**

```cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const int nmax = 2;
  const auto c = gsh::legendre::make_coefs(nmax);
  const double lat = gsh::deg2rad(30.0);
  std::vector<double> P(gsh::tri_size(nmax)), dP(P.size()), ddP(P.size());
  gsh::legendre::sfc(std::sin(lat), std::cos(lat), c, P);
  gsh::legendre::dphi(nmax, P, dP);                        // dP/dphi
  gsh::legendre::dphi(nmax, dP, ddP);                      // d2P/dphi2
  std::printf("%.7f %.7f\n", dP[gsh::idx_tri(1, 0, nmax)], dP[gsh::idx_tri(1, 1, nmax)]);  // 1.5000000 -0.8660254
  std::printf("%.7f %.7f\n", dP[gsh::idx_tri(2, 0, nmax)], ddP[gsh::idx_tri(1, 0, nmax)]); // 2.9047375 -0.8660254
}
```

```cpp
const int nmax = 10, n = 4, m = 2;
const double lat = 0.5;  // geocentric latitude [rad]
const auto c = gsh::legendre::make_coefs(nmax);
std::vector<double> P(gsh::tri_size(nmax)), dP(P.size()), ddP(P.size());
gsh::legendre::compute(std::sin(lat), std::cos(lat), c, P);
gsh::legendre::dphi(nmax, P, dP);
gsh::legendre::dphi(nmax, dP, ddP);
const double valor = P[gsh::idx_tri(n, m, nmax)];
```

`sfc`, `msfc`, `era` and `compute` throw `invalid_argument` if `P` is smaller than `tri_size(nmax)`. The internal function `check` performs that test.

## Poles and the numerical limit near them

**Analytic pole shortcut.** `sfc`, `msfc` and `era` check `|(|t| − 1)| < 1e-15` and `|u| < 1e-15`. If both hold, they fill the triangle directly instead of running the recurrence:

$$
\bar{P}_{n,0}(\pm 1) = (\pm 1)^n \sqrt{2n+1}, \qquad \bar{P}_{n,m}(\pm 1) = 0 \quad (m \ge 1).
$$

The condition on `u` is an addition to the one on `t`. `|t| = 1 − 1e-15` already means `u = 4.5e-8` (about 0.3 m from the pole), where $\bar{P}_{n,1}$ is still of order 1e-3 at n = 2700, so the pole values are used only when `u` itself is negligible. `cos(π/2) = 6e-17` passes. The shortcut makes the sum rules exact at φ = ±90° (before it, the relative error was 3.5e-11 at M = 1800 and 9.0e-11 at M = 2700).

**Example.** At the north pole the shortcut fills the triangle exactly, and the south pole flips the sign of the odd-n zonal terms:

```cpp
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const auto c = gsh::legendre::make_coefs(2);
  std::vector<double> P(gsh::tri_size(2));
  gsh::legendre::sfc(1.0, 0.0, c, P);                       // north pole: t = 1, u = 0
  std::printf("%.7f %.7f %.7f\n", P[gsh::idx_tri(0, 0, 2)], P[gsh::idx_tri(1, 0, 2)], P[gsh::idx_tri(2, 0, 2)]); // 1 1.7320508 2.2360680
  std::printf("%.1f %.1f\n", P[gsh::idx_tri(1, 1, 2)], P[gsh::idx_tri(2, 1, 2)]);  // 0.0 0.0 (m >= 1 vanish)
  gsh::legendre::sfc(-1.0, 0.0, c, P);                      // south pole
  std::printf("%.7f\n", P[gsh::idx_tri(1, 0, 2)]);          // -1.7320508
}
```

**Numerical limit for 86.5° ≤ |φ| < 89.5°.** Close to the pole, but not at it, the zonal column is still computed with the three-term recurrence $\bar P_n = a\,t\,\bar P_{n-1} - b\,\bar P_{n-2}$. At $t = \pm 1$ its characteristic equation has a double root, so near the pole the rounding error of each step is carried almost linearly through the remaining steps and the accumulated relative error grows as $\mathcal{O}(M^2 \epsilon)$, with $\epsilon = 2^{-52}$. The slight loss of precision in this strip (sum-rule errors of 1e-12 to 6e-12) is that mathematical accumulation, not an instability of the kernel: values stay finite and bounded, and the measured error follows the $M^2$ law with $\text{err}/(M^2\epsilon) \le 0.0051$ for every derivative order $d = 0 \ldots 3$. A plain double-precision recurrence reproduces the same numbers, so the cause is the recurrence and not the implementation. The acceptance test (`benchmark/alf/test_alf_identities.cpp`) allows `max(base, 0.02·M²·ε)` for 85° ≤ |φ| < 90° (`get_allowed_tolerance`); the base thresholds are 1e-12 for M < 1800, 1e-10 for 1800 ≤ M < 2190 and 1e-9 for M ≥ 2190. The kernels of the fast synthesis (`kernels/column_kernel.inl`) do not use this shortcut.

## Putting it together

Compute one triangle at 30 degrees geocentric latitude, take the first and second latitude derivatives, and check the sum rule for two degrees:

```cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const int nmax = 20;
  const double phi = gsh::deg2rad(30.0);                   // geocentric latitude
  const auto c = gsh::legendre::make_coefs(nmax);
  std::vector<double> P(gsh::tri_size(nmax)), dP(P.size()), ddP(P.size());
  gsh::legendre::compute(std::sin(phi), std::cos(phi), c, P);
  gsh::legendre::dphi(nmax, P, dP);                        // dP/dphi
  gsh::legendre::dphi(nmax, dP, ddP);                      // d2P/dphi2

  // Sum rule: for each n, sum over m of Pbar_nm^2 = 2n+1
  for (int n : {10, 20})
  {
    double s = 0.0;
    for (int m = 0; m <= n; ++m)
      s += P[gsh::idx_tri(n, m, nmax)] * P[gsh::idx_tri(n, m, nmax)];
    std::printf("n=%d sum=%.12f expected=%d\n", n, s, 2 * n + 1); // sum close to 21 (n = 10) and 41 (n = 20)
  }
  // Pbar_10 = sqrt(3) sin(phi) = 0.8660254; dP = sqrt(3) cos(phi) = 1.5; ddP = -sqrt(3) sin(phi) = -0.8660254
  std::printf("%.7f %.7f %.7f\n", P[gsh::idx_tri(1, 0, nmax)], dP[gsh::idx_tri(1, 0, nmax)], ddP[gsh::idx_tri(1, 0, nmax)]);
}
```

## Batched extraction for many points

`compute()` returns the whole triangle for one point. At degree 5540 one triangle is about 123 MB, so the triangles of many points cannot be stored. `extract()` (C++) and `legendre_points()` (Python) compute each point's triangle in a buffer owned by one thread, reuse that buffer for every point of the thread, and copy out only the requested orders `m` and degrees `n`.

**Memory of the result.** Let `count` be the number of stored columns:

$$
count = \sum_{m \in \text{orders}} \left(nmax - \max(m, nmin) + 1\right)
$$

Then the result takes `bytes = 8 * npoints * count` for `P`, and twice that with `derivative`. Orders are stored in ascending `m`. Output larger than 2^31 doubles (16 GB) is refused with `std::length_error` (C++) or `ValueError` (Python); the message says how to reduce the request. The limit is the constant `extract_max_doubles` in `legendre.hpp`.

**Selection rules.** `orders` must be unique and lie in `0..nmax`; empty means all orders. `0 <= nmin <= nmax <= coefs.nmax`, and `nmax = -1` means `coefs.nmax`. Latitudes are geocentric, in radians. The result is bit-identical for any thread count.

**C++ example**

```cpp
#include <cstdio>
#include <vector>
#include "legendre.hpp"
int main()
{
  const auto c = gsh::legendre::make_coefs(5540);        // build once per degree
  const std::vector<double> lat = {0.0, 0.5, -1.2};      // geocentric latitudes [rad]
  gsh::legendre::Selection sel;
  sel.orders = {0, 1, 100};                              // m = 0, 1 and 100
  sel.nmin = 2;                                          // degrees max(m, 2)..5540
  sel.derivative = true;                                 // also dP/dphi
  const auto ex = gsh::legendre::extract(lat, c, sel, 0); // threads = 0: OpenMP default
  // count = (5540-2+1) + (5540-2+1) + (5540-100+1) = 16519
  std::printf("count %lld\n", static_cast<long long>(ex.count));
  std::printf("P(point 1, n=%d, m=%d) = %.12e\n", ex.n[0], ex.m[0], ex.P[1 * ex.count + 0]);
  std::printf("dP(point 1, n=%d, m=%d) = %.12e\n", ex.n[0], ex.m[0], ex.dP[1 * ex.count + 0]);
}
```

**Python example**

```python
import math
import numpy as np
import gshtools as gsh

c = gsh.legendre.make_coefs(5540)                  # optional: reuse across calls
lat = np.radians([0.0, 30.0, -70.0])               # geocentric latitudes [rad]
r = gsh.legendre_points(lat, 5540, orders=[0, 1, 100], nmin=2, derivative=True, coefs=c)
print(r.P.shape)                                   # (3, count)
print(r.n[:3], r.m[:3])                            # int32 arrays of length count
print(r.dP.shape)                                  # same shape as P; None without derivative=True
```

The returned object is a `LegendrePoints` with the attributes `P`, `dP`, `n` and `m`. The values equal `compute()` and `dphi()` on the full triangle, entry by entry.
