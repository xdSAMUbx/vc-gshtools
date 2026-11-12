# vc-gshtools documentation

Reference for the modules included in this repository. The benchmark is documented in
[benchmark/shs/results/Summary.md](../benchmark/shs/results/Summary.md).

| Module | Header(s) | Content |
|---|---|---|
| [core](modules/core.md) | `gsh_utils/core.hpp` | Constants, triangular indices, aligned memory, coordinates |
| [model](modules/model.md) | `gsh_utils/model.hpp` | Coefficient model and the `.gfc` reader |
| [table_io](modules/table_io.md) | `gsh_utils/table_io.hpp` | Delimited-table reader |
| [ellipsoid](modules/ellipsoid.md) | `gsh_utils/ellipsoid.hpp` | Reference ellipsoid, normal gravity, difference coefficients |
| [errors](modules/errors.md) | `gsh_expansion/errors.hpp`, `error_report.hpp` | Commission and omission errors, error report |
| [height](modules/height.md) | `gsh_utils/height.hpp` | Topographic height source |
| [legendre](modules/legendre.md) | `gsh_legendre/legendre.hpp` | Legendre functions: SFC, MSFC, ERA, and their derivatives |
| [kernels](modules/kernels.md) | `gsh_legendre/kernels/` | Internal SIMD kernels of the CPU backend |
| [shs](modules/shs.md) | `gsh_expansion/shs.hpp` | Synthesis at scattered points and the reference oracle |
| [grid](modules/grid.md) | `gsh_expansion/grid.hpp` | Synthesis on latitude-longitude grids with FFTW |
| [functionals](modules/functionals.md) | `gsh_expansion/functionals.hpp` | Functionals of the gravity field |


Build and test instructions are in the [README](../README.md).
