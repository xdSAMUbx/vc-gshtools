# Examples

Small programs that show the public API of the `gsh` library (`src/cpp/gsh_utils`, `src/cpp/gsh_legendre` and
`src/cpp/gsh_expansion`). Each example is one `.cpp` file in this folder, with a header comment that says what it
shows, the command line and what the output means. They are built with the other targets of the project
(see the main README); the executables go to `build/<preset>/bin`.

Most examples read a fully normalised ICGEM `.gfc` model. The model file is not part of the repository: any
fully normalised model works, and EGM2008 (degree 2190, with formal errors in the file) was used for development.
`ex_errors` needs the formal-error columns (`sigma C`, `sigma S`), which EGM2008 has.

## Index

| Example | What it shows | How to run | Runtime (estimate, not measured) | Module |
|---|---|---|---|---|
| `ex_puntual.cpp` | Fully normalised Legendre functions (MSFC) at one latitude and degree | `ex_puntual [lat_deg=10] [nmax=2699]` | seconds | `gsh_legendre` |
| `ex_legendre.cpp` | sfc, msfc, era and compute() agreement, the sum-rule identity, and dP/dphi checked by finite differences | `ex_legendre [lat_deg=45] [nmax=200]` (nmax 2..1800, no model needed) | under 1 s at the defaults | `gsh_legendre` |
| `ex_functionals_table.cpp` | Table of points read with `read_table`, then T, zeta_ell, Dg_sa, dg and V_xx, V_zz, V_xz with `compute_functionals`, written to CSV | `ex_functionals_table <model.gfc> examples/data/points_demo.csv <out.csv> [nmax]` | dominated by reading the model (tens of seconds for degree 2190) | `gsh_expansion` (functionals), `gsh_utils` (table_io) |
| `ex_grid.cpp` | Gauss-Legendre and regional grids, `synth_grid` timing, V = GM/r S_0 and T of the disturbing model, regional CSV | `ex_grid <model.gfc> <regional.csv> [nmax=90]` | under a few seconds plus model reading | `gsh_expansion` (grid, shs) |
| `ex_errors.cpp` | Commission error from formal sigmas, Kaula and Tscherning-Rapp omission errors, and the error report | `ex_errors <model.gfc> [nmax=2190]` | dominated by reading the model (read twice) | `gsh_expansion` (errors, error_report) |
| `ex_ellipsoid_tides.cpp` | GRS80 normal gravity (Somigliana against spherical harmonics), coordinate conversions, and tide-system shifts of C20, N, h and gravity | `ex_ellipsoid_tides <model.gfc>` (reads degree 2 only) | dominated by reading the file | `gsh_utils` (ellipsoid, coordinates, tides) |
| `ex_model_inspect.cpp` | Header of a `.gfc` model, degree variances, and C20, C21, S21 | `ex_model_inspect <model.gfc> [nmax]` | dominated by reading the model | `gsh_utils` (model), `gsh_expansion` (errors) |

The new examples (all except `ex_puntual`) print a usage line for bad arguments and return exit code 2. A failed
computation prints the library's error message and returns 1.

## Input data

- `data/points_demo.csv` (in this folder, `examples/data/`): five points with the header
  `id,lat_deg,lon_deg,h_m`. Latitudes and longitudes are geodetic on GRS80, heights are ellipsoidal in metres.
  Ids are numeric. The file must not contain rows at the exact poles when `dg` is requested.

## Notes

- Angles in the examples are degrees unless the name says otherwise. Latitudes given to the spherical-harmonic
  routines (for example the grid rows) are geocentric (spherical).
- Units are SI in the CSV files. The screen output of `ex_functionals_table` shows mGal and Eotvos.
- Timings are not in these files; the run-time column above is an estimate, not a measurement.
