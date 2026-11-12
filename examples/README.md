# Examples

Small programs that show the public API of the `gsh` library (`src/cpp/gsh_utils`, `src/cpp/gsh_legendre` and
`src/cpp/gsh_expansion`). Each example is one file in `cpp/` (C++) or `py/` (Python, same workflows through the `gshtools` package), with a header comment that says what it
shows, the command line and what the output means. They are built with the other targets of the project
(see the main README); the executables go to `build/<preset>/bin`.

Most examples read a fully normalised ICGEM `.gfc` model. The model file is not part of the repository: any
fully normalised model works, and EGM2008 (degree 2190, with formal errors in the file) was used for development.
`ex_errors` needs the formal-error columns (`sigma C`, `sigma S`), which EGM2008 has.

## Folders

- `cpp/`: the C++ examples, built with the other targets of the project (CMake). Executables go to `build/<preset>/bin`.
- `py/`: the same workflows in Python, through the `gshtools` package (pybind11 wrapper of the C++ library). Run them
  with `python`; see "Running the Python examples" below. Each Python example has the same arguments and prints the
  same quantities as its C++ counterpart (checked with EGM2008 at small nmax).
- `data/`: input files shared by the examples.

## C++ index (`cpp/`)

| Example | What it shows | How to run | Runtime (estimate, not measured) | Module |
|---|---|---|---|---|
| `cpp/ex_puntual.cpp` | Fully normalised Legendre functions (MSFC) at one latitude and degree | `ex_puntual [lat_deg=10] [nmax=2699]` | seconds | `gsh_legendre` |
| `cpp/ex_legendre.cpp` | sfc, msfc, era and compute() agreement, the sum-rule identity, and dP/dphi checked by finite differences | `ex_legendre [lat_deg=45] [nmax=200]` (nmax 2..1800, no model needed) | under 1 s at the defaults | `gsh_legendre` |
| `cpp/ex_functionals_table.cpp` | Table of points read with `read_table`, then T, zeta_ell, Dg_sa, dg and V_xx, V_zz, V_xz with `compute_functionals`, written to CSV | `ex_functionals_table <model.gfc> examples/data/points_demo.csv <out.csv> [nmax]` | dominated by reading the model (tens of seconds for degree 2190) | `gsh_expansion` (functionals), `gsh_utils` (table_io) |
| `cpp/ex_grid.cpp` | Gauss-Legendre and regional grids, `synth_grid` timing, V = GM/r S_0 and T of the disturbing model, regional CSV | `ex_grid <model.gfc> <regional.csv> [nmax=90]` | under a few seconds plus model reading | `gsh_expansion` (grid, shs) |
| `cpp/ex_errors.cpp` | Commission error from formal sigmas, Kaula and Tscherning-Rapp omission errors, and the error report | `ex_errors <model.gfc> [nmax=2190]` | dominated by reading the model (read twice) | `gsh_expansion` (errors, error_report) |
| `cpp/ex_ellipsoid_tides.cpp` | GRS80 normal gravity (Somigliana against spherical harmonics), coordinate conversions, and tide-system shifts of C20, N, h and gravity | `ex_ellipsoid_tides <model.gfc>` (reads degree 2 only) | dominated by reading the file | `gsh_utils` (ellipsoid, coordinates, tides) |
| `cpp/ex_model_inspect.cpp` | Header of a `.gfc` model, degree variances, and C20, C21, S21 | `ex_model_inspect <model.gfc> [nmax]` | dominated by reading the model | `gsh_utils` (model), `gsh_expansion` (errors) |

The new examples (all except `ex_puntual`) print a usage line for bad arguments and return exit code 2. A failed
computation prints the library's error message and returns 1.

## Python index (`py/`)

| Example | What it shows | How to run |
|---|---|---|
| `py/ex_legendre.py` | sfc, msfc, era and compute() agreement, the sum-rule identity, and dP/dphi checked by finite differences (C++: `cpp/ex_legendre.cpp`) | `python examples/py/ex_legendre.py [lat_deg=45] [nmax=200]` (nmax 2..1800, no model needed) |
| `py/ex_functionals_table.py` | Table of points read with `read_table`, then T, zeta_ell, Dg_sa, dg and V_xx, V_zz, V_xz with `compute_functionals`, written to CSV (C++: `cpp/ex_functionals_table.cpp`) | `python examples/py/ex_functionals_table.py <model.gfc> examples/data/points_demo.csv <out.csv> [nmax]` |
| `py/ex_grid.py` | Gauss-Legendre and regional grids, `synth_grid` timing, V and T on the grids, regional CSV; optional figure with `--plot` (C++: `cpp/ex_grid.cpp`) | `python examples/py/ex_grid.py <model.gfc> <regional.csv> [nmax=90] [--plot out.png]` |
| `py/ex_errors.py` | Commission error from formal sigmas, Kaula and Tscherning-Rapp omission errors, and the error report (C++: `cpp/ex_errors.cpp`) | `python examples/py/ex_errors.py <model.gfc> [nmax=2190]` |
| `py/ex_ellipsoid_tides.py` | GRS80 normal gravity (Somigliana against spherical harmonics), coordinate conversions, and tide-system shifts (C++: `cpp/ex_ellipsoid_tides.cpp`) | `python examples/py/ex_ellipsoid_tides.py <model.gfc>` (reads degree 2 only) |
| `py/ex_model_inspect.py` | Header of a `.gfc` model, degree variances, and C20, C21, S21 (C++: `cpp/ex_model_inspect.cpp`) | `python examples/py/ex_model_inspect.py <model.gfc> [nmax]` |

### Running the Python examples

The examples need the `gshtools` package (Python 3.10 or later, numpy). matplotlib is needed only for
`ex_grid.py --plot`, and it is imported only then.

- In-tree build, without installation: build the extension in `build/py` (see `src/python/README.md`) and point
  `PYTHONPATH` at it. From the repository root:
  - PowerShell: `$env:PYTHONPATH = "build\py\python"`, then `python examples\py\ex_legendre.py`
  - Git Bash: `PYTHONPATH=build/py/python python examples/py/ex_legendre.py`
- Installed package: `pip install .` with the toolchain set in `CMAKE_ARGS` (see `src/python/README.md`). After that
  the examples run without `PYTHONPATH`.

The model files are not in the repository. Pass the path of a fully normalised `.gfc` file, for example
`data/cache/modelos/EGM2008/EGM2008.gfc` (not in git). `ex_errors.py` needs the formal-error columns.

## Input data

- `data/points_demo.csv` (in this folder, `examples/data/`): five points with the header
  `id,lat_deg,lon_deg,h_m`. Latitudes and longitudes are geodetic on GRS80, heights are ellipsoidal in metres.
  Ids are numeric. The file must not contain rows at the exact poles when `dg` is requested.

## Notes

- Angles in the examples are degrees unless the name says otherwise. Latitudes given to the spherical-harmonic
  routines (for example the grid rows) are geocentric (spherical).
- Units are SI in the CSV files. The screen output of `ex_functionals_table` shows mGal and Eotvos.
- Timings are not in these files; the run-time column above is an estimate, not a measurement.
