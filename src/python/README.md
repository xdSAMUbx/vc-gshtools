# gshtools: gravity-field functionals from spherical-harmonic models

`gshtools` evaluates 45 functionals of the Earth's gravity field (potential V, disturbing potential T, gravity
vector and anomalies, geoid undulation N, height anomaly zeta, deflections of the vertical, and the second-derivative
tensor in spherical and local coordinates) from a spherical-harmonic model such as EGM2008 or XGM2019e, at scattered
points or on grids, up to degree and order 5540 and beyond. It is a thin Python layer over a C++20 library that uses
AVX2 or AVX-512 when the CPU has them and all cores through OpenMP. **The C++ library is the reference
implementation**: the wrapper converts types and calls C++; no numerical logic is written in Python.

## Install

```bash
pip install gshtools
```

The package is distributed as **precompiled wheels only** (no compiler needed): Windows x86-64 and Linux x86-64 (manylinux 2.28), Python 3.10 to 3.13.
One wheel serves Intel and AMD processors: at run time it detects the CPU and uses the AVX-512 kernel, else AVX2,
else a scalar kernel. `gshtools.active_isa_name()` tells which one is in use. The only dependency is NumPy.

## Quick start

You need a model file in ICGEM `.gfc` format (for example EGM2008 or XGM2019e, free at
<https://icgem.gfz-potsdam.de/>). The model is read once; use `nmax` to keep only the degrees you need.

```python
import numpy as np
import gshtools as gsh

print(gsh.active_isa_name())                         # "avx512", "avx2" or "scalar"

model = gsh.load_gfc("EGM2008.gfc", nmax=360)        # degrees up to 360; nmax < 0 keeps all
pts = gsh.Points.from_geodetic([45.0, 38.2], [-70.0, -105.6], [0.0, 5376.0])   # lat deg, lon deg, height m

req = gsh.FunctionalRequest()
req.outputs = [gsh.Functional.g, gsh.Functional.T, gsh.Functional.V_zz]
req.nmax = 360
res = gsh.compute_functionals(model, pts, req)
g = res.values[gsh.Functional.g]                     # numpy array, m/s^2, one value per point
print(g, res.timing_s)
```

Ask for all the functionals you need in **one** call: the library plans the minimum number of synthesis passes, and
asking for many together costs much less than asking for each one alone. `gsh.all_functionals()` lists the 45.

Geoid undulation N and height anomaly zeta also need orthometric heights (a DEM):

```python
grid = gsh.load_gravsoft_grid("Colorado_dem.grd")    # GRAVSOFT grid, heights in metres
req.outputs = [gsh.Functional.N, gsh.Functional.zeta]
req.height = gsh.height_from_grid(grid, 5.0 / 60.0)  # block-average over 5 arc-minutes
```

Other entry points: `gsh.Ellipsoid.GRS80()`, `gsh.coords` (coordinate conversions), `gsh.tides` (permanent-tide
systems), `gsh.legendre` (`make_coefs`, `sfc`, `msfc`, `era`, `compute`, `dphi`), the batched
`gsh.legendre_points(lat, nmax, orders=..., nmin=..., derivative=...)` (selected orders of the associated Legendre
functions at many points), `gsh.Grid` and `gsh.synth_grid`, `gsh.prepare` and `gsh.synth_points` (raw
sums), `gsh.commission_error`, `gsh.make_error_report`, `gsh.read_table`.

## Accuracy and speed

Results agree with an independent sequential double-precision implementation to better than 1e-12
(max|d| / max|ref|) in the benchmarks of the repository, and with GrafLab, CHarm and SHTOOLS at degree 2190. Timings
and the full benchmark reports are in the repository (`benchmark/`). Speed depends on the CPU: on a 6-core AMD
Ryzen 7 with AVX-512, V, T and dg at 509 points and degree 2190 take about 0.23 s.

## Support and citation

Issues and questions: <https://github.com/xdSAMUbx/vc-gshtools/issues>. Licensed under Apache-2.0. See
`CITATION.cff` in the repository to cite the software.

## Units

| Quantity | Unit | Where |
|---|---|---|
| Geodetic latitude and longitude, `*_deg` arguments | degrees | `Points.from_geodetic`, `Points.from_spherical`, `coords`, `tides` |
| Geocentric latitude and longitude stored in `Points` (`Points.lat`, `Points.lon`) | radians | `Points`, `Grid.lat`, `legendre` |
| Geocentric radius (`Points.r`, `Grid.r`) and ellipsoidal heights `h` | metres | everywhere |
| `GM` | m^3 s^-2 | `Model.GM`, `Ellipsoid.GM` |
| Gravitational potential `V`, `W`, `T`, `dV_phi`, `dV_lambda` | m^2 s^-2 (angular derivatives per radian) | `compute_functionals` |
| `dV_r`, `g_x`, `g_y`, `g_z`, `g`, `g_sa`, `dg_sa`, `Dg_sa`, `dg` | m s^-2 | `compute_functionals` |
| Second derivatives `V_*`, `T_*`, `W_rr` | s^-2 (1 Eotvos = 1e-9 s^-2) | `compute_functionals` |
| `xi`, `eta`, `Theta` | radians | `compute_functionals` |
| `zeta_ell`, `N`, `zeta` | metres | `compute_functionals` |
| Normal gravity (`Ellipsoid.normal_gravity`, `somigliana`) | m s^-2 | `Ellipsoid` |
| Tidal gravity shifts | mGal (1 mGal = 1e-5 m s^-2) | `tides` |
| Raw sums `S0`, `SR`, ... (`synth_points`, `synth_grid`) | dimensionless | `RawFields` |
| Fully normalised coefficients `C`, `S` | dimensionless | `Model` |

The shift functions in `tides` are defined as "target minus source" (see `tides.hpp`).

## Memory and threads

- Array inputs are copied into the C++ containers they feed (the C++ API takes owning containers). Float64
  C-contiguous arrays are not converted first; other dtypes and strided views are converted.
- Outputs are zero-copy: the C++ storage is handed to NumPy.
- `compute_functionals`, `synth_points`, `synth_grid`, `commission_error`, `make_error_report`, `load_gfc` and
  the Legendre evaluators release the GIL while C++ runs.

## Errors

C++ exceptions map to Python exceptions: `std::invalid_argument` becomes `ValueError`, `std::runtime_error`
becomes `RuntimeError`, and out-of-range indices raise `IndexError`. For example, `nmax` above the model degree,
N or zeta without a height source, and `eta` at a geocentric pole all raise `ValueError`.

## Tests

```bash
python -m pytest tests/py
```

The tests use `data/cache/modelos/EGM2008/EGM2008.gfc` truncated to degree 60 and are skipped when the file is
missing.

## Build from source (developers)

Requirements: Windows with clang-cl (or MSVC), CMake 3.24 or later, Ninja, Python 3.10 or later, and the vcpkg
package `eigen3`.

```bash
# pip needs the toolchain through CMAKE_ARGS (the vcpkg toolchain and the compiler).
# PowerShell: $env:CMAKE_ARGS = "-DCMAKE_CXX_COMPILER=clang-cl.exe -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows"
pip install .
```

The build always uses Release with optimisation (`/O2`). The root `CMakeLists.txt` option `GSH_BUILD_PYTHON` is
OFF by default; `pyproject.toml` turns it on together with `GSH_BUILD_TOOLS=OFF` and `BUILD_TESTING=OFF`, so only the
library and the module are built.

In-tree build (no installation), for example in `build/py`:

```bash
cmake -S . -B build/py -G Ninja -DCMAKE_BUILD_TYPE=Release -DGSH_BUILD_PYTHON=ON -DGSH_BUILD_TOOLS=OFF -DBUILD_TESTING=OFF ...
cmake --build build/py --target _core
set PYTHONPATH=build\py\python
```

Windows DLLs: the extension needs the OpenMP and MSVC runtimes. The build copies them next to the
extension (`src/python/cmake/bundle_runtime.cmake`, run after the build and at install time). `gshtools/__init__.py`
then adds the package directory with `os.add_dll_directory`, so `import gshtools` works without PATH changes.
