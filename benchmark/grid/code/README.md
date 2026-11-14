# Grid synthesis benchmark: external-engine drivers

Drivers that synthesise the gravitational potential V and the disturbing potential T on latitude-longitude grids at
degree and order 2190 with the external engines: CHarm (C and the pyharm wrapper), SHTOOLS (Fortran and the
pyshtools wrapper) and GrafLab (MATLAB). The gsh driver is not in this folder.

**Status: written, not run.** None of these files was compiled, executed or checked against an engine's output.
Read the "Not verified" sections before using any number from them.

## Grids, quantities and output (common contract)

| case | grid | rows x columns |
|---|---|---|
| `gl` | Gauss-Legendre, roots of P_{N+1}, north to south | (N+1) x (2N+2) |
| `dh1` | Driscoll-Healy, rows every 180/(2(N+1)) deg from the north pole, equator row included | 2(N+1) x 2(N+1) |
| `dh2` | same rows as dh1 | 2(N+1) x 4(N+1) |
| `regional` | lat 45 -> 35 deg, lon -110 -> -100 deg, step 0.05 deg | 201 x 201 |

* Latitudes are spherical (geocentric) at the constant radius R of the model header. h is not used.
* Global longitudes: lon_j = j * 360/nlon, j = 0..nlon-1. Regional longitudes are the explicit range above.
* V = GM/r sum_n (R/r)^n sum_m (...) with r = R, degree 0 included (C00 of the file), no normal potential
  subtracted. Units m^2/s^2.
* T = V - U, U = GRS80 normal potential, closed form (`shs_common.normal_potential_closed`, same as
  `benchmark/shs/code/c/charm_shs_bench.c` and `benchmark/shs/code/fortran/shtools_shs_bench.f90`). Because U
  depends on latitude only at r = R, it is evaluated once per latitude row. The closed form already contains the
  GM term of the ellipsoid, so no separate delta C00 term is applied to the potential.
* Model: `benchmark/shs/data/XGM2019e_2190_tide_free.gfc` (decompress the `.xz` first). Every driver reads
  GM and R from its header and checks them against the engine's own values.

Outputs, per run (never a whole grid):

* `<out-dir>/<engine>_<case>_rows.csv`: header `row,lat_deg,lon_deg,V_m2s2,T_m2s2`, values `%.17g`.
  `row` is 0-based, 0 = northernmost. Exception: GrafLab writes `graflab_regional_<recurrence>_rows.csv`.
* `<out-dir>/timings.csv`: appended, one line per repetition:
  `engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb`.

Row selection:

* Global cases: the row whose latitude is closest to 0, 30, 60, 85 and 89.5 deg north, and to -30, -60, -85 and
  -89.5 deg, plus rows 0 and nlat-1 (about 10 rows, each with all nlon columns).
* Regional default: rows 0, 50, 100, 150 and 200 (45, 42.5, 40, 37.5 and 35 deg). This default is my choice, the
  brief did not fix it.
* `--rows i,j,k` overrides the default in every driver. Indices are 0-based.

Timing:

* The timed block is the synthesis only: CHarm `charm_shs_point`, pyharm `shs.point`, SHTOOLS `MakeGridDH`,
  pyshtools `MakeGridDH`, or GrafLab's grid call. The GM/R scaling, U and T are computed after the timed runs.
* `peak_mem_mb` is the process peak RSS (`getrusage` `ru_maxrss` in C and Python on Linux/WSL, VmHWM in the
  Fortran driver). It is cumulative over the process, not per repetition. -1 means not measured.
* `recurrence` labels: CHarm `fft` or `chebyshev` (a derived expectation, see below); SHTOOLS, pyshtools `native`;
  GrafLab `sfc`, `msfc`, `era`.

## Driver 1: CHarm (`c/charm_grid_bench.c`, `c/build_charm_grid_driver.sh`)

What it does: builds the point structure of the case with CHarm's constructors (`charm_crd_point_gl(nmax, R)`,
`charm_crd_point_dh1`, `charm_crd_point_dh2`) or, for regional, `charm_crd_point_init(CHARM_CRD_POINT_GRID, 201,
201, lat, lon, r)`, then runs `charm_shs_point()` on the whole grid, `--repeat` times. The node coordinates are read
back from the CHarm structure. The driver checks that row 0 is the northernmost row.

Build (WSL, after CHarm is built in `build/external/charm`, see `benchmark/shs/code/c/README.md`):

```
bash benchmark/grid/code/c/build_charm_grid_driver.sh
mkdir -p build/grid_out
build/external/charm-grid/charm_grid_bench --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \
    --nmax 2190 --case dh1 --threads 8 --repeat 3 --out-dir build/grid_out
```

Limits:

* The `recurrence` label (`fft` or `chebyshev`) is derived from the conditions that `charm/charm_shs.h` documents
  for the FFT path. It is not read back from CHarm, so it is an expectation, not a measurement.
* The FFT path needs CHarm built with FFTW. The build script links `-lfftw3_omp -lfftw3`; this is assumed to match
  the CHarm build on disk and was not checked.
* The DH1 equator row latitude is computed by CHarm as pi/2 - (pi/(2L)) * L, which may be a few 1e-17 rad rather
  than exactly 0. The printed `lat_deg` is the value CHarm used.
* Peak memory is the process RSS, which includes the model and the grid arrays. At nmax 2190 one grid of doubles
  is about 77 MB (gl: 2191 x 4382), 154 MB (dh1: 4382 x 4382) and 307 MB (dh2: 4382 x 8764).
* Exit codes: 0 ok, 1 CHarm or runtime error, 2 bad arguments or model mismatch.

## Driver 2: SHTOOLS Fortran (`fortran/shtools_grid_bench.f90`, `fortran/build_shtools_grid_driver.sh`)

What it does: `MakeGridDH(grid, n, cilm, nmax, norm=1, sampling=1|2, csphase=1, exitstatus=...)` for dh1 and dh2,
`--repeat` times. Then V = (GM/R) * grid and T per row.

Not run, by design (exit codes):

* `gl` exits **4**. SHTOOLS `MakeGridGLQ` places 2*lmax+1 longitudes (spacing 360/(2*lmax+1); 2*lmax+2 columns with
  `extend=1`), not the contract's 2N+2 longitudes at spacing 360/(2N+2). Its output is therefore a different grid.
  The brief asked for MakeGridGLQ on this case; I did not implement it for that reason. Possible alternatives:
  point-by-point `MakeGridPoint` (too slow at 9.6 million nodes) or interpolation (changes the values). Open
  question.
* `regional` exits **3**. SHTOOLS has no local grid synthesis.

Build (WSL):

```
bash benchmark/grid/code/fortran/build_shtools_grid_driver.sh
build/external/shtools-grid-wsl/shtools_grid_bench --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \
    --nmax 2190 --case dh2 --threads 1 --repeat 3 --out-dir build/grid_out
```

Limits:

* `MakeGridDH.f95` has no OpenMP parallel region (its only OpenMP directive is `threadprivate`), so the timed
  synthesis runs on one thread. The `--threads` value is informational and is written to the timings as given.
* The header comment of `MakeGridDH.f95` says N = 2*(LMAX+1), but the code sets n = 2*lmax+2. The driver follows
  the code; the row and column counts match the contract.
* The driver reads the `.gfc` itself (the same code as `shtools_shs_bench.f90`) and does not use pyshtools.
* The build script clones SHTOOLS 4.14.1 into `build/external/shtools-wsl`. It does not reuse the Windows/MSYS2
  tree in `build/external/shtools`.

## Driver 3: pyharm (`python/run_pyharm_grid.py`)

What it does: `pyharm.crd.PointGridGL/DH1/DH2(nmax, r=R)` for the global cases, `pyharm.crd.PointGrid.from_arrays`
for regional; `pyharm.shs.point(pnt, shcs, nmax)` returns V with shape (nlat, nlon). Model from
`pyharm.shc.Shc.from_file('gfc', ...)`. Reuses `gfc_reader.py` and `shs_common.py` from `benchmark/shs/code/python`
by path.

Run (venv with `benchmark/shs/code/python/requirements.txt`, Python 3.12):

```
python benchmark/grid/code/python/run_pyharm_grid.py --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \
    --nmax 2190 --case dh1 --threads 8 --repeat 3 --out-dir build/grid_out
```

Limits:

* `OMP_NUM_THREADS` is set before `import pyharm`, as in `run_pyharm.py`.
* Memory: T is computed as a full (nlat, nlon) array; dh2 at nmax 2190 is about 38 million doubles per array.
* The `recurrence` label uses the same derived rule as the C driver.
* Peak memory is -1 on Windows, where the `resource` module is missing.

## Driver 4: pyshtools (`python/run_pyshtools_grid.py`)

What it does: `pyshtools.backends.backend_module('shtools').MakeGridDH(cilm, lmax=nmax, norm=1, sampling=1|2,
csphase=1)` with the coefficients from `SHGravCoeffs.from_file(...).to_array(normalization='4pi', csphase=1)`.

Important: pyshtools 4.14 prefers the **ducc0** backend when it is installed. The runner calls the `shtools`
module directly, so the timed library is the Fortran SHTOOLS one. It prints the preferred backend in its last line.

Not run, by design (exit codes): `regional` exits **3**; `gl` exits **4** (same reason as the Fortran driver).

Run:

```
python benchmark/grid/code/python/run_pyshtools_grid.py --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \
    --nmax 2190 --case dh2 --threads 1 --repeat 3 --out-dir build/grid_out
```

Limits:

* The `threads` value is informational (`OMP_NUM_THREADS`), like the Fortran driver.
* The keyword names of the compiled `shtools` backend `MakeGridDH` were taken from the ducc0 wrapper signature
  (same API). The compiled module was not inspected.

## Driver 5: GrafLab (`matlab/run_graflab_grid.m`, see `matlab/README.md`)

What it does: GrafLab 2.1.4 grid mode (`point_type = 0`), `coord = 1` (spherical latitude, r = R + h with h = 0),
`Functional = [11; 5]`, `fnALFs` 1, 2 or 3 (sfc, msfc, era). Regional case only. The global cases are refused with
an error. The rows file is named `graflab_regional_<recurrence>_rows.csv`, one per mode.

Limits:

* GrafLab is not in this repository. Its folder is an argument.
* `peak_mem_mb` is -1. `threads` is MATLAB's `maxNumCompThreads`.
* Its T is GrafLab's own functional 5, not the closed-form U of the other engines (see `matlab/README.md`).

## Not verified

Everything below was not checked by compiling, running or observing an engine's output.

* No file in this folder was compiled, linked, run or checked against output. This includes the C driver, the
  Fortran driver, both Python runners and the MATLAB script.
* CHarm: the compiled `libcharm.a` and its FFTW configuration (which path CHarm takes at runtime); the timing of
  the FFT path on the grids; the pyharm 0.4.11 wheel behaviour for `PointGrid.from_arrays` with regional arrays
  (the class constructor was read from `wrap/pyharm/crd.py`, not executed).
* SHTOOLS: `MakeGridDH` results against the contract (values, sign and ordering of rows and columns); the
  `libSHTOOLS-mp.a` build on WSL; the compiled `_SHTOOLS` keyword names used by the pyshtools shtools backend.
* GrafLab: the grid output layout (columns `lat, lon, V, T`, one row per node) was read from the source, not
  observed; the scalar-step grid generation was not read in full; the GrafLab T was not compared with the closed form.
* Timing and memory are not measured anywhere in this folder.

## Sources checked

Facts below were read in files on disk; anything not listed was not read.

CHarm 0.4.11 (`build/external/charm`):

* `charm/charm_crd.h`: `charm_crd_point_malloc/calloc/init`, `charm_crd_point_free`, `charm_crd_point_gl(nmax, r)`,
  `charm_crd_point_dh1(nmax, r)`, `charm_crd_point_dh2(nmax, r)`, the `*_shape` functions, and the point-type docs.
* `charm/charm_shs.h`: `charm_shs_point(pnt, shcs, nmax, f, err)`, `f[i*nlon + j]` for grids; the FFT conditions
  (nlon > 1, (nlon-1)/2 >= nmax, constant dlon, lon[0] = 0, lon[last] + dlon = 2 pi); the formula
  mu/R sum (R/r)^{n+1} ...
* `charm/charm_shc.h`: `charm_shc_calloc(nmax, mu, r)`, `charm_shc_read_gfc(path, nmax, epoch, shcs, err)`,
  `charm_shc_free`.
* `src/crd/crd_point_quad_l.c`: L = nmax + 1. `src/crd/crd_point_{gl,dh1,dh2}_shape.c`: GL (L, 2L), DH1 (2L, 2L),
  DH2 (2L, 4L).
* `src/crd/crd_point_gl_chunk.c`: GL latitudes north first (index 0 is the largest root), longitudes i*pi/L.
* `src/crd/crd_point_dh_chunk.c` and `crd_point_dh_lat_w_chunk.c`: DH latitude pi/2 - i*pi/(2L), north first, no
  south pole; longitudes i*pi/L (DH1) and i*pi/(2L) (DH2).
* `wrap/pyharm/crd.py`: `PointGridGL(nmax, r)` (line ~752), `PointGridDH1(nmax, r)` (~831), `PointGridDH2`
  (~869), `PointGrid.from_arrays(lat, lon, r)` (~686), attributes `lat`, `lon`, `r`.
* `wrap/pyharm/shs.py`: `point(pnt, shcs, nmax)`, output shape (lat.size, lon.size) for grids.
* `tools/compare_charm/charm_driver.c` (repository): uses `charm_crd_point_gl(nmax, shcs->r)`.
* `benchmark/shs/code/c/charm_shs_bench.c` (repository): includes, `charm_shc_calloc(nmax, 1.0, 1.0)`,
  `charm_shc_read_gfc(...)`, `shcs->mu`, `shcs->r`, and the normal potential.

SHTOOLS 4.14.1 (`build/external/shtools`):

* `src/SHTOOLS.f95` lines 170-180: interface `MakeGridDH(griddh, n, cilm, lmax, norm, sampling, csphase, lmax_calc,
  extend, exitstatus)`.
* `src/MakeGridDH.f95`: n = 2*lmax+2 (line ~109), nlong = n (sampling 1) or 2n (sampling 2); latitude from 90 down
  in steps of 180/n (header description); the header notes say the latitude is geocentric; the only OpenMP
  directive is `threadprivate` (line ~105).
* `src/MakeGridGLQ.f95`: nlat = lmax+1, nlong = 2*lmax+1 (lines ~108-109), +1 with `extend`.
* `src/PreGLQ.f95`: Gauss-Legendre nodes and weights (not used).
* `pyshtools/backends/backends.py`: default backend `ducc` when available (lines ~26-31),
  `backend_module(backend)` (line ~58).
* `pyshtools/backends/shtools.py` line 30: `MakeGridDH` from the compiled `_SHTOOLS` module.
* `pyshtools/expand/__init__.py`: the expand functions are injected from the preferred backend at import.
* `pyshtools/__init__.py`: `from . import backends`.

GrafLab 2.1.4 (`GrafLab.m` in the scratchpad folder, read only):

* Function header (lines ~260-345): `coord` 0 ellipsoidal, 1 spherical; `point_type` 0 grid; `lat_min`, `lat_step`,
  `lat_max`, `lon_min`, `lon_step`, `lon_max`, `h`; output description.
* Grid branch (around line 2346): for `coord == 1`, r = R + h and the latitude is spherical (`fiG = fi`).
* Export (around lines 8926-8990): for `coord == 1` the exported latitude is the spherical one; `output = [fi; lambda; P]`,
  transposed when `nargout == 1`; returned via `varargout{1}` (line ~9590).
* Normal field for T: a zonal series of the ellipsoid, n = 0..10 (around line 2370).

Repository (read): `benchmark/shs/code/c/charm_shs_bench.c`, `benchmark/shs/code/c/build_charm_driver.sh`,
`benchmark/shs/code/fortran/shtools_shs_bench.f90`, `benchmark/shs/code/fortran/build_shtools_driver.sh`,
`benchmark/shs/code/python/{run_pyharm.py,run_pyshtools.py,shs_common.py,gfc_reader.py,requirements.txt}`,
`benchmark/shs/code/matlab/{graflab_evaluate.m,read_gfc_header.m}`, `benchmark/shs/CONTRACT.md`.
