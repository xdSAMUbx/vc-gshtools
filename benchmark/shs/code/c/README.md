# CHarm engine of the SHS benchmark (native Windows, MSYS2 UCRT64)

Driver `charm_shs_bench.c` runs CHarm (charm 0.4.11) on the spherical-harmonic synthesis benchmark of
`../../results/csv.md` (definitions), engine `charm`, recurrence `native`. It writes the result CSV
`id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal`, its companion `.md`, and one row per repetition in
`timings.csv`. Definitions are those of `../python/run_pyharm.py` and `../python/shs_common.py`, ported to C.

Nothing from CHarm or FFTW is redistributed here. CHarm is built locally from its public source; FFTW is an
MSYS2 package used as an external build dependency.

## Files

| File | Content |
|---|---|
| `charm_shs_bench.c` | the driver (C99/gnu17, OpenMP, CHarm public API) |
| `build_charm_driver.sh` | builds `libcharm.dll` (from CHarm's static archive) and `charm_shs_bench.exe` |
| `README.md` | this file |

Build products go to `build/external/` (git-ignored by the `**/build/` rule): the CHarm source and build tree in
`build/external/charm`, the DLL and the driver in `build/external/charm-driver`.

## Versions

| Component | Version | Source |
|---|---|---|
| CHarm | v0.4.11, commit `1ead5828b527b084b077daa06db811f9993414d3` | `git clone --branch v0.4.11 https://github.com/blazej-bucha/charm.git` (anonymous HTTPS) |
| pyharm (reference runner) | 0.4.11 (its wrapper of CHarm) | `.venv` of the project |
| GCC | 16.2.0 (`mingw-w64-ucrt-x86_64-gcc 16.2.0-6`) | MSYS2 UCRT64 |
| FFTW | 3.3.11 (`mingw-w64-ucrt-x86_64-fftw 3.3.11-2`), OpenMP build `libfftw3_omp` | MSYS2 UCRT64 |
| GNU make | 4.4.1 (`make 4.4.1-3`) | MSYS2 |
| OpenMP runtime | `libgomp-1.dll` shipped with GCC 16.2.0 | MSYS2 UCRT64 |

The CHarm version is the one of pyharm 0.4.11: pyharm is built from `wrap/pyharm` of the same repository, and
its version comes from the same `VERSION` file (value `0.4.11` in the v0.4.11 tag). The build uses the tag, not the
default branch (HEAD is a later commit, `1ead5828` is the tag's commit as checked with `git ls-remote`).

Other MSYS2 packages present (`cmake`, `ninja`, `gcc-fortran`, `openblas`) are not used by this build.

## Build

Run in the MSYS2 UCRT64 shell (`C:\msys64\ucrt64.exe`), in this order:

```sh
# 0) packages (only if missing): gcc, fftw, make, tar, git
pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-fftw make tar git

# 1) CHarm source (anonymous HTTPS, no credentials)
cd <repo>
GIT_TERMINAL_PROMPT=0 git clone --depth 1 --branch v0.4.11 https://github.com/blazej-bucha/charm.git build/external/charm

# 2) configure and build CHarm (double precision, OpenMP, -O3 -march=native)
cd build/external/charm
./configure --enable-double-precision --enable-openmp CFLAGS="-O3 -march=native" LDFLAGS="-fopenmp"
make -j12
# "make install" is not needed. On Windows it fails because the source tree has both INSTALL and install
# (case-insensitive file system). libtool prints "Could not determine host file name corresponding to '.libs/'"
# during the build; it is harmless here, the static archive .libs/libcharm.a is complete (220 members).

# 3) driver: libcharm.dll + charm_shs_bench.exe
cd <repo>
bash benchmark/shs/code/c/build_charm_driver.sh
```

`build_charm_driver.sh` uses FFTW's OpenMP library, `-O3 -march=native -fopenmp`, and records these flags in the
notes of every result. Compiler warnings: none with the flags used.

The driver needs `libcharm.dll`, `libfftw3_omp-3.dll`, `libfftw3-3.dll` and `libgomp-1.dll` at run time. Put
`build/external/charm-driver` and `/ucrt64/bin` on `PATH` before running it (from the MSYS2 shell).

## Run

```sh
export PATH="<repo>/build/external/charm-driver:/ucrt64/bin:$PATH"
charm_shs_bench.exe --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \
  --points benchmark/shs/data/points_colorado.csv --nmax 2190 \
  --out benchmark/shs/results/charm_native_colorado.csv --timings benchmark/shs/results/timings.csv \
  --set colorado --threads 12 --repeat 3
```

Options: `--gfc`, `--points`, `--nmax`, `--out`, `--timings`, `--set colorado|latsweep` (required), `--threads K`
(>= 1), `--repeat R` (>= 1), `--parallel points|charm` (default `charm`, the official form; `points` is a diagnostic, not run by run_all.ps1; see "Parallel behaviour").
`--recurrence native` is accepted and is the only value. The model is the decompressed
tide-free file (`xz -dk data/XGM2019e_2190_tide_free.gfc.xz`). The driver checks `nmax <= max_degree`,
`norm = fully_normalized`, `tide_system = tide_free`, and that CHarm's GM and R equal the header values.

## What is computed

Per point, in one synthesis pass:

* `V_m2s2` from `charm_shs_point` (mu/r sum, degree 0 included, no normal potential subtracted).
* `T_m2s2 = V - U`, `U` the GRS80 level-ellipsoid potential (closed form, `normal_potential_closed`).
* `dg_mgal = (|g| - gamma_SH) * 1e5`, with `g` built from `charm_shs_point_grad1` (LNOF) plus the centrifugal
  gradient (omega = 7.292115e-5 rad/s), and `gamma_SH` the project's zonal series (`normal_gravity_gamma`).

Geodetic (GRS80) to geocentric conversion is exact and identical to the python runner. The `.md` of every result
lists the exact command and parameters.

### Gradient component order (important)

`charm_shs_point_grad1` swaps the pointers of the `f` array it receives. After the call, `f[0]` is the LNOF `x`
element (latitude: `1/r dV/dphi`), `f[1]` is `y` (longitude: `-1/(r cos phi) dV/dlambda`) and `f[2]` is `z` (radial:
`dV/dr`). The driver reads the values through the swapped pointers. Reading the buffers passed in gives the wrong
dg (about 2000 mGal off in the first smoke test). This order was verified two ways: by central finite differences of
V (step 1e-7 rad and 1 m) at a point of the smoke set, and by agreement with pyharm (whose wrapper uses the same
order, see `pyharm/shs.py`, `_point_gradn`).

## Parallel behaviour (option `--parallel`)

* **CHarm itself**: `charm_shs_point` and `charm_shs_point_grad1` for scattered points run an OpenMP parallel
  region (`#pragma omp for schedule(dynamic)`) over the spherical-harmonic orders m, inside one call, when CHarm is
  built with `--enable-openmp`. The FFT (FFTW) path is for grids and is not used here.
* **`--parallel points` (DIAGNOSTIC, not in run_all.ps1)**: the driver runs OpenMP over the points,
  `#pragma omp parallel num_threads(K)` with `schedule(dynamic,4)`. Each point has its own one-point `charm_point`
  (`charm_crd_point_init`, borrowing the driver's arrays) and its own `charm_err`. The `charm_shc` is read-only and
  shared. Inside each driver thread `omp_set_num_threads(1)` is called, so CHarm's inner region over m runs with one
  thread; `omp_set_max_active_levels(1)` is also set. The total number of threads is K.
  Without the per-thread call, a 1-thread run let CHarm's inner region take the default team (12 threads), and the
  1-thread and 4-thread outputs differed in the 10th-11th digit. With it, K = 1, 4, 12 give byte-identical output
  (smoke test, 5 points, nmax 100).
* **`--parallel charm` (default, official)**: one CHarm call over all points, as pyharm calls CHarm. CHarm's OpenMP over m uses K threads
  (`omp_set_num_threads(K)`); the driver adds no parallelism.

### Timing finding (full size, points mode)

The one full-size correctness run (nmax 2190, 509 points, 12 threads, `--parallel points`) took 23.7 s wall and
240.5 CPU-s in the synthesis, and 5.5 s to load the model. pyharm, which makes one batched call, took 0.67 s wall
and 6.7 CPU-s in its existing timed rows (`results/timings.csv`, 12 threads). The points mode is therefore far
slower in CPU time than the batched CHarm call, most likely from the per-call overhead of one point per call and
unused SIMD lanes. This is an inference, not measured. `--parallel charm` was not run at full size (the run allowance
was used); it is the mode that matches pyharm's design and should be timed on an idle machine before any comparison
is drawn.

## Validation

* Smoke, nmax 100, 5 points, against `run_pyharm.py` (current runner) with the same points and nmax: max|d| V 1.0e-7
  m^2/s^2, T 9.7e-8 m^2/s^2, dg 1.1e-9 mGal. Batch and points modes agree to 1.3e-7 (V, T) and 1.2e-9 mGal.
* Full size (nmax 2190, 509 points, points mode, 12 threads), against the existing pyharm files of this repository
  (`results/pyharm_native_colorado.csv` for T and dg, `results/pyharm_native_colorado_V.csv` for V; both written
  before the current runner): max|d| V 5.0e-7 m^2/s^2 (rms 1.5e-7), T 6.9e-7 m^2/s^2 (rms 1.8e-7),
  dg 9.4e-9 mGal (rms 2.5e-9). These are about 1e-14 relative on V and T, from summation order and the different
  build of pyharm's CHarm; they are above the 1e-7 target by a factor of about 5 for V and T.

## Known limits

* The points mode is much slower in CPU time than pyharm's batched call (see the timing finding above). The
  official runs use `--parallel charm` only (library mode, the same design as pyharm).
* CHarm's warning about polar points applies to `point_grad1`: the x and y components are inaccurate at the poles.
  No Colorado point is polar.
* `make install` is not used (see the build notes); the driver links to the static archive through the local DLL.

## Licences

* CHarm: `COPYING` (GNU GPL v2) and `LICENSE` (BSD-style, (c) 2022 B. Bucha) are both in the source tree. Which
  applies is not settled here; the `COPYING` file is the one checked. The driver source contains no CHarm code.
* FFTW: GPL. External build dependency, not redistributed.
* Model: XGM2019e, CC BY 4.0 (see `../../data/DATA.md`).
