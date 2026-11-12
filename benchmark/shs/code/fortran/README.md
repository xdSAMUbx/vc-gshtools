# SHTOOLS engine for the SHS benchmark (`shtools`, recurrence `native`)

Fortran 95 driver `shtools_shs_bench.f90` for the spherical-harmonic synthesis benchmark
(definitions: `../../results/csv.md`). It evaluates the same three functionals as `../python/run_pyshtools.py`
(V, T, dg) from the same tide-free XGM2019e model at nmax 2190 on the 509 Colorado points, using the
SHTOOLS 4.14.1 Fortran library (the library behind the pyshtools 4.14.1 that the benchmark already uses).

Only the driver and the scripts in this folder are part of the repository. SHTOOLS and FFTW are external
build dependencies: they are fetched and built under `build/external/` (git-ignored) and are not
redistributed.

## Files

| File | Content |
|---|---|
| `shtools_shs_bench.f90` | the driver (reads the `.gfc`, evaluates V, T, dg, writes the CSV, `.md` and timings rows) |
| `build_shtools_driver.sh` | fetches SHTOOLS 4.14.1, builds `libSHTOOLS-mp.a` (OpenMP) and the driver, from the MSYS2 UCRT64 shell |
| `README.md` | this file |

Build products (all under `build/external/`, git-ignored):

* `build/external/shtools/` - SHTOOLS v4.14.1 sources (commit `98c15efd26b224859e1a8f8ab0c1298368657f59`),
  `lib/libSHTOOLS-mp.a`, `include/*.mod`.
* `build/external/shtools_bench/shtools_shs_bench.exe` - the driver executable.

## Packages and versions (MSYS2, UCRT64)

Verified with `pacman -Q` on this machine:

| Package | Version | Use |
|---|---|---|
| mingw-w64-ucrt-x86_64-gcc-fortran | 16.2.0-6 | gfortran (library and driver) |
| mingw-w64-ucrt-x86_64-libgomp | 16.2.0-6 | OpenMP runtime (`-fopenmp`) |
| mingw-w64-ucrt-x86_64-fftw | 3.3.11-2 | FFTW3 (GPL), linked by the build script; only the FFT-based SHTOOLS routines need it |
| mingw-w64-ucrt-x86_64-lapack | 3.12.1-2 | `-llapack` (SHTOOLS eigen/least-squares routines) |
| mingw-w64-ucrt-x86_64-blas | 3.12.1-2 | `-lblas` |
| SHTOOLS | v4.14.1 (GitHub tag) | the Fortran library; same version line as pyshtools 4.14.1 |
| make, git | MSYS2 `make 4.4.1-3`, `git 2.56.1-1` | make is not used by the build script (see below) |

Install once from the MSYS2 UCRT64 shell (wait 30 s and retry if pacman reports a lock):

```
pacman -S --needed mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-fftw \
    mingw-w64-ucrt-x86_64-lapack mingw-w64-ucrt-x86_64-blas
```

`git` is needed only to fetch SHTOOLS (the build script clones anonymously; it never prompts for credentials).
FFTW is GPL: it is an external build dependency of this engine, nothing built here is redistributed.

## Build and run

```
# from the MSYS2 UCRT64 shell, in the repository root
bash benchmark/shs/code/fortran/build_shtools_driver.sh

# run (full size: slow; the timed run is the director's, on an idle machine)
build/external/shtools_bench/shtools_shs_bench.exe \
    --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \
    --points benchmark/shs/data/points_colorado.csv \
    --nmax 2190 --set colorado --threads K --repeat R \
    --out benchmark/shs/results/shtools_native_colorado.csv \
    --timings benchmark/shs/results/timings.csv
```

Options: `--gfc`, `--points`, `--nmax`, `--out`, `--timings`, `--set colorado|latsweep`, `--threads K`,
`--repeat R` (same as the other drivers). The `.gfc` must be the decompressed tide-free file
(`xz -dk benchmark/shs/data/XGM2019e_2190_tide_free.gfc.xz`, or `prepare_data.py`).

Outputs: `<out>` (`id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal`, 15 significant digits, ES notation),
`<out>.md` (column table, command, parameters, notes) and one row per repetition appended to `<timings>`
with `engine=shtools`, `recurrence=native`.

### Build notes

* Flags: library `gfortran -fopenmp -fPIC -O3 -march=native -std=gnu`; driver `gfortran -fopenmp -O3 -march=native -std=f2008`.
  The SHTOOLS Makefile default adds `-ffast-math`; **it is deliberately not used here** (not in the brief and it
  can change the floating-point results). The library flags are written into the driver's `.md` and timings
  notes (`library_flags=`).
* `make fortran-mp` fails in this MSYS2 setup: gfortran, when started from make, reports
  `Cannot create temporary file in C:\WINDOWS\`. The build script therefore compiles the same object list
  (read from `src/Makefile`: modules, `SRCS0`, `SRCSFFTW`, `SRCSLAPACK`) directly, with the same flags, and runs
  `ar`. The temp-directory workaround is in the script.
* The driver's `LIB_FLAGS` string is replaced by the script in the build copy, so the notes record the real flags.

## Official use (read first)

The published runs use this driver serially: `code/run_all.ps1` calls it with `--threads 1` only, and its
`parallel_mode` is `none`. The OpenMP paths described below are kept in the source for diagnosis; they are not
part of any published number. SHTOOLS itself has no parallel point routine, so nothing else is added.

## What is parallel, and where

Checked in the SHTOOLS 4.14.1 sources (`src/MakeGridPoint.f95`, `src/MakeGravGridPoint.f95`,
`src/PlmBar.f95`, `src/PlmBar_d1.f95`):

* **Library:** `MakeGridPoint` and `MakeGravGridPoint` contain no OpenMP directives. One point is evaluated
  serially inside SHTOOLS. This is the case in the pyshtools wheel too (`MakeGridPoint` there is a Python
  `np.vectorize` loop over points, so each point is one Fortran call).
* **Thread safety:** the Legendre helpers used by both routines (`PlmBar`, `PlmBar_d1`) keep their work arrays
  in `save` variables that are declared `!$OMP threadprivate` in the OpenMP build. So the `libSHTOOLS-mp.a`
  build can be called from several threads at once. This is a property of the library's OpenMP build, not
  something this driver adds. The driver only uses it with the OpenMP build; the serial build is not thread-safe.
* **Driver (OpenMP):**
  * V: the points are grouped by identical radius (as `np.unique` in the python runner). Each group needs a
    radius-scaled copy of `cilm`; the groups run in parallel (dynamic schedule), each thread owning its own
    scaled copy (`2 x 2191 x 2191` doubles = 77 MB per thread at nmax 2190).
  * dg: the gradient (`MakeGravGridPoint`, omega = 0) runs in parallel over the points.
  * Model reading, the geodetic conversion, U, gamma and dg arithmetic are serial.
* Thread count (`--threads K`) is the number of OpenMP threads of the driver and is recorded in
  `timings.csv` and in the notes.

## Radius scaling cost

The library has no radius argument: `V = GM/r * sum_n (R/r)^n S_n` is obtained by scaling the degree-n
coefficients by `(R/r)^n` before `MakeGridPoint`, one scaled copy per distinct radius. This is the same as the
python runner. For the 509 Colorado points the number of distinct radii is 509 (`distinct_radii=509` in the
notes), so 509 scaled copies of the 77 MB coefficient array are made per V pass (about 39 GB of memory writes
per pass, single-threaded equivalent). The copies are inside the timed block. **This cost has not been measured
at nmax 2190 here** (the only full-size run is the correctness run); the figure is only an estimate from the
array size.

## Model load

The reader reads the header (GM, R, max_degree, norm, tide_system) and every `gfc` line of the file with a
formatted read. At nmax 2190 the 2.4 M-line file is read in full. Load time on this machine is reported as
`load_s` in the notes; the reader was not optimised.

## Validation

Definitions, formulas and constants are ported from `../python/shs_common.py` and
`../python/run_pyshtools.py` (geodetic to geocentric, the zonal gamma series with `sr`/`sp`, the closed-form
U with the flattening term, dg with the centrifugal terms, and `V = GM/r * MakeGridPoint(scaled cilm)`).

Checks done:

* Smoke test, nmax 100, 5 points (`points_colorado.csv` first 5 rows), against `run_pyshtools.py` at nmax 100
  on the same inputs: V and dg agree to the printed digits; T differs by about 4e-8 m^2/s^2 (a few ulp at V ~ 6e7,
  consistent with libm `atan`/`sqrt` differences between the MSYS2 runtime and the numpy wheel).
* Thread consistency: nmax 100, 5 points, `--threads 1` against `--threads 4`: outputs identical byte for byte.

Full-size correctness run (nmax 2190, 509 points, 12 threads), once, scratch outputs only. Reference: `results/pyshtools_native_colorado.csv` (T, dg) and `results/pyshtools_native_colorado_V.csv` (V). The committed T/dg file has no V column.

| column | max abs diff | rms diff | max abs value |
|---|---|---|---|
| V_m2s2 | 1.04e-7 m^2/s^2 | 1.07e-8 | 6.26e7 |
| T_m2s2 | 5.51e-7 m^2/s^2 | 1.12e-7 | 2.65e2 |
| dg_mgal | 4.09e-9 mGal | 5.22e-10 | 1.13e2 |

lat, lon and h echoes are identical (0). Ids and row order are identical (509 rows).

T is above the 1e-7 target. Diagnosis (point id 433, worst T): the difference is in U, not V. At that point, a high-precision (mpmath, 50 digits) U gives 62550152.0691184; the python U is 1.7e-7 from it and the Fortran U is -3.8e-7 from it. The cause is gfortran's libm `atan` (MSYS2 UCRT runtime), which is about 1 ulp off where numpy is correctly rounded; with GM/e ~ 7.5e9 that is about 1e-8 per ulp in U, and the cancellation in `q(x)` adds to it. Checked and ruled out: FMA contraction (`-ffp-contract=off` gives identical U on a test point). Not fixed here: it would need a different libm or a higher-precision atan.

The V difference (1e-7 on 6e7, about 14 ulp) comes from the same MakeGridPoint sum ordering and libm differences; it was not isolated further.

Timing of this correctness run (not a benchmark; 12 threads, shared machine): load 23.4 s, synthesis wall 17.0 s, cpu 162.4 s (`cpu_time`, all threads). The pyshtools reference in `results/timings.csv` is 241-252 s wall for the same set, serial. The director's idle-machine timing is still to be run.
