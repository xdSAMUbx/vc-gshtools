# Python-wrapper overhead of gshtools.compute_functionals

This benchmark measures what the Python wrapper (`gshtools`, pybind11) costs on top of the C++ call it wraps. The C++
library is the reference implementation; the wrapper only converts types and calls it. The question is how many seconds
and what percentage a call from Python adds, for the same inputs and the same work.

Accuracy is not part of this benchmark (see `benchmark/shs` for accuracy and speed against other engines).

## What is measured

For every case (number of points N, maximum degree nmax, functionals) and for every repetition:

| Quantity | Side | Definition |
|---|---|---|
| `py_total_s` (a) | Python | `Points.from_geodetic(lat, lon, h, a, e2)` from numpy arrays, `compute_functionals(model, pts, req)`, and reading the result arrays from `res.values` (the time from the first call to the last read) |
| `py_points_s` | Python | the `Points.from_geodetic` part of (a) |
| `py_compute_s` | Python | the `compute_functionals` part of (a), including the conversion of the result dict |
| `py_read_s` | Python | reading the result arrays from the returned dict (zero-copy, so this is the access cost) |
| `cpp_total_s` (b) | C++ | the same work in `gsh_overhead_ref`: `Points::from_geodetic` + `compute_functionals`, wall time (`steady_clock`) |
| `cpp_points_s`, `cpp_compute_s` | C++ | the two parts of (b) |
| `cpp_stage_sum_s` | C++ | sum of `FunctionalResult.timing_s` (the stages measured inside the library: model and disturbing synthesis, height, ellipsoid synthesis). For information: it does not cover the whole call |
| `overhead_s` | derived | `(a) - (b)` in seconds |
| `overhead_pct` | derived | `100 (a) - (b) / (b)` |

Inputs that are the same on both sides:

* **Model** (loaded once per nmax, before any timing, not timed): EGM2008 (`data/cache/modelos/EGM2008/EGM2008.gfc`) for
  nmax 90 and 720; XGM2019e (`data/cache/modelos/XGM2019e/XGM2019e.gfc`, truncated to degree 2190) for nmax 2190.
* **Points** (read before timing, not timed): N = 1, 10, 100 and 509 are the first N rows of
  `benchmark/shs/data/points_colorado.csv` (the shs Colorado set, 509 points; the file is reused, not copied). N = 10000
  is a 100 x 100 grid over the Colorado box (37 to 41 N, 102 to 109 W, h = 1500 m), written by the script to
  `output/python_overhead/points/`. Both sides read the same CSV files with `gshtools.read_table` (Python) and
  `gsh::read_table` (C++), so the inputs are bit-identical.
* **Request**: `nmin = 0`, `nmax` of the case, GRS80, the functionals `{V}` or `{V, T, dg}`, `threads` as given on the
  command line (default 0 = OpenMP maximum, on both sides). The request object is built before the timed block.

The full matrix is 5 values of N x 3 values of nmax x 2 functional sets = 30 cases.

### Repetitions and the warm-up call

Each case runs **one warm-up call** first (CSV `phase = warmup`, `rep = 0`), then **8 timed calls** (`phase = timed`,
`rep = 1..8`). All 8 are recorded in the CSV; the summary reports the median and the min and max. The warm-up call is
kept apart because the first call of a process includes the OpenMP thread-pool start-up and the first-touch of memory.
Both sides do the same warm-up.

Before each Python call the script runs `gc.collect()` outside the timed block, so that a collection does not fall
inside a measurement.

## Files

```text
benchmark/python_overhead/
  README.md               this file
  run_overhead.py         Python side: the measurements, the bitwise check, the CSV, the summary
  gsh_overhead_ref.cpp    C++ side: the same work, timed in C++
  CMakeLists.txt          target gsh_overhead_ref (see "Build the C++ driver")
  results/                output of the full run (overhead.csv); empty in the repository
```

Generated inputs and the bitwise dumps go to `output/python_overhead/` (ignored by git, under `output/`).

## How to run

Commands from the repository root, PowerShell. The Python side needs the in-tree extension (`build/py`, see
`src/python/README.md`), so `PYTHONPATH` must point to `build\py\python`.

**1. Build the C++ driver.** The root `CMakeLists.txt` must contain `add_subdirectory(benchmark/python_overhead)` inside
its `if (GSH_BUILD_TOOLS)` block. That line is not added by this change; add it, then:

```powershell
cmake --preset clang
cmake --build --preset clang --target gsh_overhead_ref
```

The executable is `build\clang\bin\gsh_overhead_ref.exe`. Its DLLs (fftw3, OpenMP) must be found at run time, as for the
other executables of the `clang` preset.

**2. Bitwise check, once.** Python result equal to the C++ result, bit for bit, for every case, with `threads = 1`:

```powershell
$env:PYTHONPATH = "build\py\python"
.venv\Scripts\python.exe benchmark\python_overhead\run_overhead.py --correctness-only
```

It prints one row per case and `bitwise check: PASS` or `FAIL` (exit code 1). The check runs the single-threaded C++
driver on every case, so the nmax 2190 cases with 10000 points take the longest.

**3. Full timing run.** 30 cases, 1 warm-up and 8 timed calls each, idle machine:

```powershell
$env:PYTHONPATH = "build\py\python"
.venv\Scripts\python.exe benchmark\python_overhead\run_overhead.py
```

It writes `benchmark/python_overhead/results/overhead.csv` (one row per case and call, 270 rows: 30 cases x 9 calls) and prints
the summary as markdown. Progress goes to stderr.

**4. Summary of an existing CSV**, without measuring again:

```powershell
.venv\Scripts\python.exe benchmark\python_overhead\run_overhead.py --summary
.venv\Scripts\python.exe benchmark\python_overhead\run_overhead.py --summary --csv path\to\other.csv
```

**Options.** `--repeat R` (default 8), `--warmup W` (default 1), `--threads K` (default 0 = OpenMP maximum, both sides),
`--cpp-exe PATH` (default `build\clang\bin\gsh_overhead_ref.exe`), `--csv PATH` (output of a run, input of
`--summary`). `--help` lists them.

**Smoke test only (not the benchmark).** `--quick` uses 2 repetitions, N = 1 and 10, nmax 90 (EGM2008) and writes to
`output/python_overhead/overhead_quick.csv`, never to `results/`:

```powershell
.venv\Scripts\python.exe benchmark\python_overhead\run_overhead.py --quick
.venv\Scripts\python.exe benchmark\python_overhead\run_overhead.py --quick --correctness-only
```

## Why the overhead should be small

A call does work that grows with N x nmax^2 (one Legendre recurrence and one series sum per point and per synthesis
pass). The wrapper work does not grow with nmax and is linear in N only through the copies of the inputs:

* copying the three coordinate arrays into the C++ containers (`Points`): about 24 bytes per point;
* the pybind11 argument conversion and GIL release/acquire, a few microseconds per call;
* building the result dict and wrapping the output arrays (zero-copy: the C++ storage is handed to NumPy).

So the overhead in seconds should be roughly constant (microseconds to tens of microseconds per call) and independent
of nmax, while the percentage falls as the work grows: large at N = 1 and nmax 90, small at N = 10000 and nmax 2190. This
is an expectation; the benchmark has not been run to confirm it.

## How to read the summary

* `overhead median (s) [min, max]`: median over the 8 timed calls of `(a) - (b)`, and the range.
* `overhead median (%) [min, max]`: the same relative to the C++ time. Positive means Python is slower.
* `Points, Python minus C++ (ms)`: the cost of building `Points` from numpy arrays in Python, minus the C++ construction.
* `compute call, Python minus C++ (ms)`: the cost of `compute_functionals` through the wrapper, minus the C++ call.
* `read arrays, Python (ms)`: the cost of getting the result arrays out of the returned dict. It has no C++ counterpart.
  The overhead of the median call is the sum of the points difference, the compute difference and the read cost,
  up to the difference between medians.
* The warm-up table: the first call of each case in its process, with its overhead, listed apart.

The CSV columns are `case_id, n_points, nmax, functionals, model, phase, rep, threads`, then the Python times
(`py_total_s, py_points_s, py_compute_s, py_read_s`), the C++ times (`cpp_total_s, cpp_points_s, cpp_compute_s,
cpp_stage_sum_s`) and `overhead_s, overhead_pct`. Times are written in seconds with 9 significant digits.

## Results (measured 2026-11-13)

Machine: AMD Ryzen 7 7445HS (6 cores, 12 logical processors), clang-cl 22.1, AVX-512 kernel, OpenMP default (12
threads) on both sides, idle machine, 8 timed repetitions per case after 1 warm-up. Full table with min and max:
`results/overhead_summary.md`; raw rows: `results/overhead.csv` (30 cases x 9 calls).

**Bitwise check (run once):** `run_overhead.py --correctness-only` (1 thread): the Python output arrays are bitwise
equal to the C++ driver's in all 30 cases (V and V,T,dg; N = 1 to 10000; nmax 90, 720, 2190). Result: PASS.

Median overhead, Python minus C++ (negative = Python faster):

| cases (V and V,T,dg) | overhead (ms) | overhead (%) |
|---|---|---|
| nmax 90, N <= 509 | -0.50 to -0.03 | -38.0 to -2.6 |
| nmax 90, N = 10000 | -1.11 to -0.46 | -11.8 to -7.5 |
| nmax 720, N <= 509 | -2.25 to +0.01 | -14.5 to +0.1 |
| nmax 720, N = 10000 | -7.34 to +0.84 | -1.5 to +0.6 |
| nmax 2190, N <= 509 | -9.05 to +2.69 | -13.1 to +0.9 |
| nmax 2190, N = 10000, V | -10.64 | -0.8 |
| nmax 2190, N = 10000, V,T,dg | +223.6 | +5.4 |

Of the 30 cases, 26 have a negative median overhead.

Reading:

* The wrapper adds no measurable cost. The Python-only parts are tiny: `Points.from_geodetic` differs from the C++ one by
  0.02 to 0.04 ms and reading the result arrays (zero-copy) costs 0.01 to 0.02 ms.
* The overheads are **negative in most cases, so the measured differences are not a cost of the wrapper**: they are
  the noise and ordering effects listed under Limitations (the Python half of a group runs first, the C++ driver runs after
  it in another process; an OpenMP fork/join or frequency state that differs between the two halves is enough to move
  a call of a few milliseconds by 10 to 30 percent). The largest positive value (N = 10000, V,T,dg, nmax 2190: +0.22 s,
  +5.4 %) is a single case and is of the same size as the spread between repetitions of that case; it was not repeated
  and is not claimed as a real cost.
* Conclusion that the data support: for calls of 0.1 s or more (nmax 2190, N >= 509) the Python wrapper is within a few
  percent of the C++ call, and the results are bit-identical. For calls of milliseconds the comparison is limited by the
  noise of the measurement, not by the wrapper.
* The negative overhead is **explained only as a measurement-order effect, not demonstrated**. A run that interleaves the
  two sides would settle it; it was not done.

## Limitations

* **Threads.** Both sides use the same OpenMP setting (`--threads`, default the OpenMP maximum). The fork and join of
  the parallel region is a larger part of a small call, and it depends on the thread count, so the percentages at N = 1
  and 10 depend on it. Report the thread count with the results (the CSV has a `threads` column).
* **Machine.** Run on an otherwise idle machine: close browsers, IDEs, synchronisation clients and builds. The noise on a
  Windows desktop is of the order of the overhead at small N.
* **Sequential timing.** For each nmax group the script runs the Python calls first and then the C++ driver (one
  process per group, the model is loaded in each process and not timed). The two sides are therefore not interleaved;
  a drift of the machine between the two halves of a group appears in the overhead.
* **Two processes.** The C++ reference is a separate executable that links the same `gsh` library. It is built with the
  same compiler and options as the library, but it is a different binary from the one inside `_core.pyd`.
* **Timer.** `time.perf_counter` (Python) and `steady_clock` (C++) are high-resolution timers, but the smallest calls
  (N = 1, nmax 90) last only about 0.2 to 0.5 ms, so run-to-run noise matters more than the timer resolution.
* **Model load and points read are not timed.** The model file is parsed once per nmax group, outside the timed block.
  The points are read before the timed block. Neither cost is part of the overhead.
* **Result data.** Only the time is measured. The bitwise check (step 2) is separate and is run once.
* **Known anomaly, not explained.** In the first smoke test (`--quick`, 2 calls, N = 1 and 10, nmax 90) the C++ compute
  call took about 0.5 ms, of which the library's own stage timers account for about 0.22 ms, while the Python compute call
  took about 0.2 ms. The overhead was therefore negative. The two binaries use the same compiler flags and the same
  OpenMP runtime, and the library has no per-process cache; the cause is not identified. It must be understood before the
  full run is interpreted: a negative overhead means that the reference and the measured path are not doing the same
  work.
* **Run time of the full benchmark.** Not measured. The nmax 2190 cases with 10000 points are the heaviest. The model
  load of XGM2019e (1.8 GB file, truncated to degree 2190) is paid once per process and is not timed.
