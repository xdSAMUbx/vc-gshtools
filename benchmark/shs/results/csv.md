# CSV deliverables of the spherical-harmonic synthesis benchmark

This file explains every folder and every CSV of `results/csv/`, column by column, with units and how each file
was built. It is the reference for the numbers in `Summary.md`. Definitions of the functionals (V, T, dg), the
point sets and the model are in `../data/DATA.md` and below.

## What is compared

Four engines compute three functionals at each point of two point sets, from the same tide-free model
XGM2019e (nmax 2190, see `../data/DATA.md`):

| Engine | Language / package | What it is |
|---|---|---|
| gsh | C++ (this project), OpenMP | the project engine; recurrences `fast` (production SIMD kernel), `msfc`, `sfc`, `era` |
| GrafLab 2.1.4 | MATLAB | reference: recurrences `sfc`, `msfc`, `era` |
| pyshtools 4.14.1 | Python, SHTOOLS Fortran library | serial library, run with 1 thread |
| SHTOOLS 4.14.1 driver | Fortran 95 driver linked to the SHTOOLS library | serial loop, run with 1 thread |
| pyharm 0.4.11 | Python, CHarm 0.4.11 | one batched CHarm call (`shs.point`, `point_grad1`) |
| CHarm 0.4.11 driver | C driver linked to CHarm | one batched CHarm call (`--parallel charm`) |

Every engine is run in its official, as-distributed form. The benchmark adds no parallelism of its own to a tool
that does not have it. The CHarm `--parallel points` mode and the SHTOOLS driver's internal OpenMP exist in the
driver sources only as diagnostics; they are not part of any published number.

Point sets (`data/`):

* `colorado`: 509 points in Colorado, geodetic latitude/longitude on GRS80, ellipsoidal height in m.
* `latsweep`: 341 points on the meridian 10 deg E, latitude -85 to +85 deg in steps of 0.5 deg, height 0 m.
  It stresses the Legendre recurrences at high latitude.

Functionals (all at the point, with its own ellipsoidal height; one synthesis pass per run gives all three):

* `V_m2s2`: gravitational potential of the model, V = GM/r * sum_{n=0}^{2190} (R/r)^n sum_m (C cos m lambda + S sin m
  lambda) P_nm, degree 0 included with C00 of the file, no normal potential subtracted. Unit m^2/s^2.
* `T_m2s2`: disturbing potential T = V - U, U the GRS80 level-ellipsoid normal potential (closed form). Unit m^2/s^2.
* `dg_mgal`: gravity disturbance dg = |g| - gamma, g from the gradient of V including the centrifugal term
  (omega = 7.292115e-5 rad/s), gamma the spherical-harmonic series of the GRS80 ellipsoid. Unit mGal (1 m/s^2 = 1e5 mGal).

Legendre recurrences (gsh and GrafLab): `sfc` standard forward column; `msfc` modified forward column; `era`
extended-range arithmetic; `fast` the production SIMD kernel of gsh. Engines with one algorithm carry the
recurrence label `native`. The two CHarm modes carry the labels `native-charm` (`--parallel charm`) and
`native-points` (diagnostic, not in the published runs).

## Folder layout

```
results/
  Summary.md          written by hand from the tables of raw/summary_tables.md (not generated)
  csv.md              this file
  csv/
    graflab/    res_value.csv, res_time.csv
    charm/      res_value.csv, res_time.csv, comp_value.csv
    shtools/    res_value.csv, res_time.csv, comp_value.csv
    gshtools/   res_value.csv, res_time.csv, comp_value.csv
  raw/                 per-run outputs of the scripts (not the deliverable)
    <engine>_<rec>[_<mode>]_<set>[_t<threads>].csv    one value file per run (id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal)
    <same name>.md                                    companion of each run: columns, command, parameters
    timings.csv                                       one row per repetition of every run
    summary_tables.md                                 the tables pasted into Summary.md (written by compare.py)
    run_all.log                                       transcript of the run
img/                   figures written by compare.py
```

The folders are the groups of engines that a reader compares: `graflab` holds the reference, `gshtools` holds
gsh (compared with GrafLab), `charm` holds pyharm and CHarm, `shtools` holds pyshtools and the SHTOOLS driver.
Both point sets are in the same file of each folder; the `set` column tells them apart.

## How the files are built

1. `code/run_all.ps1` runs every engine configuration of the matrix (below) for each point set, writes one value
   file per run to `raw/`, and appends one row per repetition to `raw/timings.csv`. Repetitions: 3 for every
   configuration.
2. `code/python/consolidate.py` reads `raw/` and writes `res_value.csv` and `res_time.csv` of each folder.
   * Value files: for each (engine, recurrence, set) the run with the largest thread count is the production run
     and its values are written. The other thread counts are used only in the thread-consistency check of
     `summary_tables.md`.
   * Timings: only the rows tagged `functional=V_T_dg` (all three functionals in one run) are used. The median,
     minimum and maximum over the repetitions are taken per (engine, recurrence, set, threads).
3. `code/python/compare.py` reads `raw/` and writes `comp_value.csv` of the three folders `charm`, `shtools`,
   `gshtools`, `raw/summary_tables.md` and the figures in `img/`.

The run matrix (per point set; 3 repetitions each):

| Engine (label in the files) | Configurations | parallel_mode |
|---|---|---|
| gsh | recurrence `fast` at 12 and 1 threads; `msfc`, `sfc`, `era` at 12 threads | `library` (OpenMP in the gsh library) |
| graflab | `sfc`, `msfc`, `era` with MATLAB's default thread count (not set by the script) | `library` (MATLAB's built-in multithreading) |
| pyshtools | 1 thread | `none` (serial library; the value of `threads` is the one requested) |
| shtools | 1 thread (serial loop of the driver) | `none` |
| pyharm | 12 and 1 threads | `library` (CHarm's OpenMP inside one batched call) |
| charm (`native-charm`) | `--parallel charm` at 12 and 1 threads | `library` |

`parallel_mode` (column of `res_time.csv`): `library` = the parallelism is inside a library call (the tool's own
OpenMP or built-in threads); `none` = no parallelism, the value of `threads` is the number of threads requested
(1). The diagnostic CHarm mode `--parallel points` (`native-points`, parallel_mode `driver_points`) is not in the
published matrix; if its rows exist in `raw/` they are consolidated, but no comparison uses them.

## Column definitions

### `res_value.csv` (folders `graflab`, `charm`, `shtools`, `gshtools`)

All points of all configurations of the folder, production run of each configuration (largest thread count).
One row per engine x recurrence x set x point.

| column | unit | meaning |
|---|---|---|
| engine | - | gsh, graflab, pyshtools, shtools, pyharm or charm |
| recurrence | - | fast, msfc, sfc, era (gsh, GrafLab); native (pyshtools, shtools, pyharm); native-charm (CHarm) |
| set | - | colorado or latsweep |
| id | - | point id, as in `data/points_<set>.csv` |
| lat_deg | deg | geodetic latitude on GRS80 (echoed from the input) |
| lon_deg | deg | longitude (echoed) |
| h_m | m | ellipsoidal height on GRS80 (echoed) |
| V_m2s2 | m^2/s^2 | gravitational potential V of the model at the point (see above) |
| T_m2s2 | m^2/s^2 | disturbing potential T = V - U |
| dg_mgal | mGal | gravity disturbance |

The values are written with 15 significant digits, exactly as the run wrote them. GrafLab has no thread tag in its
file name; its single value file is the one used.

### `res_time.csv` (one file per folder)

One row per engine x recurrence x set x threads. Times are of the synthesis only (model loading is not included;
`load_s_median` gives the loading time from the notes of each run).

| column | unit | meaning |
|---|---|---|
| engine | - | as in `res_value.csv` |
| recurrence | - | as in `res_value.csv` |
| set | - | colorado or latsweep |
| threads | - | threads requested by the run (for GrafLab: MATLAB's `maxNumCompThreads`, as recorded by the run) |
| parallel_mode | - | `library` or `none`, see above |
| n_repeats | - | number of repetitions the medians are taken over (3) |
| wall_s_median | s | median of the wall time of the synthesis, over the repetitions |
| wall_s_min | s | minimum over the repetitions |
| wall_s_max | s | maximum over the repetitions |
| wall_ms_median | ms | `wall_s_median` in milliseconds |
| cpu_s_median | s | median of the CPU time of the synthesis (all threads together; for `none` about equal to the wall time) |
| load_s_median | s | median of the model-loading time from the run's notes (empty when not recorded) |
| notes | - | what the engine's wall time contains (text, per engine) |

For GrafLab, the wall time is the total time of a `graflab_evaluate` call minus the time of a separate 1-point call
through the same pipeline (which measures the model loading). Because the load time is measured separately, the
subtraction is only meaningful when the synthesis dominates the call (true at nmax 2190 with 509 points, not at
the smoke-test sizes).

### `comp_value.csv` (folders `charm`, `shtools`, `gshtools`)

Pairwise differences, matched by point id. For each pair, point set and quantity: d = first - second at each point.

| column | unit | meaning |
|---|---|---|
| pair | - | `<engine a>_vs_<engine b>`; d = a - b |
| recurrence_a | - | recurrence (or CHarm mode) of a |
| recurrence_b | - | recurrence of b |
| set | - | colorado or latsweep |
| quantity | - | V, T or dg |
| unit | - | m^2/s^2 (V, T) or mGal (dg) |
| n | - | number of common points used |
| mean_d | unit of quantity | mean of d over the n points |
| rms_d | unit of quantity | root mean square of d |
| min_d | unit of quantity | smallest d |
| max_d | unit of quantity | largest d |
| max_abs_d | unit of quantity | largest |d| (the accuracy figure of the benchmark) |

Pairs in each folder:

* `charm/comp_value.csv`: `pyharm_vs_charm` (native vs native-charm); `pyharm_vs_graflab` (era);
  `charm_vs_graflab` (native-charm vs era).
* `shtools/comp_value.csv`: `pyshtools_vs_shtools` (native vs native); `pyshtools_vs_graflab` (era);
  `shtools_vs_graflab` (era).
* `gshtools/comp_value.csv`: `gsh_vs_graflab` with recurrence pairs sfc-sfc, msfc-msfc, era-era (gsh recurrence
  against GrafLab of the same recurrence) and fast-era (gsh fast against GrafLab era).

Where a comparison is against GrafLab of the same recurrence, the reference is named in `recurrence_b`. The
single-recurrence engines are compared with GrafLab era, the reference recurrence.

## Checks written with the files

* `raw/summary_tables.md`: differences of V, T, dg for every pair; run-time medians, minima and maxima per
  configuration; the thread-consistency table (the same configuration run with another thread count, compared with
  the production run; should be zero or round-off); warnings (a 1-thread run whose cpu/wall time exceeds 1.2, which
  would mean the thread count was not honoured).
* `img/`: `accuracy_max_abs_diff_colorado.png` (max|d| against GrafLab), `timing_bars_colorado.png` (run times),
  `latitude_sweep_V_gsh_vs_graflab.png` (|d| of V along the meridian), `speed_1_vs_max_threads_colorado.png`
  (1 thread against the largest thread count). Numbers printed on the figures are computed from the data.

## Stale files in `results/` (not part of this layout)

The files directly in `results/` other than `Summary.md` and `csv.md` (`*_colorado*.csv`, `*_colorado*.md`,
`timings.csv`, `res_values_colorado.csv`, `res_times_colorado.csv`) come from the earlier layout. They are not
read by the scripts of this layout and must be removed or moved before publication. The publication check at the
end of `run_all.ps1` fails while any `results/*.md` contains a local path.

## Reproducibility and checks (added with the run matrix)

* Results are reproducible up to round-off, not bit for bit. In the smoke test, gsh `fast` at 1 and 12 threads
  differed in dg at the 1e-10 mGal level, and pyharm at 1 and 12 threads differed in V at the level of the
  15-significant-digit output (about 1e-7 m^2/s^2 at V ~ 6e7). The thread-consistency table in `summary_tables.md`
  reports these differences.
* The `-OnlyGsh` regression check compares the new gsh max|d| against the stored `gshtools/comp_value.csv` with
  tolerance max(1e-9 relative, 1e-9 in the quantity unit). Values below this floor are not treated as changes.
* A 1-thread run whose cpu/wall ratio exceeds 1.2 is flagged in `summary_tables.md`. At smoke-test size the ratio is
  dominated by timer noise; at full size it is a real check that the run was serial (shtools, pyshtools, pyharm and
  charm at 1 thread). Any flag at full size needs a look before the timings are published.
* CHarm's diagnostic `--parallel points` is not in the matrix and has no comparison.
* The CHarm driver's default mode is now `charm`. The compiled binary predates this change and must be rebuilt
  (`code/c/build_charm_driver.sh`) before the next run. The official run passes `--parallel charm` explicitly, so
  the old binary still gives the official behaviour.
