# Grid benchmark

Synthesis of the gravitational potential V, the disturbing potential T and the exact gravity disturbance dg on
latitude-longitude grids, for the gsh engine. The gsh side is described below. The definitions are the same for every
engine: spherical (geocentric) latitudes in radians internally, constant radius per row, and the same functional formulas.

## Definitions

Latitudes are SPHERICAL (geocentric). Every grid uses a constant radius r = R, the reference radius of the gravity model
(R in the `.gfc` header). Longitudes are in degrees east in the output.

| case | rows | columns | definition |
|---|---|---|---|
| `gl` | N + 1 | 2N + 2 | Gauss-Legendre latitudes (roots of P_{N+1}), full circle, uniform longitudes from 0 |
| `dh1` | 2(N + 1) | 2(N + 1) | Driscoll-Healy, first variant: latitudes 90 - i·180/(2(N+1)) deg, the equator row exactly 0; full circle (grid.hpp, Driscoll and Healy 1994) |
| `dh2` | 2(N + 1) | 4(N + 1) | the same latitudes as dh1 with twice the longitudes (grid.hpp) |
| `regional` | 201 | 201 | latitudes 45 to 35 deg (north to south) and longitudes -110 to -100 deg, step 0.05 deg, both ends included |

The regional grid is not a full circle; gsh takes its direct-summation path (grid.cpp).

## gsh side

Files (all under `benchmark/grid/`):

- `gsh_grid_bench.cpp`: the benchmark driver (target `gsh_grid_bench`). `grid_assembly.hpp` holds the shared
  assembly of V, T and dg from the raw sums, grid and row-selection helpers.
- `gsh_grid_check.cpp`: the accuracy check of the grid path against the point path (target `gsh_grid_check`).
- `tests/cpp/test_grid_highdegree.cpp` (in `tests/cpp/`): GoogleTest at nmax = 360.

### How V, T and dg are obtained

There is no grid functional API. The driver calls `gsh::synth_grid` and assembles the functionals with the formulas of
`src/cpp/gsh_expansion/functionals.cpp`:

- V = GM/r · S_0 of the model (`synth_grid` with grad 0);
- T = GM/r · S_0 of the disturbing model `disturbing_model(model, GRS80)` (`synth_grid` with grad 1 when dg is requested,
  which also yields S_0);
- dg = (2 γ·δ + |δ|²)/(|γ+δ| + γ), with δ = (GM/r² S_P, -GM/(r² cos φ) S_L, -GM/r² S_R) and γ = (γ_x, 0, γ_z) from
  `Ellipsoid::normal_gravity` at the spherical latitude of the row. dg is written as nan on rows with cos φ below the
  library pole threshold (5e-3, `FunctionalRequest::pole_cos_threshold`).

### CLI: `gsh_grid_bench`

```
gsh_grid_bench --gfc <file.gfc> --nmax N --case gl|dh1|dh2|regional --recurrence fast|msfc --threads T
               --repeat K --out-dir D --rows auto|i,j,k --fields V,T,dg
```

- `--rows auto` (global grids): the rows closest to 0, ±30, ±60, ±85 and ±89.5 deg, plus rows 0 and nlat - 1 (the poles
  for dh1 and dh2). For `regional`, `auto` selects rows 0, (nlat-1)/4, (nlat-1)/2, 3(nlat-1)/4 and nlat-1 (45, 42.5,
  40, 37.5 and 35 deg). Indices can be given as a comma-separated list.
- `--recurrence msfc` is rejected: `synth_grid` has only the fast kernel tables. MSFC exists only in the per-point
  oracle (`reference_point`), which is too slow for grids of this size. The MSFC grid is not measured.
- Each repetition calls `synth_grid` once for the model (V) and once for the disturbing model (T, dg). The fields not
  named in `--fields` are skipped.
- The full grid is never written to disk. Only the selected rows are assembled and written.

Outputs in `--out-dir`:

- `gsh_<case>_rows.csv`: columns `row,lat_deg,lon_deg,V_m2s2,T_m2s2,dg_ms2`; one line per longitude of each selected row,
  written with `%.17g`. `row` is the 0-based index in the grid, `lat_deg` the spherical latitude of that row. Fields
  not requested are `nan`.
- `timings.csv` (appended, one line per repetition):
  `engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,lumps_thread_s,fft_thread_s,peak_mem_mb`.
  `wall_s` is the sum of `GridTiming.total` over the `synth_grid` calls of the repetition (the assembly is not timed).
  `lumps_thread_s` and `fft_thread_s` are the matching sums of the thread-seconds. `peak_mem_mb` is the process peak
  working set (GetProcessMemoryInfo, Windows), which is cumulative over the process.

### CLI: `gsh_grid_check`

```
gsh_grid_check --gfc <file.gfc> [--out-dir D] [--cases gl,regional,dh1] [--nmax-high 2190] [--nmax-light 720]
               [--print-ref-quad]
```

For each case the full grid is synthesised with `synth_grid`, and the selected rows are assembled. The same rows are then
evaluated with `grid_to_points` (on the selected rows only) and `compute_functionals` with outputs V, T and dg. Per field
the check reports max |grid - points| / max |points| against the criterion 1e-12 (fixed beforehand), and the row of the
largest difference. Rows at an exact pole are excluded, because `compute_functionals` rejects dg there. Rows near a pole
(cos φ < 5e-3) are included, since both paths use the oracle there.

`--print-ref-quad` prints, without running it, the `ref_quad` command line for the selected Gauss-Legendre rows
(`benchmark/meridian/ref_quad`, built and run under WSL, `--engine-rad --lat-list`). ref_quad evaluates at longitude 0
only.

Output: a Markdown table on stdout and in `<out-dir>/grid_vs_points.md`.

### Accuracy (measured)

Measured with `gsh_grid_check` on the EGM2008 model (degree 2190) in a single run. This is an accuracy run, not a
timing run.

| case | nmax | rows | points | field | max abs diff | max abs value | relative diff | row of max diff (lat deg) | verdict |
|---|---|---|---|---|---|---|---|---|---|
| gl | 2190 | 11 | 48202 | V | 7.302e-07 | 6.253e+07 | 1.168e-14 | row 1095 (0.000000) | PASS |
| gl | 2190 | 11 | 48202 | T | 9.436e-12 | 1.011e+03 | 9.335e-15 | row 1095 (0.000000) | PASS |
| gl | 2190 | 11 | 48202 | dg | 9.634e-16 | 3.512e-03 | 2.743e-13 | row 1095 (0.000000) | PASS |
| regional | 2190 | 5 | 1005 | V | 2.086e-07 | 6.250e+07 | 3.338e-15 | row 200 (35.000000) | PASS |
| regional | 2190 | 5 | 1005 | T | 3.382e-12 | 2.716e+02 | 1.245e-14 | row 150 (37.500000) | PASS |
| regional | 2190 | 5 | 1005 | dg | 2.662e-17 | 1.222e-03 | 2.178e-14 | row 150 (37.500000) | PASS |
| dh1 | 720 | 10 | 14420 | V | 4.545e-07 | 6.253e+07 | 7.268e-15 | row 721 (0.000000) | PASS |
| dh1 | 720 | 10 | 14420 | T | 7.276e-12 | 1.011e+03 | 7.199e-15 | row 961 (-29.958391) | PASS |
| dh1 | 720 | 10 | 14420 | dg | 2.069e-16 | 2.449e-03 | 8.448e-14 | row 721 (0.000000) | PASS |

Wall time of each side in this single accuracy run (not a timing benchmark):

| case | grid side (s) | point side (s) |
|---|---|---|
| gl | 1.5 | 5344.2 |
| regional | 0.5 | 0.5 |
| dh1 | 0.2 | 107.3 |

Criterion max|grid - points| / max|points| <= 1e-12: all nine rows pass (worst 2.7e-13, dg on the Gauss-Legendre grid).
The Driscoll-Healy case was checked at degree 720 (the point side at 2190 costs hours); gl and regional at 2190.
Point-side wall times of that run: gl 5344 s, dh1 107 s, regional 0.5 s; grid side 1.5, 0.2 and 0.5 s. Full output:
`results/grid_vs_points.md`.

### Timing (measured 2026-11-13)

gsh `fast`, EGM2008 degree 2190, V + T + dg, 12 threads (6 cores), AVX-512, idle machine, 8 repetitions
(`gsh_grid_bench --repeat 8`; the first repetition is included). Per-repetition rows: `results/timings_<case>.csv`;
the rows used for the accuracy check: `results/gsh_<case>_rows.csv`.

| case | grid | nodes | wall median (s) [min, max] | lumps thread-s | FFT thread-s | peak memory (MB) |
|---|---|---|---|---|---|---|
| gl | 2191 x 4382 | 9.60 M | 1.309 [1.224, 1.529] | 5.14 | 6.25 | 728 |
| dh1 | 4382 x 4382 | 19.20 M | 3.269 [3.185, 3.785] | 15.33 | 12.85 | 1093 |
| regional | 201 x 201 | 0.04 M | 0.316 [0.288, 0.364] | 1.26 | 1.63 | 283 |

`lumps` and `FFT` are thread-seconds summed over the threads (not wall time). A global Gauss-Legendre grid of 9.6 million
nodes with V, T and dg takes about 1.3 s; the Driscoll-Healy grid, twice as large, about 3.3 s. The point-path
equivalent costs 5344 s for only 48 202 points (single accuracy run), which is why the check uses selected rows.

### Not measured

- Timing was measured for gsh only (above). No other engine was timed on grids.
- The `msfc` recurrence on grids: not implemented on the grid path (see above).
- The external engines (MATLAB/GrafLab, CHarm, ref_quad) have not been run for this grid.
