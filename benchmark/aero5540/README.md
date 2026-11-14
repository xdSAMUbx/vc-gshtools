# Aero 5540: GRAV-D flight, 45 functionals at degree/order 5540

Official large-scale demonstration of `gsh_aero5540`: every one of the 45 functionals of `gsh::compute_functionals`,
degree and order 5540, XGM2019e, on a real airborne flight. Anyone with the three input files can rerun it.

## What is computed

- **Flight:** GRAV-D (NGS) line over Colorado, **283 716 points** (1 Hz, about 5.4-6 km high). Only the coordinates are
  used; the GRAV-D data are not redistributed.
- **Model:** XGM2019e, degree/order 5540. The program reads the tide system from the file header; for this file it
  reports `zero_tide`.
- **Heights for N and zeta:** Colorado DEM (5' mean heights). Points outside the DEM give NaN for those two functionals.
- **Machine of the results below:** AMD Ryzen 7 7445HS (6 cores, 12 threads), Windows 11, clang-cl, OpenMP,
  active kernel **AVX-512** (`active_isa_name()` = `avx512`). The machine was idle during every timing.

## Protocol

| Item | Protocol |
|---|---|
| Timing, one functional at a time | **1 run**, each functional requested alone (so the time includes the synthesis pass it needs) |
| Timing, fused pass A/B | 8 repetitions, all recorded, the first included |
| Accuracy | 1 run, against the sequential double-precision oracle (`reference_point`) at the same degree |

A functional alone pays for its whole synthesis pass, and functionals of one group repeat that pass, so the 45
per-functional times are **not additive**. They show which functionals force which pass. The real cost of "everything
at once" is the all-together run (next section).

## Result 1: all 45 together, and the effect of fusing the model and disturbing passes

2013 flight points (every k-th point of the flight), nmax 5540, 8 repetitions, `results/h2_fuse/`:

| Mode | Median total (s) | Model pass | Disturbing pass | Ellipsoid pass (N, zeta) |
|---|---:|---:|---:|---:|
| Two passes (`fuse_model_and_disturbing = false`) | 17.37 | 6.69 s (grad 2) | 6.69 s (grad 2) | 3.70 s |
| Fused (default) | **10.85** | not separate | 6.79 s | 3.71 s |

The disturbing field differs from the full model in 11 zonal coefficients only. The fused path synthesises the
disturbing set once at the highest gradient order needed and adds a short zonal series (degree at most 20) to get the
model sums; the sums are added, never subtracted, so there is no cancellation. **37.5 % less time** with the same
peak memory (1953 MB at 2013 points).

Fused vs two-pass results on the same points, max|d| / max|ref| per functional: worst case **2.2e-14**
(`results/h2_fuse/fuse_vs_unfused.csv`), far below the 1e-12 criterion.

**Supposition (not measured):** if the time stays proportional to the number of points, the 45 together on the whole
flight would take roughly 25 minutes (10.85 s for 2013 points). This run was not made at full size.

## Result 2: assumed time of each functional on the whole flight

**These are assumptions, not measurements of the whole flight.** Each functional was requested alone, once, on a
10 133-point subsample at degree/order 5540; the table assumes that the time grows in proportion to the number of points
(283 716). Each functional pays for its whole synthesis pass, so the times must not be added to get the cost of
several functionals: for that, request them in one call (Result 1).

Pass types: grad 0 (V, W, T, zeta_ell), grad 1 (first derivatives, g, dg, xi, eta, Theta, N, zeta) and grad 2 (second
derivatives); they cost about 4, 9 and 15 minutes on the whole flight under this assumption.

| Functional | Assumed time, whole flight (min) |
|---|---:|
| V | 4 |
| dV_r | 9 |
| dV_phi | 9 |
| dV_lambda | 9 |
| W | 4 |
| g_x | 8 |
| g_y | 8 |
| g_z | 8 |
| g | 8 |
| g_sa | 8 |
| W_rr | 15 |
| V_rr | 15 |
| V_rphi | 15 |
| V_rlambda | 15 |
| V_phiphi | 15 |
| V_philambda | 15 |
| V_lambdalambda | 15 |
| V_xx | 15 |
| V_xy | 15 |
| V_xz | 15 |
| V_yy | 15 |
| V_yz | 15 |
| V_zz | 15 |
| T | 4 |
| dg_sa | 8 |
| Dg_sa | 8 |
| T_rr | 15 |
| dg | 8 |
| T_rphi | 15 |
| T_rlambda | 15 |
| T_phiphi | 15 |
| T_philambda | 15 |
| T_lambdalambda | 15 |
| T_xx | 15 |
| T_xy | 15 |
| T_xz | 15 |
| T_yy | 15 |
| T_yz | 15 |
| T_zz | 15 |
| xi | 8 |
| eta | 8 |
| Theta | 8 |
| zeta_ell | 4 |
| N | 8 |
| zeta | 8 |

Raw measurements: `results/subsample_10133/timings_alone.csv` (all 45, 10 133 points) and
`results/full_flight_283716/timings_alone.csv` (16 functionals that were also run on the whole flight before the run
was stopped).

## Result 3: accuracy at degree 5540

200 flight points (201 used), V and the six local tensor components against the oracle. Metric max|d| / max|ref|; the
criterion of `benchmark/meridian` is 1e-12.

| Quantity | max abs diff | max abs ref | Relative |
|---|---:|---:|---:|
| V | 9.46e-07 | 6.25e+07 | 1.5e-14 |
| V_xx | 1.50e-18 | 1.59e-06 | 9.5e-13 |
| V_xy | 8.16e-21 | 3.22e-08 | 2.5e-13 |
| V_xz | 1.88e-20 | 7.56e-08 | 2.5e-13 |
| V_yy | 9.59e-19 | 1.59e-06 | 6.0e-13 |
| V_yz | 1.58e-20 | 8.29e-08 | 1.9e-13 |
| V_zz | 2.10e-18 | 3.18e-06 | 6.6e-13 |

V is at 1.5e-14. The tensor components are at most 9.5e-13 (V_xx): inside the criterion, with less than a factor 2
of margin on that one, so do not assume the criterion holds for every functional at every point. The oracle is the
library's own sequential code, so this tests the optimized kernels, not the formulas. A binary128 check at degree 5540
was not done.

## Not measured

- All 45 together on the whole flight (only the supposition in Result 1).
- Most of the functionals alone on the whole flight (Result 2 is an assumption).
- Accuracy of the other 38 functionals at 5540 (the tensor and V share the kernels with them; not checked one by one).

## Inputs (not included in the repository)

| Input | Where to get it | Format |
|---|---|---|
| Model `XGM2019e.gfc` | ICGEM, <https://icgem.gfz-potsdam.de/> (XGM2019e, degree 5540) | ICGEM `.gfc` |
| Flight | NGS GRAV-D airborne data, <https://www.ngs.noaa.gov/GRAV-D/> (the file used here is the median-debiased 1 Hz product) | comma-separated, one header line; columns 4, 5 and 6 (counting from 1) are geodetic latitude (deg), longitude (deg) and ellipsoidal height (m) |
| DEM `Colorado_dem.grd` | any DEM block-averaged to a GRAVSOFT grid in metres (the one used here: SRTM, 5' means) | GRAVSOFT grid, read by `load_gravsoft_grid` |

Any other flight file in that column layout works; N and zeta need the DEM only for those two functionals.

## Run

```
build/clang/bin/gsh_aero5540 <model.gfc> <flight.txt> <dem.grd> <out_dir> [nmax=5540] [reps=8] [sub_points=10000]     [acc_points=200] [full_reps=reps] [full_points=0=all] [resume=0] [ab=0]
```

| Argument | Meaning |
|---|---|
| `reps` | repetitions of each functional alone |
| `sub_points` | points of the per-functional stage; every k-th point. `0` skips the stage. `1000000` = every point |
| `acc_points` | points of the accuracy stage. `0` skips it |
| `full_reps` | repetitions of the all-together run; `0` skips it (no `timings_full.csv`) |
| `full_points` | points of the all-together run (`0` = all) |
| `resume` | `1` keeps the rows already in `timings_alone.csv` and computes only the missing functionals |
| `ab` | `1` repeats the all-together run with the fused pass on and off and writes `fuse_vs_unfused.csv` |

Official full-flight, one run per functional, nothing else:

```
gsh_aero5540 data/cache/modelos/XGM2019e/XGM2019e.gfc   data/cache/tablas/GRAVD_ms05_median_debiased_1hz/GRAVD_ms05_median_debiased_1hz.txt   data/cache/tablas/Colorado_dem/Colorado_dem.grd out/aero5540 5540 1 1000000 200 0 0
```

The fused-pass comparison of Result 1:

```
gsh_aero5540 <model> <flight> <dem> out/h2 5540 1 0 0 8 2000 0 1
```

Memory: about 2 GB at 2000 points; keep other programs closed while timing. Output files: `timings_alone.csv`
(`functional,rep,points,wall_s`), `timings_full.csv` (per stage, with a `fuse` column), `accuracy.csv`.
