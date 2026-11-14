# GrafLab grid runner (MATLAB), regional case

`run_graflab_grid.m` runs GrafLab 2.1.4 in **grid mode** (`point_type = 0`) on the regional grid of the grid
benchmark: latitude 45 to 35 deg, longitude -110 to -100 deg, step 0.05 deg, 201 x 201 nodes, at the
constant radius R of the model header. It computes V (functional 11) and T (functional 5) in one call per
repetition, for one of the three Legendre-function options:

| recurrence | GrafLab `fnALFs` | meaning |
|---|---|---|
| `sfc` | 1 | standard forward column |
| `msfc` | 2 | modified forward column |
| `era` | 3 | extended-range arithmetic |

The global cases (`gl`, `dh1`, `dh2`) are **refused** with an error: GrafLab is not run on the 9.6 million-node
grids (memory and hours).

## GrafLab is not in this repository

`GrafLab.m` (Bucha and Janak, 2013, Computers and Geosciences 56, 186-196,
doi:10.1016/j.cageo.2013.03.012) is not copied here. Pass its folder as the last argument of the function. Check
its terms of use at its source before redistributing anything.

## Call

```matlab
% from the repository root, MATLAB with GrafLab 2.1.4 on the path folder given below
mkdir('build/grid_out')
addpath('benchmark/grid/code/matlab')
run_graflab_grid('benchmark/shs/data/XGM2019e_2190_tide_free.gfc', 2190, 'regional', 1, ...
                 'build/grid_out', 1, [], 'C:\path\to\GrafLab_folder')
```

Arguments: `gfc_path, nmax, grid_case, fnALFs, out_dir, nrep, rows, graflab_dir`. `rows` (0-based, empty for the
default `[0 50 100 150 200]`) selects the latitude rows written; row 0 is 45 deg north.

The `.gfc` must be the decompressed tide-free file (`xz -dk benchmark/shs/data/XGM2019e_2190_tide_free.gfc.xz`).

## Coordinate type (checked in GrafLab.m)

The common benchmark definition is a **spherical (geocentric) latitude** at radius R. GrafLab supports this
directly, so no mismatch with the other engines arises:

* the input description (GrafLab.m, function header, around lines 282-322): `coord` 0 = ellipsoidal,
  1 = spherical; for `coord == 1` the reference surface is a sphere of radius `R`, and `h` is the height above it
  in the radial direction. The script uses `coord = 1, h = 0`, so r = R.
* grid branch (GrafLab.m, around line 2346): for `coord == 1` the radius is `r = (R + h)` and the latitude is
  the spherical one (`fiG = fi`).
* export (GrafLab.m, around lines 8912-8930): for `coord == 1` the exported latitude is the spherical latitude
  (`fi = fiG`), so the returned latitude column is the same geocentric latitude as the benchmark node.

For `coord == 0` (ellipsoidal) GrafLab would export geodetic latitude; this runner does not use it.

## Output and mapping

GrafLab returns, for `nargout == 1` in grid mode, a matrix with one row per node and the columns
`[lat_deg, lon_deg, V, T]` (GrafLab.m, grid branch: `output = [fi; lambda; P]`, transposed). The script places
each returned value by its own latitude and longitude (row index `round((45 - lat)/0.05) + 1`, column index
`round((lon + 110)/0.05) + 1`), so GrafLab's internal node order is not assumed. It errors if a node is missing.

Files written in `out_dir` (never a whole grid):

* `graflab_regional_<recurrence>_rows.csv`: header `row,lat_deg,lon_deg,V_m2s2,T_m2s2`, %.17g, row 0-based.
  The file name carries the recurrence because the three GrafLab modes would otherwise overwrite each other.
  This differs from the `<engine>_<case>_rows.csv` pattern of the other drivers (noted in
  `../README.md`).
* `timings.csv`: one line per repetition, `graflab,regional,<recurrence>,<threads>,<nmax>,201,201,<rep>,<wall_s>,-1`.

## Known limits

* `threads` in the timings is `maxNumCompThreads` of MATLAB. GrafLab itself is MATLAB code; its parallelism is
  whatever MATLAB's BLAS and the vectorised code use. It is recorded, not controlled.
* `peak_mem_mb` is -1: the script does not measure memory.
* GrafLab's `T` (functional 5) uses GrafLab's own normal field, while the other engines use the closed-form
  `U` of `shs_common.py`. The SHS benchmark already reports the two as agreeing to about 1e-6 m^2/s^2
  (`shs_common.py` docstring). The difference is not removed here.
* The latitude and longitude written are the nominal grid values (`45 - i*0.05`, `-110 + j*0.05`), not the
  values GrafLab echoes, which agree to floating-point accuracy.

## Not verified

* This file was **written but never run**: GrafLab was not executed, no output was checked against GrafLab,
  and the call was not tested on the installed GrafLab version.
* The argument order was copied from `../../../shs/code/matlab/graflab_evaluate.m` (point mode), and the grid
  fields were checked against the function header of the GrafLab.m in the scratchpad folder. The positional list
  itself was not executed.
* The grid row order (fi fastest in the flattened meshgrid) was read from the source but not observed; the
  script does not depend on it.
