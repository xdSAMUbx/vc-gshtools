# GrafLab runner (MATLAB)

Runs GrafLab 2.1.4 headless to compute, in one call per repetition, the gravitational potential
V (functional 11), the disturbing potential T (functional 5) and the gravity disturbance dg
(functional 20) at the benchmark points (`Functional = [11; 5; 20]`), with each of its three
Legendre-function options:

| recurrence | GrafLab `fnALFs` | meaning |
|---|---|---|
| `sfc` | 1 | standard forward column |
| `msfc` | 2 | modified forward column |
| `era` | 3 | extended-range arithmetic |

## GrafLab is not in this repository

`GrafLab.m` (Bucha and Janak, 2013, Computers and Geosciences 56, 186-196,
doi:10.1016/j.cageo.2013.03.012) is not copied here. Its header has no licence text, so
redistribution is not assumed to be allowed. Obtain it from its author (the header refers to
the `graflab-cookbook` repository by Blazej Bucha) and check the terms of use there.

Point `run_graflab.m` at the folder that contains `GrafLab.m`, either:

- as the 8th argument of `run_graflab`, or
- with the environment variable `GRAFLAB_DIR`.

## Requirements

- MATLAB (tested target: R2026b, `C:\Program Files\MATLAB\R2026b\bin\matlab.exe`).
- The tide-free ICGEM `.gfc` (`benchmark/shs/data/XGM2019e_2190_tide_free.gfc`).
- A points CSV with header `id,lat_deg,lon_deg,h_m` (GRS80 geodetic latitude, ellipsoidal height).

## Usage

From PowerShell, in `benchmark\shs\code\matlab`:

```powershell
& "C:\Program Files\MATLAB\R2026b\bin\matlab.exe" -batch "run_graflab('C:\...\XGM2019e_2190_tide_free.gfc', 'C:\...\points_colorado.csv', 360, 'sfc', 'colorado', 'C:\...\benchmark\shs\results', 1, 'C:\...\GrafLab_folder')"
```

Arguments, in order: `gfc_path`, `points_csv`, `nmax`, `recurrence` (`sfc|msfc|era`),
`set_name`, `out_dir` (the `results` folder), `repeat_count`, `graflab_dir` (optional if
`GRAFLAB_DIR` is set).

Outputs:

- `results/graflab_<recurrence>_<set>.csv` with header `id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal`.
- `results/graflab_<recurrence>_<set>.md` companion (command, parameters, limitations, timing).
- One row per repetition appended to `results/timings.csv`.

## Files

- `run_graflab.m`: entry point, timing, output.
- `graflab_evaluate.m`: one headless GrafLab call and the mapping of its arguments.
- `read_gfc_header.m`: reads GM, R, degree, `norm` and `tide_system` from the `.gfc` header.
- `read_points_csv.m`: reads the points file.
- `append_timing_row.m`: appends to `timings.csv`.
- `write_result_files.m`: writes the result CSV and its `.md`.

See `graflab_evaluate.m` for the GrafLab-specific conventions (GRS80, `nmin = 0`, point-wise mode).
