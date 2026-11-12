# Spherical-harmonic synthesis benchmark (SHS)

Compares the spherical-harmonic synthesis of the gravity field model XGM2019e (nmax 2190) at two point sets, with
four engines and their official distributions: **gsh** (this project, C++), **GrafLab 2.1.4** (MATLAB),
**pyshtools 4.14.1** and **SHTOOLS 4.14.1** (serial), **pyharm 0.4.11** and **CHarm 0.4.11** (one batched call).
Each run computes the gravitational potential V, the disturbing potential T and the gravity disturbance dg.

Definitions, units and the layout of every CSV: [`results/csv.md`](results/csv.md). Inputs and the origin of the
model and of the points: [`data/DATA.md`](data/DATA.md). The summary of the results is `results/Summary.md`.

## Contents

```
benchmark/shs/
  README.md            this file
  data/                points (colorado, latsweep), SHA256SUMS.txt, DATA.md (the model .xz is a Release asset)
  code/
    run_all.ps1        the whole benchmark, one command
    run_latsweep.ps1   the latitude-sweep set only (wrapper of run_all.ps1)
    cpp/               gsh driver (CMake target gsh_shs_bench)
    c/                 CHarm driver (README and build script inside)
    fortran/           SHTOOLS driver (README and build script inside)
    python/            run_pyshtools.py, run_pyharm.py, consolidate.py, compare.py, shs_common.py,
                       requirements.txt, make_latitude_sweep.py, prepare_data.py (provenance only)
    matlab/            run_graflab.m (+ helpers); GrafLab.m is NOT included
  results/
    csv.md             definitions of the deliverable CSVs
    csv/               deliverable CSVs, one folder per group of engines
    raw/               per-run outputs and timings.csv (produced by run_all.ps1)
    Summary.md         written by hand from raw/summary_tables.md
  img/                 figures written by compare.py
```

## 1. Prerequisites

Each engine is built or installed separately. Versions matter; the ones below are those of the published runs.

| Engine | What to install | How |
|---|---|---|
| gsh | Visual Studio with clang-cl (preset `clang`), CMake | from the repository root: `cmake --preset clang` and `cmake --build --preset clang --target gsh_shs_bench` produce `build/clang/bin/gsh_shs_bench.exe` (the `gsh_shs_bench` target) |
| pyshtools, pyharm (and numpy, matplotlib) | Python 3.12 | `python -m venv .venv`, then `.venv\Scripts\python.exe -m pip install -r benchmark\shs\code\python\requirements.txt` (pins: numpy 2.5.4, matplotlib 3.11.2, pyshtools 4.14.1, pyharm 0.4.11) |
| CHarm driver | MSYS2 UCRT64 (GCC, FFTW, make) | see `code/c/README.md`; builds `build/external/charm-driver/charm_shs_bench.exe` and `libcharm.dll` |
| SHTOOLS driver | MSYS2 UCRT64 (gfortran, FFTW, LAPACK, BLAS) | see `code/fortran/README.md`; builds `build/external/shtools_bench/shtools_shs_bench.exe` |
| GrafLab | MATLAB R2026b and GrafLab 2.1.4 (`GrafLab.m`) | **obtained by the user, not redistributed**. Put the folder that contains `GrafLab.m` in the environment variable `GRAFLAB_DIR` (or pass `-GrafLabDir`) |

Environment variables read by the scripts (none of them has a value that belongs to a particular machine):

| Variable | Used for | Default when unset |
|---|---|---|
| `GRAFLAB_DIR` | folder with `GrafLab.m` (required for the GrafLab stage) | none: the script stops with a message |
| `MATLAB_EXE` | `matlab.exe` | the `matlab` found on `PATH` |
| `MSYS2_UCRT64_BIN` | the `ucrt64\bin` folder of MSYS2 (runtime DLLs of the CHarm and SHTOOLS drivers) | none: required for those two engines |

The machine: the published runs use 12 logical threads (6 cores). The run needs the machine otherwise idle: close
browsers, IDEs, synchronisation clients and any build. `run_all.ps1` prints the busiest processes at start.

## 2. Input data

The model is distributed compressed, as an asset of the GitHub Release `data-v1` of the public repository (it is not
in the git tree; see `data/DATA.md`). Download it into `benchmark/shs/data/`. The derived tide-free model is needed
before any run:

```
xz -dk benchmark\shs\data\XGM2019e_2190_tide_free.gfc.xz
```

Check it against `data/SHA256SUMS.txt` (the first line refers to the original zero-tide model, which is not in the
repository; check lines 2 and 3 only):

```
Get-FileHash -Algorithm SHA256 benchmark\shs\data\XGM2019e_2190_tide_free.gfc
```

The result must match the second line of `SHA256SUMS.txt`. The original zero-tide `XGM2019e_2159.gfc` is from ICGEM;
see `data/DATA.md`.

## 3. Run

From the repository root, PowerShell 7:

```
pwsh benchmark/shs/code/run_all.ps1 -Sets colorado,latsweep -GrafLabDir <folder with GrafLab.m>
```

Useful options:

| Option | Meaning |
|---|---|
| `-Sets` | `colorado`, `latsweep`, or both (default colorado) |
| `-Repeat` | timed repetitions per configuration (default 3) |
| `-SkipEngines` | any of `gsh, graflab, pyshtools, pyharm, shtools, charm` |
| `-OnlyGsh` | only the gsh stages; the frozen references already in `results/raw` are kept; the new gsh numbers are checked against the stored comparison (`results/csv/gshtools/comp_value.csv`) |
| `-DryRun` | print the plan and every command line; run nothing |
| `-Append` | keep an existing `results/raw/timings.csv` and add to it |

`code/run_latsweep.ps1` runs the latsweep set only and forwards the other options.

The run matrix (per point set, 3 repetitions each):

* gsh: `fast` at 12 and 1 threads; `msfc`, `sfc`, `era` at 12 threads.
* GrafLab: `sfc`, `msfc`, `era` with MATLAB's default threads (recorded in the timings).
* pyshtools: 1 thread. SHTOOLS driver: 1 thread. Both are serial libraries.
* pyharm: 12 and 1 threads. CHarm (`--parallel charm`): 12 and 1 threads.

The script stops at the first failing command (non-zero exit code) and prints the elapsed time after every run.
At the end it runs `consolidate.py` and `compare.py`, then a publication check that fails if `results/csv`,
`results/*.md` or `img` contain a local path or a personal name.

Cost: the stages are sequential. Most of the time goes to GrafLab (each repetition reloads the model), pyshtools
(serial, about four minutes per repetition at nmax 2190 with 509 points on the reference machine) and the SHTOOLS
and CHarm drivers. The `-Repeat` value multiplies every stage. The script prints the time of each run, so the
duration of the full run is known only when it finishes.

### Outputs

* `results/raw/`: one value file and one companion `.md` per run (`<engine>_<rec>[_<mode>]_<set>[_t<threads>].csv`),
  `timings.csv` (one row per repetition), `summary_tables.md`, `run_all.log` (transcript).
* `results/csv/<folder>/`: `res_value.csv`, `res_time.csv` and, for `charm`, `shtools`, `gshtools`, `comp_value.csv`.
  Column definitions in `results/csv.md`.
* `img/`: four figures (accuracy, run times, latitude sweep, 1 thread against the largest thread count).

To rebuild only the tables and figures from `results/raw` (no engine is run):

```
.venv\Scripts\python.exe benchmark\shs\code\python\consolidate.py
.venv\Scripts\python.exe benchmark\shs\code\python\compare.py
```

## 4. Reading the results

* `results/csv.md` explains every file and column.
* `results/raw/summary_tables.md` holds the tables of differences (V, T, dg), the run times and the
  thread-consistency check. `results/Summary.md` is the written summary.
* The differences are d = first - second at each point, matched by point id; the reference for accuracy is GrafLab
  (same recurrence for gsh; era for the single-recurrence engines).

## 5. Known limits

* `XGM2019e_2159.gfc` (zero tide) is not in the repository; the derived model is the input of every engine.
* The GrafLab wall time is obtained by subtracting a separately measured model load (see `results/csv.md`). It is not
  meaningful at the small smoke-test sizes.
* `code/python/prepare_data.py` is a provenance script. It imports the private project package and is not runnable from the
  public tree.
