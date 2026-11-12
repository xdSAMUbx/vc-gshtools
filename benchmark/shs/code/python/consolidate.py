"""Consolidate the raw benchmark outputs into the tidy CSV folders of results/csv.

Input (default benchmark/shs/results/raw, written by code/run_all.ps1):

  <engine>_<recurrence>[_<mode>]_<set>[_t<threads>].csv
      value file, header id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal (one run = V, T and dg together).
      engine in gsh | graflab | pyshtools | shtools | pyharm | charm
      recurrence in fast | msfc | sfc | era (gsh, graflab) or native (all others)
      mode in points | charm (charm only: the value of --parallel)
      _t<threads> absent for GrafLab (its thread count is MATLAB's own, see timings.csv)
      Legacy names of the Colorado results (<engine>_<recurrence>_colorado.csv) are read too.
  timings.csv
      engine,recurrence,set,nmax,n_points,threads,repeat,wall_s,cpu_s,notes
      only rows tagged functional=V_T_dg are used. For charm the mode is read from the notes.

Output (default benchmark/shs/results/csv), one folder per group of engines:

  graflab/   res_value.csv, res_time.csv
  charm/     res_value.csv, res_time.csv   (engines pyharm and charm; charm recurrence = native-points | native-charm)
  shtools/   res_value.csv, res_time.csv   (engines pyshtools and shtools)
  gshtools/  res_value.csv, res_time.csv   (engine gsh)

comp_value.csv of each folder is written by compare.py, which imports the readers of this module.
Both point sets go into the same files, distinguished by the `set` column.

Usage (repository root):
    .venv\\Scripts\\python.exe benchmark\\shs\\code\\python\\consolidate.py
    .venv\\Scripts\\python.exe benchmark\\shs\\code\\python\\consolidate.py --raw <dir> --csv <dir>
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

SHS_DIR = Path(__file__).resolve().parents[2]          # benchmark/shs
DEFAULT_RAW = SHS_DIR / "results" / "raw"
DEFAULT_CSV = SHS_DIR / "results" / "csv"

RESULT_COLUMNS = ["id", "lat_deg", "lon_deg", "h_m", "V_m2s2", "T_m2s2", "dg_mgal"]
VALUE_COLUMNS = ["V_m2s2", "T_m2s2", "dg_mgal"]
FUNCTIONAL_TAG = "functional=V_T_dg"

ENGINES = ("gsh", "graflab", "pyshtools", "shtools", "pyharm", "charm")
ENGINE_RANK = {e: i for i, e in enumerate(ENGINES)}
SETS = ("colorado", "latsweep")
SET_RANK = {s: i for i, s in enumerate(SETS)}
RECURRENCES = ("fast", "msfc", "sfc", "era", "native", "native-points", "native-charm")
REC_RANK = {r: i for i, r in enumerate(RECURRENCES)}
FOLDER_OF = {"gsh": "gshtools", "graflab": "graflab", "pyharm": "charm", "charm": "charm",
             "pyshtools": "shtools", "shtools": "shtools"}
FOLDERS = ("graflab", "gshtools", "charm", "shtools")

# Parallel mode of each engine (and of each CHarm mode). Meaning, per engine, is written in csv.md.
PARALLEL_MODE = {
    ("gsh", None): "library",           # OpenMP inside the gsh library, threads set by --threads
    ("graflab", None): "library",       # MATLAB built-in multithreading (maxNumCompThreads), not set by the script
    ("pyshtools", None): "none",        # serial library, run with 1 thread
    ("pyharm", None): "library",        # CHarm OpenMP inside one batched shs.point call
    ("shtools", None): "none",          # serial driver loop, run with 1 thread (driver OpenMP not used)
    ("charm", "points"): "driver_points",  # DIAGNOSTIC only (not in run_all.ps1): OpenMP over points in the driver
    ("charm", "charm"): "library",      # official: one batched CHarm call, CHarm OpenMP over orders m
}

NOTES = {
    "gsh": "gsh_shs_bench; wall = synthesis only; threads = --threads of the run",
    "graflab": "GrafLab in MATLAB; wall = total call minus a 1-point probe (model load); model reloaded in every "
               "repetition; threads = maxNumCompThreads of MATLAB",
    "pyshtools": "pyshtools 4.14.1 point routines, serial (1 thread); V needs one radius-scaled coefficient copy per "
                 "distinct radius (inside the timed block)",
    "pyharm": "pyharm 0.4.11 shs.point + point_grad1 (one batched call); OMP_NUM_THREADS set before import",
    "shtools": "SHTOOLS 4.14.1 Fortran driver, serial loop (1 thread); wall = synthesis only; model load excluded",
    "charm": "CHarm 0.4.11 C driver, one batched call (--parallel charm); wall = synthesis only; load excluded",
}

VALUE_NAME = re.compile(
    r"^(?P<engine>gsh|graflab|pyshtools|shtools|pyharm|charm)"
    r"_(?P<rec>sfc|msfc|era|fast|native)"
    r"(?:_(?P<mode>points|charm))?"
    r"_(?P<set>colorado|latsweep)"
    r"(?:_t(?P<threads>[0-9]+))?$")
MODE_RE = re.compile(r"parallel=(points|charm)")
RECURRENCE_OK = {"gsh": {"fast", "msfc", "sfc", "era"}, "graflab": {"sfc", "msfc", "era"}}


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def recurrence_label(engine: str, rec: str, mode: Optional[str]) -> str:
    """Label used in the csv files: CHarm modes are part of the recurrence column."""
    if engine == "charm":
        return f"{rec}-{mode}"
    return rec


def parallel_mode(engine: str, mode: Optional[str]) -> str:
    return PARALLEL_MODE[(engine, mode if engine == "charm" else None)]


# ----------------------------------------------------------------------------------------
# Value files
# ----------------------------------------------------------------------------------------
@dataclass
class ValueFile:
    engine: str
    rec: str                      # recurrence label (see recurrence_label)
    set_name: str
    threads: Optional[int]        # None for GrafLab (no threads tag)
    path: Path
    ids: list
    lat: list
    lon: list
    h: list
    V: list
    T: list
    dg: list
    rows: list = field(repr=False)  # raw text rows, column order RESULT_COLUMNS

    @property
    def key(self) -> tuple:
        return (self.engine, self.rec, self.set_name)

    def column(self, name: str) -> list:
        return {"V_m2s2": self.V, "T_m2s2": self.T, "dg_mgal": self.dg}[name]


def parse_value_name(stem: str):
    m = VALUE_NAME.match(stem)
    if m is None:
        return None
    engine, rec, mode, set_name = m["engine"], m["rec"], m["mode"], m["set"]
    if engine == "charm" and mode is None:
        return None
    if engine != "charm" and mode is not None:
        return None
    if engine in RECURRENCE_OK and rec not in RECURRENCE_OK[engine]:
        return None
    if engine not in RECURRENCE_OK and rec != "native":
        return None
    threads = int(m["threads"]) if m["threads"] is not None else None
    return engine, recurrence_label(engine, rec, mode), set_name, threads


def read_value_file(path: Path) -> Optional[ValueFile]:
    meta = parse_value_name(path.stem)
    if meta is None:
        warn(f"skipping {path.name}: not a value file name (see the module docstring)")
        return None
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = [c.strip() for c in next(reader, [])]
        if header != RESULT_COLUMNS:
            warn(f"skipping {path.name}: header {header} differs from {RESULT_COLUMNS}")
            return None
        rows = [r for r in reader if r and any(c.strip() for c in r)]
    engine, rec, set_name, threads = meta
    return ValueFile(engine, rec, set_name, threads, path,
                     [r[0].strip() for r in rows], [float(r[1]) for r in rows],
                     [float(r[2]) for r in rows], [float(r[3]) for r in rows],
                     [float(r[4]) for r in rows], [float(r[5]) for r in rows],
                     [float(r[6]) for r in rows], rows)


def load_value_files(raw: Path) -> list:
    files = []
    for p in sorted(raw.glob("*.csv")):
        if p.name == "timings.csv":
            continue
        vf = read_value_file(p)
        if vf is not None:
            files.append(vf)
    return files


def select_values(files: list) -> dict:
    """One value file per (engine, recurrence, set): the run with the most threads (the production run).

    Other thread counts of the same configuration are kept only for the thread-consistency check.
    """
    chosen = {}
    for vf in files:
        cur = chosen.get(vf.key)
        if cur is None or (vf.threads or 0) > (cur.threads or 0):
            chosen[vf.key] = vf
    return chosen


# ----------------------------------------------------------------------------------------
# Timings
# ----------------------------------------------------------------------------------------
@dataclass
class TimingRow:
    engine: str
    rec: str
    set_name: str
    threads: int
    mode: Optional[str]
    repeat: int
    wall: float
    cpu: float
    load: Optional[float]


def _load_from_notes(notes: str) -> Optional[float]:
    m = re.search(r"(?:^|[;\s])load_s=([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)", notes)
    return float(m.group(1)) if m else None


def load_timings(path: Path) -> tuple:
    """(rows, n_untagged): every repetition tagged functional=V_T_dg; older rows are counted and skipped."""
    if not path.exists():
        warn(f"no timings file {path}")
        return [], 0
    rows, skipped = [], 0
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if FUNCTIONAL_TAG not in r["notes"]:
                skipped += 1
                continue
            engine = r["engine"]
            mode = None
            if engine == "charm":
                m = MODE_RE.search(r["notes"])
                if m is None:
                    warn(f"charm timing row without parallel= tag skipped: {r['notes'][:80]}")
                    continue
                mode = m.group(1)
            rows.append(TimingRow(engine, recurrence_label(engine, r["recurrence"], mode), r["set"],
                                  int(r["threads"]), mode, int(r["repeat"]),
                                  float(r["wall_s"]), float(r["cpu_s"]), _load_from_notes(r["notes"])))
    return rows, skipped


def summarize_times(rows: list) -> list:
    groups = {}
    for t in rows:
        groups.setdefault((t.engine, t.rec, t.set_name, t.threads, t.mode), []).append(t)
    out = []
    for (engine, rec, set_name, threads, mode), g in groups.items():
        walls = [t.wall for t in g]
        loads = [t.load for t in g if t.load is not None]
        out.append({
            "engine": engine, "rec": rec, "set": set_name, "threads": threads,
            "mode": mode, "parallel_mode": parallel_mode(engine, mode), "n": len(g),
            "wall_med": statistics.median(walls), "wall_min": min(walls), "wall_max": max(walls),
            "cpu_med": statistics.median(t.cpu for t in g),
            "load_med": statistics.median(loads) if loads else None,
        })
    return sorted(out, key=lambda x: (ENGINE_RANK[x["engine"]], REC_RANK[x["rec"]],
                                      SET_RANK[x["set"]], x["threads"]))


TIME_COLUMNS = ["engine", "recurrence", "set", "threads", "parallel_mode", "n_repeats",
                "wall_s_median", "wall_s_min", "wall_s_max", "wall_ms_median", "cpu_s_median",
                "load_s_median", "notes"]


def time_row_out(s: dict) -> list:
    load = "" if s["load_med"] is None else f"{s['load_med']:.3f}"
    return [s["engine"], s["rec"], s["set"], s["threads"], s["parallel_mode"], s["n"],
            f"{s['wall_med']:.4f}", f"{s['wall_min']:.4f}", f"{s['wall_max']:.4f}",
            f"{s['wall_med'] * 1e3:.1f}", f"{s['cpu_med']:.4f}", load, NOTES[s["engine"]]]


# ----------------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------------
VALUE_OUT_COLUMNS = ["engine", "recurrence", "set", "id", "lat_deg", "lon_deg", "h_m",
                     "V_m2s2", "T_m2s2", "dg_mgal"]


def write_folder_outputs(csv_dir: Path, value_rows: dict, time_rows: dict) -> list:
    """value_rows / time_rows: folder -> list of rows. Writes res_value.csv and res_time.csv per folder."""
    written = []
    for folder in FOLDERS:
        vr = value_rows.get(folder, [])
        tr = time_rows.get(folder, [])
        if not vr and not tr:
            warn(f"no data for folder {folder}: nothing written")
            continue
        out = csv_dir / folder
        out.mkdir(parents=True, exist_ok=True)
        with (out / "res_value.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(VALUE_OUT_COLUMNS)
            w.writerows(vr)
        with (out / "res_time.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(TIME_COLUMNS)
            w.writerows(tr)
        written.append(folder)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--raw", default=str(DEFAULT_RAW), help="raw outputs folder")
    ap.add_argument("--csv", default=str(DEFAULT_CSV), help="output folder (results/csv)")
    args = ap.parse_args(argv)
    raw, out = Path(args.raw), Path(args.csv)

    files = load_value_files(raw)
    chosen = select_values(files)
    timing_rows, skipped = load_timings(raw / "timings.csv")
    if skipped:
        print(f"note: {skipped} timings row(s) without {FUNCTIONAL_TAG} skipped")

    # value rows of the chosen (production) run of each configuration
    value_rows = {}
    for key in sorted(chosen, key=lambda k: (ENGINE_RANK[k[0]], REC_RANK[k[1]], SET_RANK[k[2]])):
        vf = chosen[key]
        folder = FOLDER_OF[vf.engine]
        for raw_row in vf.rows:
            r = [c.strip() for c in raw_row]
            value_rows.setdefault(folder, []).append(
                [vf.engine, vf.rec, vf.set_name, r[0], r[1], r[2], r[3], r[4], r[5], r[6]])

    # timing rows: one per (engine, recurrence, set, threads, mode)
    time_rows = {}
    summaries = summarize_times(timing_rows)
    for s in summaries:
        time_rows.setdefault(FOLDER_OF[s["engine"]], []).append(time_row_out(s))

    # configurations without timings, or timings without values, are reported
    value_keys = {k for k in chosen}
    time_keys = {(s["engine"], s["rec"], s["set"]) for s in summaries}
    for k in sorted(value_keys - time_keys):
        warn(f"values without timings: {k}")
    for k in sorted(time_keys - value_keys):
        warn(f"timings without values: {k}")

    written = write_folder_outputs(out, value_rows, time_rows)
    print(f"read {len(files)} value file(s), {len(chosen)} configuration(s), {len(summaries)} timing group(s)")
    print(f"wrote folders {', '.join(written)} under {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
