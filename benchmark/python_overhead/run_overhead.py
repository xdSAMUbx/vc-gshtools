"""run_overhead.py: what the Python wrapper costs on top of one C++ call of gshtools.compute_functionals.

Definitions, inputs, the table layout and the limits are in README.md (same folder). In short, for every case
(points N, nmax, functionals) the script measures, with the same inputs:

  (a) Python: wall time of Points.from_geodetic(numpy arrays) + compute_functionals() + reading the result arrays
      from the returned dict (the arrays are zero-copy, so this is the access cost). The three parts are also timed
      separately: py_points_s, py_compute_s, py_read_s.
  (b) C++: the same work timed inside gsh_overhead_ref (gsh_overhead_ref.cpp): Points::from_geodetic + compute_
      functionals, wall time. The library's own stage timers (FunctionalResult.timing_s) are reported too.

overhead_s = (a) - (b), overhead_pct = 100 (a) - (b) / (b). Each case has one warm-up call (phase "warmup", kept
apart: it includes the OpenMP start-up) and REPEAT timed calls (phase "timed", REPEAT = 8).

Modes:
  python run_overhead.py                      full run, writes results/overhead.csv and prints the summary
  python run_overhead.py --quick              smoke test: 2 repetitions, N = 1, 10, nmax 90 (EGM2008) only
  python run_overhead.py --correctness-only   one-time bitwise check Python == C++ (threads = 1), no timing
  python run_overhead.py --summary            prints the summary of an existing CSV (results/overhead.csv)

The script needs the gshtools package (PYTHONPATH=build/py/python in-tree, or pip install) and the C++ driver
(build target gsh_overhead_ref, default build/clang/bin/gsh_overhead_ref.exe, or --cpp-exe).
"""

import argparse
import csv
import gc
import math
import os
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import gshtools as gsh

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EGM2008 = ROOT / "data" / "cache" / "modelos" / "EGM2008" / "EGM2008.gfc"
XGM2019E = ROOT / "data" / "cache" / "modelos" / "XGM2019e" / "XGM2019e.gfc"
POINTS_SOURCE = ROOT / "benchmark" / "shs" / "data" / "points_colorado.csv"  # the shs Colorado points (509)
WORK = ROOT / "output" / "python_overhead"  # generated inputs and dumps (under output/, not in git)
RESULTS_CSV = HERE / "results" / "overhead.csv"
QUICK_CSV = WORK / "overhead_quick.csv"

N_VALUES = (1, 10, 100, 509, 10000)
NMAX_VALUES = (90, 720, 2190)
FUNCTIONALS = ("V", "VTdg")
MODEL_FOR_NMAX = {90: EGM2008, 720: EGM2008, 2190: XGM2019E}  # EGM2008 is degree 2190; XGM2019e file is 5540
REPEAT = 8
WARMUP = 1
GRID_N = 10000  # generated points: a 100 x 100 grid over the Colorado box, no poles
GRID_BOX = {"lat": (37.0, 41.0), "lon": (-109.0, -102.0), "h": 1500.0}

CSV_COLUMNS = [
    "case_id", "n_points", "nmax", "functionals", "model", "phase", "rep", "threads",
    "py_total_s", "py_points_s", "py_compute_s", "py_read_s",
    "cpp_total_s", "cpp_points_s", "cpp_compute_s", "cpp_stage_sum_s",
    "overhead_s", "overhead_pct",
]


def default_cpp_exe():
    name = "gsh_overhead_ref.exe" if os.name == "nt" else "gsh_overhead_ref"
    return ROOT / "build" / "clang" / "bin" / name


@dataclass(frozen=True)
class Case:
    n: int
    nmax: int
    functionals: str

    @property
    def case_id(self):
        return f"N{self.n}_nmax{self.nmax}_{self.functionals}"


def outputs_of(name):
    if name == "V":
        return [gsh.Functional.V]
    return [gsh.Functional.V, gsh.Functional.T, gsh.Functional.dg]


def make_request(nmax, name, threads):
    req = gsh.FunctionalRequest()
    req.outputs = outputs_of(name)
    req.nmin = 0
    req.nmax = nmax
    req.ellipsoid = gsh.Ellipsoid.GRS80()
    opt = gsh.SynthOptions()
    opt.threads = threads
    req.synth = opt
    return req


def points_path(n):
    return WORK / "points" / f"points_n{n}.csv"


def prepare_points(ns):
    """Writes one points CSV per N (header id,lat_deg,lon_deg,h_m). Both sides read these same files."""
    (WORK / "points").mkdir(parents=True, exist_ok=True)
    source_lines = POINTS_SOURCE.read_text(encoding="utf-8").splitlines()
    for n in ns:
        out = points_path(n)
        if n == GRID_N:
            side = math.isqrt(n)
            if side * side != n:
                raise ValueError("GRID_N must be a square")
            lat = np.linspace(*GRID_BOX["lat"], side)
            lon = np.linspace(*GRID_BOX["lon"], side)
            LA, LO = np.meshgrid(lat, lon, indexing="ij")
            lines = ["id,lat_deg,lon_deg,h_m"]
            for i, (la, lo) in enumerate(zip(LA.ravel(), LO.ravel())):
                lines.append(f"{i + 1},{la:.12f},{lo:.12f},{GRID_BOX['h']:.1f}")
        else:
            if n > len(source_lines) - 1:
                raise ValueError(f"{POINTS_SOURCE} has fewer than {n} points")
            lines = source_lines[: n + 1]
        out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_points(n):
    """Latitude, longitude (degrees) and height (metres) as numpy arrays, read with gshtools.read_table."""
    tab = gsh.read_table(str(points_path(n)), sep=",", header_lines=1, cols=[1, 2, 3])
    if len(tab.skipped_lines):
        raise RuntimeError(f"non-numeric line in {points_path(n)}")
    lat, lon, h = tab.columns
    return lat, lon, h


def python_call(model, lat, lon, h, ell, req):
    """One timed call. Returns (total, points, compute, read) in seconds."""
    gc.collect()  # outside the timed block: keep the collector from running inside it
    t0 = time.perf_counter()
    pts = gsh.Points.from_geodetic(lat, lon, h, ell.a, ell.e2)
    t1 = time.perf_counter()
    res = gsh.compute_functionals(model, pts, req)
    t2 = time.perf_counter()
    _ = [res.values[f] for f in req.outputs]  # the result arrays handed to Python
    t3 = time.perf_counter()
    return t3 - t0, t1 - t0, t2 - t1, t3 - t2


def run_cpp(exe, model_path, nmax, jobs, threads, warmup, repeat, dump_dir=None):
    """Runs gsh_overhead_ref once for a list of jobs (name, points path). Returns {job index: [records]}."""
    cmd = [str(exe), "--gfc", str(model_path), "--nmax", str(nmax), "--threads", str(threads),
           "--warmup", str(warmup), "--repeat", str(repeat)]
    for name, path in jobs:
        cmd += ["--job", name, str(path)]
    if dump_dir is not None:
        dump_dir.mkdir(parents=True, exist_ok=True)
        cmd += ["--dump-dir", str(dump_dir)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"gsh_overhead_ref failed ({proc.returncode}): {proc.stderr.strip()}")
    records = {}
    for line in proc.stdout.splitlines():
        if not line.startswith("RUN,"):
            continue
        _, job, phase, rep, points_s, compute_s, stage_s = line.split(",")
        records.setdefault(int(job), []).append({
            "phase": phase, "rep": int(rep), "points_s": float(points_s),
            "compute_s": float(compute_s), "stage_sum_s": float(stage_s),
        })
    return records


def group_by_nmax(cases):
    groups = {}
    for case in cases:
        groups.setdefault(case.nmax, []).append(case)
    return dict(sorted(groups.items()))


def build_cases(quick):
    if quick:
        return [Case(n, 90, f) for n in (1, 10) for f in FUNCTIONALS]
    return [Case(n, nmax, f) for n in N_VALUES for nmax in NMAX_VALUES for f in FUNCTIONALS]


def timing_run(cases, repeat, warmup, threads, exe):
    """Full measurement. Per nmax group: Python calls first, then the C++ driver. Returns CSV rows."""
    ell = gsh.Ellipsoid.GRS80()
    rows = []
    for nmax, group in group_by_nmax(cases).items():
        model_path = MODEL_FOR_NMAX[nmax]
        print(f"[load] {model_path.name} nmax={nmax} (not timed)", file=sys.stderr)
        model = gsh.load_gfc(str(model_path), nmax)

        py_calls = {}
        for case in group:
            lat, lon, h = read_points(case.n)
            req = make_request(nmax, case.functionals, threads)
            py_calls[case] = [python_call(model, lat, lon, h, ell, req) for _ in range(warmup + repeat)]
            print(f"[py ] {case.case_id} done", file=sys.stderr)
        del model

        jobs = [(case.functionals, points_path(case.n)) for case in group]
        cpp = run_cpp(exe, model_path, nmax, jobs, threads, warmup, repeat)
        print(f"[cpp] nmax={nmax} done", file=sys.stderr)

        for j, case in enumerate(group):
            recs = cpp.get(j, [])
            calls = py_calls[case]
            if len(recs) != len(calls):
                raise RuntimeError(f"{case.case_id}: {len(recs)} C++ calls, {len(calls)} Python calls")
            for (total, p_s, c_s, r_s), rec in zip(calls, recs):
                cpp_total = rec["points_s"] + rec["compute_s"]
                overhead = total - cpp_total
                rows.append({
                    "case_id": case.case_id, "n_points": case.n, "nmax": nmax, "functionals": case.functionals,
                    "model": model_path.name, "phase": rec["phase"], "rep": rec["rep"], "threads": threads,
                    "py_total_s": total, "py_points_s": p_s, "py_compute_s": c_s, "py_read_s": r_s,
                    "cpp_total_s": cpp_total, "cpp_points_s": rec["points_s"], "cpp_compute_s": rec["compute_s"],
                    "cpp_stage_sum_s": rec["stage_sum_s"],
                    "overhead_s": overhead, "overhead_pct": 100.0 * overhead / cpp_total,
                })
    return rows


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        for row in rows:
            out = {}
            for k in CSV_COLUMNS:
                v = row[k]
                out[k] = f"{v:.9e}" if isinstance(v, float) else v
            w.writerow(out)


def correctness_run(cases, exe, threads=1):
    """Python result == C++ result, bit for bit, for every case (threads = 1). Returns True when all match."""
    ell = gsh.Ellipsoid.GRS80()
    all_ok = True
    print("| case | outputs | points | bitwise equal |")
    print("|---|---|---|---|")
    for nmax, group in group_by_nmax(cases).items():
        model = gsh.load_gfc(str(MODEL_FOR_NMAX[nmax]), nmax)
        py_values = {}
        for case in group:
            lat, lon, h = read_points(case.n)
            req = make_request(nmax, case.functionals, threads)
            pts = gsh.Points.from_geodetic(lat, lon, h, ell.a, ell.e2)
            res = gsh.compute_functionals(model, pts, req)
            py_values[case] = [np.ascontiguousarray(res.values[f], dtype=np.float64) for f in req.outputs]
        del model
        dump_dir = WORK / "dump" / f"nmax{nmax}"
        jobs = [(case.functionals, points_path(case.n)) for case in group]
        run_cpp(exe, MODEL_FOR_NMAX[nmax], nmax, jobs, threads, warmup=1, repeat=0, dump_dir=dump_dir)
        for j, case in enumerate(group):
            outs = outputs_of(case.functionals)
            dump = np.fromfile(dump_dir / f"job{j}_threads{threads}_{case.functionals}.bin", dtype="<f8")
            cpp_values = dump.reshape(len(outs), case.n)
            equal = all(np.array_equal(py.view(np.uint64), cpp_values[k].view(np.uint64))
                        for k, py in enumerate(py_values[case]))
            all_ok &= equal
            names = ",".join(gsh.functional_name(f) for f in outs)
            print(f"| {case.case_id} | {names} | {case.n} | {'yes' if equal else 'NO'} |")
            sys.stdout.flush()
    return all_ok


def _median_min_max(values):
    return statistics.median(values), min(values), max(values)


def summarize(path):
    """Markdown summary of a CSV written by run_overhead.py (median and min/max over the timed repetitions)."""
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    order = []
    for r in rows:
        if r["case_id"] not in order:
            order.append(r["case_id"])
    lines = [
        f"# Python overhead of gshtools.compute_functionals ({path.name})",
        "",
        "Overhead = Python call - C++ call (same inputs, same work). Positive: Python is slower. Medians over the timed",
        "repetitions; [min, max] in brackets. The warm-up call is listed apart (it includes the OpenMP start-up).",
        "",
        "| case | N | nmax | outputs | C++ median (s) | Python median (s) | overhead median (s) [min, max] | "
        "overhead median (%) [min, max] | Points, Python minus C++ (ms) | compute call, Python minus C++ (ms) | "
        "read arrays, Python (ms) |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    warm_lines = []
    for cid in order:
        timed = [r for r in rows if r["case_id"] == cid and r["phase"] == "timed"]
        warm = [r for r in rows if r["case_id"] == cid and r["phase"] == "warmup"]
        if not timed:
            continue
        g = lambda key: [float(r[key]) for r in timed]  # noqa: E731
        cpp_med = statistics.median(g("cpp_total_s"))
        py_med = statistics.median(g("py_total_s"))
        ov_med, ov_min, ov_max = _median_min_max(g("overhead_s"))
        pct_med, pct_min, pct_max = _median_min_max(g("overhead_pct"))
        pts_diff = statistics.median([a - b for a, b in zip(g("py_points_s"), g("cpp_points_s"))]) * 1e3
        cmp_diff = statistics.median([a - b for a, b in zip(g("py_compute_s"), g("cpp_compute_s"))]) * 1e3
        read_med = statistics.median(g("py_read_s")) * 1e3
        first = timed[0]
        lines.append(
            f"| {cid} | {first['n_points']} | {first['nmax']} | {first['functionals']} | {cpp_med:.6g} | {py_med:.6g} | "
            f"{ov_med:.3e} [{ov_min:.3e}, {ov_max:.3e}] | {pct_med:.3f} [{pct_min:.3f}, {pct_max:.3f}] | "
            f"{pts_diff:.4f} | {cmp_diff:.4f} | {read_med:.4f} |")
        for r in warm:
            warm_lines.append(
                f"| {cid} | {float(r['cpp_total_s']):.6g} | {float(r['py_total_s']):.6g} | "
                f"{float(r['overhead_s']):.3e} | {float(r['overhead_pct']):.3f} |")
    lines += ["", "Warm-up call (first call of the process for this case, not in the statistics above):", "",
              "| case | C++ (s) | Python (s) | overhead (s) | overhead (%) |", "|---|---|---|---|---|"] + warm_lines
    reps = sum(1 for r in rows if r["case_id"] == order[0] and r["phase"] == "timed") if order else 0
    lines += ["", f"Timed repetitions per case: {reps}."]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Overhead of the gshtools Python wrapper over one C++ compute_functionals call (see README.md).")
    parser.add_argument("--quick", action="store_true",
                        help="smoke test only: 2 repetitions, N = 1 and 10, nmax 90; writes output/python_overhead/")
    parser.add_argument("--repeat", type=int, default=REPEAT, help=f"timed calls per case (default {REPEAT})")
    parser.add_argument("--warmup", type=int, default=WARMUP,
                        help=f"warm-up calls per case, recorded apart (default {WARMUP})")
    parser.add_argument("--threads", type=int, default=0,
                        help="OpenMP threads for both sides; 0 = OpenMP default (default 0)")
    parser.add_argument("--cpp-exe", type=Path, default=None,
                        help="gsh_overhead_ref executable (default build/clang/bin/gsh_overhead_ref.exe)")
    parser.add_argument("--csv", type=Path, default=None,
                        help="CSV to write (run) or read (--summary); default results/overhead.csv, "
                             "or output/python_overhead/overhead_quick.csv with --quick")
    parser.add_argument("--summary", action="store_true", help="print the summary of --csv and exit")
    parser.add_argument("--correctness-only", action="store_true",
                        help="one-time bitwise check Python == C++ with threads = 1 (no timing, nothing written "
                             "to results/)")
    args = parser.parse_args()

    csv_path = args.csv or (QUICK_CSV if args.quick else RESULTS_CSV)
    if args.summary:
        print(summarize(csv_path))
        return 0
    if args.repeat < 1 or args.warmup < 0 or args.threads < 0:
        parser.error("--repeat must be >= 1, --warmup >= 0, --threads >= 0")

    exe = args.cpp_exe or default_cpp_exe()
    if not exe.exists():
        parser.error(f"C++ driver not found: {exe}. Build the target gsh_overhead_ref (see README.md) or pass --cpp-exe")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    repeat = 2 if args.quick else args.repeat
    cases = build_cases(args.quick)
    prepare_points(sorted({c.n for c in cases}))
    for case in cases:
        if not MODEL_FOR_NMAX[case.nmax].exists():
            parser.error(f"model file missing: {MODEL_FOR_NMAX[case.nmax]}")

    if args.correctness_only:
        ok = correctness_run(cases, exe)
        print("\nbitwise check:", "PASS" if ok else "FAIL")
        return 0 if ok else 1

    print("Run on an otherwise idle machine. Progress goes to stderr.", file=sys.stderr)
    rows = timing_run(cases, repeat, args.warmup, args.threads, exe)
    write_csv(csv_path, rows)
    print(f"rows written: {len(rows)} -> {csv_path}", file=sys.stderr)
    print(summarize(csv_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
