"""Compare the spherical-harmonic synthesis results, write the comparison CSVs, tables and figures.

Reads the raw outputs of code/run_all.ps1 (results/raw, see consolidate.py for the naming) and writes:

  * results/csv/<folder>/comp_value.csv for the folders charm, shtools and gshtools
      columns pair,recurrence_a,recurrence_b,set,quantity,unit,n,mean_d,rms_d,min_d,max_d,max_abs_d
      d = first - second at each point, matched on the id (coordinates are checked too).
      quantity V, T (m^2/s^2) or dg (mGal). The comparisons vs GrafLab (era for the single-recurrence
      engines, same recurrence for gsh sfc/msfc/era) are in the folder of the engine being compared.
  * results/raw/summary_tables.md: markdown tables of the V, T and dg differences, the run times (median,
      min, max over the repetitions) and the thread-consistency check. Pasted by hand into Summary.md.
  * PNG figures in img/ (dpi 150, light background, titles and axis units; numbers printed on a figure are
      computed from the data):
        accuracy_max_abs_diff_colorado.png     max|d| against GrafLab, Colorado set, V, T and dg
        timing_bars_colorado.png               wall time per engine (log scale), Colorado, the largest thread count
        latitude_sweep_V_gsh_vs_graflab.png    |d| of V per latitude, gsh sfc/msfc/era against GrafLab (latsweep)
        speed_1_vs_max_threads_colorado.png    wall time with 1 thread and with the largest thread count (Colorado)

No REPORT*.md file is written. The reference of the accuracy comparisons is GrafLab era, except where a
recurrence of gsh is compared with the same GrafLab recurrence (see COMPARISONS).

Usage (repository root, project venv):
    .venv\\Scripts\\python.exe benchmark\\shs\\code\\python\\compare.py
    .venv\\Scripts\\python.exe benchmark\\shs\\code\\python\\compare.py --raw <dir> --csv <dir> --img <dir>
"""
from __future__ import annotations

import argparse
import datetime as dt
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from consolidate import (DEFAULT_CSV, DEFAULT_RAW, ENGINE_RANK, REC_RANK, SET_RANK, SETS,  # noqa: E402
                         load_timings, load_value_files, select_values, summarize_times, warn)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SHS_DIR = Path(__file__).resolve().parents[2]
DEFAULT_IMG = SHS_DIR / "img"
ID_COORD_TOL_DEG = 1e-5
ID_HEIGHT_TOL_M = 1e-2
CPU_SINGLE_RATIO = 1.2          # threads=1 rows with cpu/wall above this are flagged (threads not honoured)
QUANTITIES = [("V_m2s2", "V", "m^2/s^2"), ("T_m2s2", "T", "m^2/s^2"), ("dg_mgal", "dg", "mGal")]
COMP_COLUMNS = ["pair", "recurrence_a", "recurrence_b", "set", "quantity", "unit", "n",
                "mean_d", "rms_d", "min_d", "max_d", "max_abs_d"]
ENGINE_COLOR = {"gsh": "#1f77b4", "graflab": "#333333", "pyshtools": "#ff7f0e", "shtools": "#d62728",
                "pyharm": "#2ca02c", "charm": "#9467bd"}

# folder -> list of (engine_a, rec_a, engine_b, rec_b); d = a - b
COMPARISONS = {
    "charm": [
        ("pyharm", "native", "charm", "native-charm"),
        ("pyharm", "native", "graflab", "era"),
        ("charm", "native-charm", "graflab", "era"),
    ],
    "shtools": [
        ("pyshtools", "native", "shtools", "native"),
        ("pyshtools", "native", "graflab", "era"),
        ("shtools", "native", "graflab", "era"),
    ],
    "gshtools": [
        ("gsh", "sfc", "graflab", "sfc"),
        ("gsh", "msfc", "graflab", "msfc"),
        ("gsh", "era", "graflab", "era"),
        ("gsh", "fast", "graflab", "era"),
    ],
}
FOLDER_OF_ENGINE = {"gsh": "gshtools", "graflab": "graflab", "pyharm": "charm", "charm": "charm",
                    "pyshtools": "shtools", "shtools": "shtools"}


def md_table(headers: list, rows: list) -> list:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines + [""]


def fe(x: float) -> str:
    return f"{x:.4e}" if math.isfinite(x) else "nan"


def ff(x: float, nd: int = 4) -> str:
    return f"{x:.{nd}f}" if math.isfinite(x) else "nan"


def field_of(vf, column: str) -> np.ndarray:
    return np.asarray(vf.column(column), dtype=float)


def align(a, b):
    """Index arrays (ia, ib) of the common ids, or None when the coordinates of a common id differ."""
    pos = {pid: k for k, pid in enumerate(b.ids)}
    ia = [k for k, pid in enumerate(a.ids) if pid in pos]
    if not ia:
        warn(f"{a.path.name} and {b.path.name} share no point ids")
        return None
    ib = [pos[a.ids[k]] for k in ia]
    ia, ib = np.array(ia), np.array(ib)
    la, lb = np.asarray(a.lat)[ia], np.asarray(b.lat)[ib]
    oa, ob = np.asarray(a.lon)[ia], np.asarray(b.lon)[ib]
    ha, hb = np.asarray(a.h)[ia], np.asarray(b.h)[ib]
    if not (np.all(np.abs(la - lb) <= ID_COORD_TOL_DEG) and np.all(np.abs(oa - ob) <= ID_COORD_TOL_DEG)
            and np.all(np.abs(ha - hb) <= ID_HEIGHT_TOL_M)):
        warn(f"{a.path.name} and {b.path.name}: same ids but different coordinates: not compared")
        return None
    if ia.size != len(a.ids) or ia.size != len(b.ids):
        warn(f"{a.path.name} vs {b.path.name}: only {ia.size} common points")
    return ia, ib


def stats_of(d: np.ndarray) -> dict:
    finite = d[np.isfinite(d)]
    if finite.size != d.size:
        warn(f"{d.size - finite.size} non-finite difference(s) excluded")
    if finite.size == 0:
        return {"n": 0, "mean": math.nan, "rms": math.nan, "min": math.nan, "max": math.nan, "maxabs": math.nan}
    return {"n": int(finite.size), "mean": float(finite.mean()), "rms": float(np.sqrt(np.mean(finite ** 2))),
            "min": float(finite.min()), "max": float(finite.max()), "maxabs": float(np.abs(finite).max())}


def compare_values(a, b, column: str):
    """(stats, d, reference, latitude) of a - b over the common points, or None."""
    al = align(a, b)
    if al is None:
        return None
    ia, ib = al
    va = field_of(a, column)[ia]
    vb = field_of(b, column)[ib]
    d = va - vb
    return stats_of(d), d, vb, np.asarray(a.lat)[ia]


# ----------------------------------------------------------------------------------------
# comp_value.csv, summary tables
# ----------------------------------------------------------------------------------------
def build_comparisons(chosen: dict) -> list:
    out = []
    for folder, pairs in COMPARISONS.items():
        for set_name in SETS:
            for ea, ra, eb, rb in pairs:
                a = chosen.get((ea, ra, set_name))
                b = chosen.get((eb, rb, set_name))
                if a is None or b is None:
                    warn(f"comparison {ea} {ra} vs {eb} {rb} ({set_name}): a file is missing, skipped")
                    continue
                for column, short, unit in QUANTITIES:
                    res = compare_values(a, b, column)
                    if res is None:
                        continue
                    st, d, ref, lat = res
                    out.append({"folder": folder, "pair": f"{ea}_vs_{eb}", "ra": ra, "rb": rb, "ea": ea,
                                "eb": eb, "set": set_name, "column": column, "short": short, "unit": unit,
                                "st": st, "d": d, "ref": ref, "lat": lat})
    return out


def write_comp_csvs(comps: list, csv_dir: Path) -> None:
    for folder in COMPARISONS:
        rows = [c for c in comps if c["folder"] == folder]
        out = csv_dir / folder
        out.mkdir(parents=True, exist_ok=True)
        with (out / "comp_value.csv").open("w", encoding="utf-8", newline="") as fh:
            fh.write(",".join(COMP_COLUMNS) + "\n")
            for c in rows:
                st = c["st"]
                vals = [c["pair"], c["ra"], c["rb"], c["set"], c["short"], c["unit"], st["n"],
                        f"{st['mean']:.10e}", f"{st['rms']:.10e}", f"{st['min']:.10e}",
                        f"{st['max']:.10e}", f"{st['maxabs']:.10e}"]
                fh.write(",".join(str(v) for v in vals) + "\n")
        print(f"wrote {out / 'comp_value.csv'} ({len(rows)} rows)")


def thread_checks(files: list, chosen: dict) -> list:
    """Each configuration run with other thread counts, against the production run (largest thread count)."""
    rows = []
    for vf in files:
        ref = chosen.get(vf.key)
        if ref is None or ref is vf or vf.threads is None or ref.threads is None:
            continue
        if vf.threads == ref.threads:
            continue
        al = align(vf, ref)
        if al is None:
            continue
        ia, ib = al
        maxes = []
        for column, _, _ in QUANTITIES:
            d = field_of(vf, column)[ia] - field_of(ref, column)[ib]
            maxes.append(float(np.nanmax(np.abs(d))) if d.size else math.nan)
        rows.append((vf.engine, vf.rec, vf.set_name, vf.threads, ref.threads, maxes))
    return sorted(rows, key=lambda r: (ENGINE_RANK[r[0]], REC_RANK[r[1]], SET_RANK[r[2]], r[3]))


def build_summary_md(comps, summaries, checks, timing_rows, raw: Path) -> list:
    lines = ["# Summary tables", "",
             f"Generated by `code/python/compare.py` on {dt.datetime.now().isoformat(timespec='seconds')} "
             f"from `{raw}`. Do not edit by hand: rerun the script. Differences are d = first - second at each "
             "point, matched on the id. Times: synthesis wall time (load excluded), median, min and max over the "
             "repetitions.", ""]
    for column, short, unit in QUANTITIES:
        lines += [f"## Differences of {short} [{unit}]", ""]
        rows = []
        for c in comps:
            if c["column"] != column:
                continue
            st = c["st"]
            rows.append([c["pair"], c["ra"], c["rb"], c["set"], st["n"], fe(st["mean"]), fe(st["rms"]),
                         fe(st["min"]), fe(st["max"]), fe(st["maxabs"])])
        lines += md_table(["pair", "recurrence a", "recurrence b", "set", "n", "mean(d)", "RMS(d)",
                           "min(d)", "max(d)", "max abs(d)"], rows)
    lines += ["## Run times [s] (synthesis wall time, one row per engine, recurrence, set and thread count)", ""]
    rows = []
    for s in summaries:
        rows.append([s["engine"], s["rec"], s["set"], s["threads"], s["parallel_mode"], s["n"],
                     ff(s["wall_med"]), ff(s["wall_min"]), ff(s["wall_max"]), ff(s["cpu_med"]),
                     "" if s["load_med"] is None else ff(s["load_med"], 3)])
    lines += md_table(["engine", "recurrence", "set", "threads", "parallel_mode", "repeats", "wall median",
                       "wall min", "wall max", "cpu median", "load median"], rows)
    lines += ["## Thread check (same configuration, other thread count; d = run - production run, largest threads)", ""]
    rows = [[e, r, s, t, ref, fe(m[0]), fe(m[1]), fe(m[2])] for e, r, s, t, ref, m in checks]
    lines += md_table(["engine", "recurrence", "set", "threads", "production threads",
                       "max abs(d) V [m^2/s^2]", "max abs(d) T [m^2/s^2]", "max abs(d) dg [mGal]"], rows)
    lines += ["## Warnings", ""]
    worst = {}   # one warning per configuration: the largest cpu/wall of its repetitions
    for t in timing_rows:
        if t.threads == 1 and t.cpu > 0 and t.wall > 0:
            key = (t.engine, t.rec, t.set_name)
            worst[key] = max(worst.get(key, 0.0), t.cpu / t.wall)
    warnings = [f"{e} {r} ({st}) with 1 thread: cpu/wall = {ratio:.2f} (threads not honoured?)"
                for (e, r, st), ratio in sorted(worst.items()) if ratio > CPU_SINGLE_RATIO]
    lines += [f"- {w}" for w in warnings] if warnings else ["No warnings."]
    return lines + [""]


# ----------------------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------------------
def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {path}")


def fig_accuracy(comps: list, img: Path) -> None:
    items = [c for c in comps if c["set"] == "colorado" and c["eb"] == "graflab"]
    if not items:
        warn("no Colorado comparison against GrafLab: accuracy figure skipped")
        return
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=False)
    for ax, (column, short, unit) in zip(axes, QUANTITIES):
        sel = [c for c in items if c["column"] == column]
        labels = [f"{c['ea']} {c['ra']} vs GrafLab {c['rb']}" + ("  (max|d| = 0)" if c["st"]["maxabs"] == 0 else "")
                  for c in sel]
        vals = [c["st"]["maxabs"] for c in sel]
        colors = [ENGINE_COLOR[c["ea"]] for c in sel]
        y = np.arange(len(sel))
        ax.barh(y, [v if v > 0 else np.nan for v in vals], color=colors)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel(f"max |d| of {short} [{unit}] (log scale)")
        ax.set_title(f"{short}: maximum absolute difference against GrafLab (Colorado, {len(sel)} comparisons)",
                     fontsize=10)
        ax.grid(axis="x", alpha=0.3)
        for yi, v in zip(y, vals):
            if v > 0:
                ax.text(v, yi, f" {v:.2e}", va="center", fontsize=7)
    fig.suptitle("Accuracy of the engines against GrafLab (matched by point id)", fontsize=11)
    _save(fig, img / "accuracy_max_abs_diff_colorado.png")


def _times(summaries: list, set_name: str) -> list:
    return [s for s in summaries if s["set"] == set_name]


def fig_timing(summaries: list, img: Path) -> None:
    rows = _times(summaries, "colorado")
    if not rows:
        warn("no Colorado timings: timing figure skipped")
        return
    # one row per engine configuration: its largest thread count (serial engines have only 1)
    best = {}
    for s in rows:
        k = (s["engine"], s["rec"], s["parallel_mode"])
        if k not in best or s["threads"] > best[k]["threads"]:
            best[k] = s
    sel = sorted(best.values(), key=lambda s: s["wall_med"])
    fig, ax = plt.subplots(figsize=(9, 0.5 + 0.42 * len(sel)))
    y = np.arange(len(sel))
    lo = [s["wall_med"] - s["wall_min"] for s in sel]
    hi = [s["wall_max"] - s["wall_med"] for s in sel]
    ax.barh(y, [s["wall_med"] for s in sel], xerr=[lo, hi], color=[ENGINE_COLOR[s["engine"]] for s in sel],
            capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s['engine']} {s['rec']} [{s['parallel_mode']}], {s['threads']} thr" for s in sel],
                       fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("wall time of the synthesis, median of the repetitions [s] (log scale; bars: min to max)")
    ax.set_title("Run time per engine configuration, Colorado set (nmax 2190, V, T and dg)\n"
                 "largest thread count of each configuration (GrafLab: MATLAB default threads)", fontsize=10)
    ax.grid(axis="x", alpha=0.3)
    for yi, s in zip(y, sel):
        ax.text(s["wall_max"] * 1.08, yi, f"{s['wall_med']:.3g} s", va="center", fontsize=7)
    _save(fig, img / "timing_bars_colorado.png")


def fig_latitude(chosen: dict, img: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    drawn = 0
    for rec in ("sfc", "msfc", "era"):
        a = chosen.get(("gsh", rec, "latsweep"))
        b = chosen.get(("graflab", rec, "latsweep"))
        if a is None or b is None:
            warn(f"latitude sweep: gsh {rec} or GrafLab {rec} missing, not drawn")
            continue
        res = compare_values(a, b, "V_m2s2")
        if res is None:
            continue
        st, d, _, lat = res
        order = np.argsort(lat)
        absd = np.abs(d)[order]
        lat_s = lat[order]
        pos = absd > 0
        ax.semilogy(lat_s[pos], absd[pos], "o-", ms=3, lw=0.8,
                    label=f"gsh {rec} - GrafLab {rec} (max |d| = {st['maxabs']:.2e} m^2/s^2)")
        drawn += 1
    if not drawn:
        plt.close(fig)
        return
    ax.set_xlabel("geodetic latitude [deg] (latsweep meridian, nmax 2190)")
    ax.set_ylabel("|d| of V [m^2/s^2] (log scale)")
    ax.set_title("Latitude sweep: gsh against GrafLab, potential V (exact zeros not drawn)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    _save(fig, img / "latitude_sweep_V_gsh_vs_graflab.png")


def fig_speed(summaries: list, img: Path) -> None:
    rows = [s for s in _times(summaries, "colorado") if s["engine"] != "graflab"]
    if not rows:
        warn("no Colorado timings: speed figure skipped")
        return
    top = max(s["threads"] for s in rows)
    groups = {}
    for s in rows:
        groups.setdefault((s["engine"], s["rec"], s["parallel_mode"]), {})[s["threads"]] = s
    pairs = [(k, v[1], v[top]) for k, v in groups.items() if 1 in v and top in v]
    if not pairs:
        warn("no configuration with both 1 and the largest thread count: speed figure skipped")
        return
    pairs.sort(key=lambda p: p[0])
    fig, ax = plt.subplots(figsize=(9, 0.6 + 0.6 * len(pairs)))
    y = np.arange(len(pairs))
    h = 0.38
    t1 = [p[1]["wall_med"] for p in pairs]
    tm = [p[2]["wall_med"] for p in pairs]
    ax.barh(y - h / 2, t1, height=h, color="#9ecae1", label="1 thread")
    ax.barh(y + h / 2, tm, height=h, color="#3182bd", label=f"{top} threads")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{k[0]} {k[1]} [{k[2]}]" for k, _, _ in pairs], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("wall time of the synthesis, median [s] (log scale)")
    ax.set_title("Single thread against the largest thread count, Colorado set", fontsize=10)
    for yi, a, b in zip(y, t1, tm):
        ax.text(max(a, b), yi, f"  speed-up x{a / b:.1f}", va="center", fontsize=7)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    _save(fig, img / "speed_1_vs_max_threads_colorado.png")


# ----------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--raw", default=str(DEFAULT_RAW))
    ap.add_argument("--csv", default=str(DEFAULT_CSV))
    ap.add_argument("--img", default=str(DEFAULT_IMG))
    args = ap.parse_args(argv)
    raw, csv_dir, img = Path(args.raw), Path(args.csv), Path(args.img)

    files = load_value_files(raw)
    chosen = select_values(files)
    timing_rows, _ = load_timings(raw / "timings.csv")
    summaries = summarize_times(timing_rows)

    comps = build_comparisons(chosen)
    write_comp_csvs(comps, csv_dir)
    checks = thread_checks(files, chosen)

    md = build_summary_md(comps, summaries, checks, timing_rows, raw)
    (raw / "summary_tables.md").write_text("\n".join(md), encoding="utf-8")
    print(f"wrote {raw / 'summary_tables.md'}")

    fig_accuracy(comps, img)
    fig_timing(summaries, img)
    fig_latitude(chosen, img)
    fig_speed(summaries, img)
    return 0


if __name__ == "__main__":
    sys.exit(main())
