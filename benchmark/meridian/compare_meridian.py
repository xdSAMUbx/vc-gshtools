"""Compare the gsh meridian CSV against the GrafLab meridian CSV (see benchmark/meridian/README.md).

Usage (from the repository root):
    .venv\\Scripts\\python.exe benchmark\\meridian\\compare_meridian.py

Inputs (default in benchmark/meridian/results/):
    gsh_meridian.csv         written by gsh_meridian.exe, SI units
    graflab_era_meridian.csv written by matlab/run_graflab_meridian.m, GrafLab units (mGal, Eotvos)

Outputs:
    comparison.csv           one row per functional with the difference metrics (SI)
    comparison.md            the same table as markdown, with the conventions and the diagnostics

Conventions (fixed before running; not tuned to the outcome):
    d = gsh - GrafLab at every point, both finite. Points are matched by position (both files list lat -90..+90).
    max_abs        max |d| over the points
    scale          max |GrafLab| over the points (the functional's own magnitude)
    max_abs/scale  the normalised difference used for the success criterion (1e-12 target, rounding level)
    max_rel_pt     max |d| / |GrafLab| over the points where |GrafLab| > 0 (ill-defined near zeros, so secondary)
    Units of the GrafLab side are converted to SI: mGal -> m s^-2 (x 1e-5), Eotvos -> s^-2 (x 1e-9).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
GSH_CSV = RESULTS / "gsh_meridian.csv"
GL_CSV = RESULTS / "graflab_era_meridian.csv"
OUT_CSV = RESULTS / "comparison.csv"
OUT_MD = RESULTS / "comparison.md"
MODES = [("era", "graflab_era_meridian.csv", 3), ("msfc", "graflab_msfc_meridian.csv", 2),
         ("sfc", "graflab_sfc_meridian.csv", 1)]  # GrafLab fnALFs: 3 = ERA, 2 = modified forward column, 1 = standard
OUT_MODES_CSV = RESULTS / "comparison_modes.csv"
OUT_MODES_MD = RESULTS / "comparison_modes.md"

TARGET = 1e-12  # success criterion stated for the benchmark (rounding level, 1e-12..1e-14)

MGAL_TO_SI = 1e-5
EOTVOS_TO_SI = 1e-9


def read_table(path: Path) -> dict[str, np.ndarray]:
    """Read a CSV with a header line into a dict of float arrays (the tokens 'nan' and 'NaN' are accepted)."""
    with path.open() as f:
        header = f.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1, dtype=float)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {name: data[:, i] for i, name in enumerate(header)}


def metrics(gsh: np.ndarray, ref: np.ndarray) -> dict[str, float]:
    ok = np.isfinite(gsh) & np.isfinite(ref)
    n_ok = int(ok.sum())
    if n_ok == 0:
        return {"n": 0, "n_nonfinite": int(gsh.size), "max_abs": np.nan, "scale": np.nan,
                "max_abs_over_scale": np.nan, "max_rel_pt": np.nan, "rms": np.nan}
    g, r = gsh[ok], ref[ok]
    d = g - r
    scale = float(np.max(np.abs(r)))
    nz = np.abs(r) > 0.0
    rel_pt = np.abs(d[nz]) / np.abs(r[nz]) if nz.any() else np.array([np.nan])
    max_abs = float(np.max(np.abs(d)))
    return {
        "n": n_ok,
        "n_nonfinite": int(gsh.size - n_ok),
        "max_abs": max_abs,
        "scale": scale,
        "max_abs_over_scale": max_abs / scale if scale > 0 else np.nan,
        "max_rel_pt": float(np.nanmax(rel_pt)),
        "rms": float(np.sqrt(np.mean(d * d))),
    }


def lat_of_max(lat: np.ndarray, gsh: np.ndarray, ref: np.ndarray) -> float:
    ok = np.isfinite(gsh) & np.isfinite(ref)
    d = np.where(ok, np.abs(gsh - ref), -1.0)
    return float(lat[int(np.argmax(d))])


def modes_report(gsh: dict[str, np.ndarray], pairs: list) -> None:
    """All GrafLab Legendre modes (era, msfc, sfc): gsh against each, and each mode against era.

    The mode-against-era block measures the spread of GrafLab itself. A gsh-GrafLab difference of the same size as
    that spread is not evidence against gsh. Same metric and same criterion as the main table; nothing is tuned.
    """
    tables = {}
    for name, fname, _ in MODES:
        path = RESULTS / fname
        if path.exists():
            tables[name] = read_table(path)
    if "era" not in tables or len(tables) < 2:
        return
    lat = gsh["lat_deg"]
    rows = []
    for label, gcol, rcol, conv, unit, _ in pairs:
        row = {"functional": label, "unit": unit}
        for name, t in tables.items():
            m = metrics(gsh[gcol], t[rcol] * conv)
            row[f"gsh_vs_{name}"] = m["max_abs_over_scale"]
            row[f"gsh_vs_{name}_lat"] = lat_of_max(lat, gsh[gcol], t[rcol] * conv)
        for name, t in tables.items():
            if name == "era":
                continue
            # Both sides are GrafLab: compare on the points where both are finite.
            m = metrics(t[rcol] * conv, tables["era"][rcol] * conv)
            row[f"{name}_vs_era"] = m["max_abs_over_scale"]
        rows.append(row)

    names = list(tables)
    others = [n for n in names if n != "era"]
    cols = ["functional", "unit"] + [f"gsh_vs_{n}" for n in names] + [f"{n}_vs_era" for n in others]
    with OUT_MODES_CSV.open("w", encoding="utf-8") as f:
        f.write(",".join(cols) + chr(10))
        for row in rows:
            f.write(",".join(row[c] if isinstance(row[c], str) else f"{row[c]:.17g}" for c in cols) + chr(10))

    def fmt(x: float) -> str:
        return "nan" if not np.isfinite(x) else f"{x:.2e}"

    def cell(x: float) -> str:
        return f"{fmt(x)} ({'yes' if np.isfinite(x) and x <= TARGET else 'no'})"

    head = ["functional"] + [f"gsh vs GrafLab {n}" for n in names] + [f"GrafLab {n} vs era" for n in others]
    lines = [
        "# Meridian benchmark: all GrafLab Legendre modes",
        "",
        "EGM2008, degree and order 2190, GRS80, longitude 0 deg, 181 points, h = 0 m. gsh is the same single run "
        "(`compute_functionals`) in every column. GrafLab 2.1.4 is run with `fnALFs` = 3 (era), 2 (msfc) and 1 (sfc).",
        "",
        f"Metric: max |d| / max |reference|, criterion <= {TARGET:g}; the verdict is in brackets. "
        "The last columns compare two GrafLab runs with each other (spread of GrafLab itself).",
        "",
        "| " + " | ".join(head) + " |",
        "|" + "---|" * len(head),
    ]
    for row in rows:
        cells = [row["functional"]] + [cell(row[f"gsh_vs_{n}"]) for n in names] +                 [cell(row[f"{n}_vs_era"]) for n in others]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "dg excludes the two exact poles for the gsh columns (gsh rejects `1/cos(phi)` there); "
              "in the GrafLab-against-GrafLab columns all 181 points are used.", ""]
    OUT_MODES_MD.write_text(chr(10).join(lines), encoding="utf-8")
    print(f"wrote {OUT_MODES_CSV}")
    print(f"wrote {OUT_MODES_MD}")


def main() -> int:
    if not GSH_CSV.exists():
        print(f"missing {GSH_CSV}: run gsh_meridian first", file=sys.stderr)
        return 2
    gsh = read_table(GSH_CSV)
    lat = gsh["lat_deg"]

    if not GL_CSV.exists():
        OUT_MD.write_text(
            "# Meridian benchmark: gsh against GrafLab\n\n"
            "Status: **GrafLab reference not available**. `results/graflab_era_meridian.csv` does not exist; "
            "run `matlab/run_graflab_meridian.m` with GrafLab.m available (see README.md). No numbers are reported.\n",
            encoding="utf-8")
        print(f"{GL_CSV.name} not found: wrote the status note only ({OUT_MD})")
        return 1

    gl = read_table(GL_CSV)

    # (label, gsh column, GrafLab column, conversion to SI, unit label, GrafLab column meaning)
    pairs = [
        ("T", "T_m2s2", "T_m2s2", 1.0, "m^2 s^-2", "functional 5"),
        ("zeta_ell", "zeta_ell_m", "zeta_ell_m", 1.0, "m", "functional 22 (height anomaly ell)"),
        ("Dg_sa (gravity anomaly, sa)", "Dg_sa_ms2", "Dg_sa_mgal", MGAL_TO_SI, "m s^-2", "functional 19"),
        ("dg (gravity disturbance)", "dg_ms2", "dg_mgal", MGAL_TO_SI, "m s^-2", "functional 20"),
        ("V_xx", "V_xx_s2", "Vxx_E", EOTVOS_TO_SI, "s^-2", "functional 14, column 1"),
        ("V_yy", "V_yy_s2", "Vyy_as_exported_E", EOTVOS_TO_SI, "s^-2",
         "functional 14, column 2 (exported as -Vyy, see README)"),
        ("V_zz", "V_zz_s2", "Vzz_E", EOTVOS_TO_SI, "s^-2", "functional 14, column 3"),
        ("V_xz", "V_xz_s2", "Vxz_E", EOTVOS_TO_SI, "s^-2", "functional 15, column 2"),
    ]

    rows = []
    for label, gcol, rcol, conv, unit, meaning in pairs:
        g = gsh[gcol]
        r = gl[rcol] * conv
        m = metrics(g, r)
        m.update({"functional": label, "unit": unit, "graflab_column": meaning,
                  "lat_of_max_abs_deg": lat_of_max(lat, g, r)})
        rows.append(m)

    # Diagnostic 1: sign of the GrafLab Vyy column. The internal GrafLab code 14 exports -Vyy; check which sign
    # satisfies the Laplace equation Vxx + Vyy + Vzz = 0 with the GrafLab columns (a consistency test, not a fit).
    vxx_g = gl["Vxx_E"] * EOTVOS_TO_SI
    vyy_g = gl["Vyy_as_exported_E"] * EOTVOS_TO_SI
    vzz_g = gl["Vzz_E"] * EOTVOS_TO_SI
    trace_gl_exported = vxx_g + vyy_g + vzz_g
    trace_gl_flipped = vxx_g - vyy_g + vzz_g
    trace_gsh = gsh["V_xx_s2"] + gsh["V_yy_s2"] + gsh["V_zz_s2"]
    scale_v = float(np.max(np.abs(gsh["V_xx_s2"])))
    diag = {
        "gsh": float(np.max(np.abs(trace_gsh))) / scale_v,
        "graflab_as_exported": float(np.max(np.abs(trace_gl_exported))) / scale_v,
        "graflab_sign_flipped": float(np.max(np.abs(trace_gl_flipped))) / scale_v,
    }

    # Diagnostic 2: the same gsh V_yy against the GrafLab column with the sign flipped.
    flipped = metrics(gsh["V_yy_s2"], -vyy_g)

    fields = ["functional", "unit", "n", "n_nonfinite", "max_abs", "scale", "max_abs_over_scale", "max_rel_pt",
              "rms", "lat_of_max_abs_deg", "graflab_column"]
    with OUT_CSV.open("w", encoding="utf-8") as f:
        f.write(",".join(fields) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) if not isinstance(row[k], float) else f"{row[k]:.17g}" for k in fields) + "\n")

    def fmt(x: float) -> str:
        return "nan" if not np.isfinite(x) else f"{x:.3e}"

    lines = [
        "# Meridian benchmark: gsh against GrafLab 2.1.4 (era)",
        "",
        "EGM2008 (`data/cache/modelos/EGM2008/EGM2008.gfc`), degree and order 2190, GRS80, one meridian at "
        "longitude 0 deg, 181 points of geodetic latitude -90..+90 deg, h = 0 m.",
        "",
        f"Success criterion: max |d| / max |GrafLab| <= {TARGET:g} (rounding level, stated by the user). "
        "Each row is reported as it is; no tolerance was changed.",
        "",
        "| functional | unit | n | max abs(d) | max abs(GrafLab) | max abs(d)/max abs(GrafLab) | max rel(pt) | rms(d) | "
        "lat of max abs(d) [deg] | meets 1e-12 |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        ok = np.isfinite(row["max_abs_over_scale"]) and row["max_abs_over_scale"] <= TARGET
        lines.append(
            f"| {row['functional']} | {row['unit']} | {row['n']} | {fmt(row['max_abs'])} | {fmt(row['scale'])} | "
            f"{fmt(row['max_abs_over_scale'])} | {fmt(row['max_rel_pt'])} | {fmt(row['rms'])} | "
            f"{row['lat_of_max_abs_deg']:.1f} | {'yes' if ok else 'no'} |")
    lines += [
        "",
        "## Conventions and mapping",
        "",
        "- gsh: `compute_functionals` (`Functional::T, zeta_ell, Dg_sa, dg, V_xx, V_yy, V_zz, V_xz`), SI units.",
        "- GrafLab: `Functional = [5; 22; 19; 20]` and `[14; 15]`, `fnALFs = 3` (extended-range arithmetic). "
        "Units as GrafLab returns them: mGal and Eotvos, converted here to SI.",
        "- dg at the exact poles: gsh does not evaluate it (1/cos(phi)); those two rows are excluded (n = 179).",
        "- zeta: GrafLab functional 23 needs a DTM file, which this benchmark does not have. The comparison uses "
        "zeta_Ell (GrafLab 22) against gsh `zeta_ell` = T/gamma. The full zeta is not compared.",
        "- Dg (gravity anomaly): GrafLab has only the spherical-approximation anomaly (functional 19); the gsh "
        "counterpart is `Dg_sa`. The exact gravity anomaly is not a GrafLab output.",
        "- V_yy: GrafLab code 14 writes `[Vxx, -Vyy, Vzz]` (GrafLab.m, `Pg=[Pg -Vyy(:)]`). The comparison is made "
        "with the column as exported; the sign-flipped comparison is reported below as a diagnostic.",
        "- max rel(pt) is the largest |d|/|GrafLab| at points with a non-zero GrafLab value. It is not meaningful "
        "where the functional passes through zero (e.g. near the equator), so the normalised column is the criterion.",
        "",
        "## Diagnostics",
        "",
        f"- Laplace residual max|Vxx+Vyy+Vzz| / max|Vxx|: gsh {fmt(diag['gsh'])}, GrafLab as exported "
        f"{fmt(diag['graflab_as_exported'])}, GrafLab with the Vyy sign flipped {fmt(diag['graflab_sign_flipped'])}.",
        f"- gsh V_yy against -GrafLab Vyy (sign flipped): max abs(d) {fmt(flipped['max_abs'])}, "
        f"normalised {fmt(flipped['max_abs_over_scale'])}.",
        "",
        "Raw numbers: `comparison.csv`. Inputs: `gsh_meridian.csv`, `graflab_era_meridian.csv`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_MD}")
    modes_report(gsh, pairs)
    for row in rows:
        print(f"{row['functional']:<30} max_abs={fmt(row['max_abs'])}  scale={fmt(row['scale'])}  "
              f"norm={fmt(row['max_abs_over_scale'])}  n={row['n']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
