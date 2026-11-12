"""Error of each engine against the binary128 reference (benchmark/meridian/ref_quad).

Usage (repository root):
    .venv\\Scripts\\python.exe benchmark\\meridian\\compare_quad.py

Inputs in results/: gsh_meridian.csv, graflab_{era,msfc,sfc}_meridian.csv, ref_quad_meridian.csv (exact decimal
latitude) and, if present, ref_quad_meridian_engine_rad.csv (latitude as the double radian value the engines use).
Outputs: comparison_quad.md and comparison_quad.csv.

Metric (same as compare_meridian.py, fixed beforehand): max |engine - reference| / max |reference| over the points
where both are finite. dg is also given in ULP of g (1 ULP of 9.8 m/s^2 = 1.78e-15 m/s^2), because dg is the
difference of two numbers of that size.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = HERE / "results"
MGAL, EOT = 1e-5, 1e-9
ULP_G = 2.0 ** -49  # spacing of doubles in [8, 16): 1.776e-15

# (label, reference column, gsh column, GrafLab column, GrafLab conversion)
ROWS = [
    ("dg", "dg_ms2", "dg_ms2", "dg_mgal", MGAL),
    ("V_xz", "V_xz_s2", "V_xz_s2", "Vxz_E", EOT),
    ("V_xx (control)", "V_xx_s2", "V_xx_s2", "Vxx_E", EOT),
    ("V_zz (control)", "V_zz_s2", "V_zz_s2", "Vzz_E", EOT),
]
ENGINES = ["gsh", "GrafLab era", "GrafLab msfc", "GrafLab sfc"]


def read(path: Path) -> dict[str, np.ndarray]:
    with path.open() as f:
        head = f.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1, dtype=float)
    return {k: data[:, i] for i, k in enumerate(head)}


def stat(eng: np.ndarray, ref: np.ndarray, lat: np.ndarray) -> tuple[float, float, float, float]:
    ok = np.isfinite(eng) & np.isfinite(ref)
    d = np.abs(eng - ref)[ok]
    scale = float(np.max(np.abs(ref[ok])))
    k = int(np.argmax(d))
    return float(d.max()), float(d.max() / scale), float(lat[ok][k]), float(np.sqrt(np.mean(d * d)))


def table(ref: dict[str, np.ndarray], engines: dict[str, dict[str, np.ndarray]], title: str) -> tuple[list[str], list[dict]]:
    lat = ref["lat_deg"]
    out, rows = [f"### {title}", "", "| functional | " + " | ".join(ENGINES) + " |", "|---|" + "---|" * len(ENGINES)], []
    for label, rcol, gcol, lcol, conv in ROWS:
        cells = []
        for name in ENGINES:
            e = engines[name]
            eng = e[gcol] if name == "gsh" else e[lcol] * conv
            mx, norm, lat_max, rms = stat(eng, ref[rcol], lat)
            extra = f", {mx / ULP_G:.1f} ULP(g)" if label == "dg" else ""
            cells.append(f"{norm:.2e} ({'yes' if norm <= 1e-12 else 'no'}) @{lat_max:+.0f}{extra}")
            rows.append({"variant": title, "functional": label, "engine": name, "max_abs": mx, "norm": norm,
                         "lat_of_max": lat_max, "rms": rms})
        out.append(f"| {label} | " + " | ".join(cells) + " |")
    return out + [""], rows


def main() -> int:
    need = ["gsh_meridian.csv", "graflab_era_meridian.csv", "graflab_msfc_meridian.csv", "graflab_sfc_meridian.csv",
            "ref_quad_meridian.csv"]
    miss = [n for n in need if not (R / n).exists()]
    if miss:
        print("missing:", ", ".join(miss), file=sys.stderr)
        return 1
    engines = {"gsh": read(R / "gsh_meridian.csv"), "GrafLab era": read(R / "graflab_era_meridian.csv"),
               "GrafLab msfc": read(R / "graflab_msfc_meridian.csv"), "GrafLab sfc": read(R / "graflab_sfc_meridian.csv")}
    lines = [
        "# Meridian benchmark: engines against the binary128 reference",
        "",
        "Reference: `ref_quad/ref_quad.cpp` (IEEE binary128, ~1e-34), EGM2008 n = 2190, GRS80, lon 0, h = 0, 181 points. "
        "Coefficients, GM, R and GRS80 constants are rounded to double first (as the engines read them).",
        "",
        "Cell = max|engine - reference| / max|reference|, verdict against 1e-12, latitude of the maximum, and for dg "
        "the maximum in ULP of g (1 ULP = 1.78e-15 m/s^2). dg excludes the two poles (gsh rejects 1/cos(phi) there); "
        "V_xx excludes the poles (no d2P closed form there in the reference).",
        "",
    ]
    allrows: list[dict] = []
    for fname, title in [("ref_quad_meridian.csv", "Position: exact decimal latitude"),
                         ("ref_quad_meridian_engine_rad.csv", "Position: latitude as the double radian used by the engines")]:
        p = R / fname
        if not p.exists():
            continue
        t, rows = table(read(p), engines, title)
        lines += t
        allrows += rows
    (R / "comparison_quad.md").write_text("\n".join(lines), encoding="utf-8")
    cols = ["variant", "functional", "engine", "max_abs", "norm", "lat_of_max", "rms"]
    with (R / "comparison_quad.csv").open("w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for r in allrows:
            f.write(",".join(str(r[c]) if isinstance(r[c], str) else f"{r[c]:.17g}" for c in cols) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
