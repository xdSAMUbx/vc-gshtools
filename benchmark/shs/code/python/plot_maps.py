"""Draw the final spatial comparison on the Colorado grid: GrafLab (sfc, msfc, era) next to gsh fast.

Reads raw/maps/{graflab_<rec>_grid.csv, gsh_fast_grid_t12.csv} (made by run_maps.ps1) and writes
img/map_colorado_grid_graflab_vs_gsh.png: rows V, T and dg; one colour scale per row; a last row shows the
difference gsh fast - GrafLab era (its own scale).
"""
import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QUANT = [("V_m2s2", "V", "m^2/s^2"), ("T_m2s2", "T", "m^2/s^2"), ("dg_mgal", "dg", "mGal")]


def load(path: pathlib.Path):
    rows = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    return rows[~np.isnan(rows["V_m2s2"].astype(float))] if rows.size else rows


def to_grid(rows, column, lats, lons):
    g = np.full((len(lats), len(lons)), np.nan)
    ia = np.rint((rows["lat_deg"] - lats[0]) / 0.2).astype(int)
    ja = np.rint((rows["lon_deg"] - lons[0]) / 0.2).astype(int)
    g[ia, ja] = rows[column]
    return g


def main() -> int:
    here = pathlib.Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=str(here / "results" / "raw" / "maps"))
    ap.add_argument("--img", default=str(here / "img"))
    args = ap.parse_args()
    raw, img = pathlib.Path(args.raw), pathlib.Path(args.img)
    files = [("GrafLab sfc", "graflab_sfc_grid.csv"), ("GrafLab msfc", "graflab_msfc_grid.csv"),
             ("GrafLab era", "graflab_era_grid.csv"), ("gsh fast", "gsh_fast_grid_t12.csv")]
    data = [(label, load(raw / f)) for label, f in files if (raw / f).exists()]
    if len(data) < 2:
        print("need at least two engines in", raw, file=sys.stderr)
        return 1
    lats = np.unique(np.round(data[0][1]["lat_deg"], 6))
    lons = np.unique(np.round(data[0][1]["lon_deg"], 6))
    ref = dict(data).get("GrafLab era")
    fast = dict(data).get("gsh fast")
    nrow = 3 + (1 if ref is not None and fast is not None else 0)
    fig, axes = plt.subplots(nrow, len(data), figsize=(3.7 * len(data), 3.1 * nrow), squeeze=False,
                             layout="constrained")
    ext = [lons[0], lons[-1], lats[0], lats[-1]]
    for r, (col, short, unit) in enumerate(QUANT):
        grids = [to_grid(rows, col, lats, lons) for _, rows in data]
        lo, hi = min(np.nanmin(g) for g in grids), max(np.nanmax(g) for g in grids)
        im = None
        for c, ((label, _), g) in enumerate(zip(data, grids)):
            ax = axes[r][c]
            im = ax.imshow(g, origin="lower", extent=ext, vmin=lo, vmax=hi, cmap="viridis", aspect="auto")
            ax.set_title(f"{label}: {short}", fontsize=9)
            ax.tick_params(labelsize=7)
            if c == 0:
                ax.set_ylabel("latitude [deg]", fontsize=8)
        fig.colorbar(im, ax=axes[r].tolist(), label=f"{short} [{unit}]", shrink=0.85)
    if nrow == 4:
        for c, (col, short, unit) in enumerate(QUANT):
            d = to_grid(fast, col, lats, lons) - to_grid(ref, col, lats, lons)
            m = np.nanmax(np.abs(d)) or 1.0
            ax = axes[3][c]
            im = ax.imshow(d, origin="lower", extent=ext, vmin=-m, vmax=m, cmap="coolwarm", aspect="auto")
            ax.set_title(f"fast - era, {short}: max |d| {m:.1e} {unit}", fontsize=8)
            ax.tick_params(labelsize=7)
            fig.colorbar(im, ax=ax, shrink=0.85)
        for c in range(3, len(data)):
            axes[3][c].axis("off")
    for c in range(len(data)):
        axes[nrow - 1 if nrow == 3 else 2][c].set_xlabel("longitude [deg]", fontsize=8)
    fig.suptitle("Colorado grid 0.2 deg, ellipsoidal height 0, nmax 2190, XGM2019e: GrafLab and gsh fast", fontsize=10)
    out = img / "map_colorado_grid_graflab_vs_gsh.png"
    img.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor="white")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
