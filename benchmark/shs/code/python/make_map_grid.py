"""Create data/points_grid.csv: a regular 0.2 degree grid over Colorado (35-40 N, 110-102 W), ellipsoidal height 0.

Purpose: the Colorado benchmark points are scattered survey tracks, so they cannot be drawn as a continuous map.
This grid (26 x 41 = 1066 points) is synthesised by GrafLab and by gsh only to draw the final spatial comparison
(img/map_colorado_grid_graflab_vs_gsh.png). It is not part of the timed benchmark.

Run from the repository root:  python benchmark/shs/code/python/make_map_grid.py
"""
import pathlib

out = pathlib.Path(__file__).resolve().parents[2] / "data" / "points_grid.csv"
lats = [35.0 + 0.2 * i for i in range(26)]    # 35.0 ... 40.0
lons = [-110.0 + 0.2 * j for j in range(41)]  # -110.0 ... -102.0
with open(out, "w", encoding="utf-8", newline="\n") as f:
    f.write("id,lat_deg,lon_deg,h_m\n")
    k = 0
    for la in lats:
        for lo in lons:
            k += 1
            f.write(f"{k},{la:.6f},{lo:.6f},0.0\n")
print(f"wrote {out} ({k} points)")
