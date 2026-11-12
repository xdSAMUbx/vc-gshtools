"""Create data/points_latsweep.csv: a meridian of points from -85 to +85 degrees latitude.

Purpose: the Legendre recurrences fail at a latitude that depends on the degree (SFC underflows
when cos(latitude) is small, see Holmes and Featherstone 2002; measured in benchmark/alf). The Colorado points lie at
35-40 degrees, where nothing breaks, so this stress set samples every 0.5 degrees from pole to pole
(excluding the last 5 degrees) at a fixed longitude with ellipsoidal height 0.

Run from the repository root:  python benchmark/shs/code/python/make_latitude_sweep.py
"""
import pathlib

out = pathlib.Path(__file__).resolve().parents[2] / "data" / "points_latsweep.csv"
lats = [-85.0 + 0.5 * i for i in range(341)]  # -85.0 ... +85.0
with open(out, "w", encoding="utf-8", newline="\n") as f:
    f.write("id,lat_deg,lon_deg,h_m\n")
    for i, lat in enumerate(lats, start=1):
        f.write(f"{i},{lat:.6f},10.000000,0.0\n")
print(f"wrote {out} ({len(lats)} points)")
