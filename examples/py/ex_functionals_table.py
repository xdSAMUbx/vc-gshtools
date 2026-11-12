"""ex_functionals_table.py: geopotential functionals at a table of points (gshtools), read and written as CSV.

Python equivalent of examples/cpp/ex_functionals_table.cpp. Same arguments, same CSV columns, same screen table.

Shows: reading a delimited table with gshtools.read_table(), building Points from geodetic coordinates
(Points.from_geodetic), requesting several functionals in one compute_functionals() call, and writing the results.

Usage: python examples/py/ex_functionals_table.py <model.gfc> <points.csv> <out.csv> [nmax]
  points.csv: one header line, then rows "id,lat_deg,lon_deg,h_m" with numeric ids and geodetic coordinates on GRS80.
              The demo table is examples/data/points_demo.csv.
  nmax:       highest degree, >= 2 (default: the degree of the model). Use a low value, for example 200, for a quick run.
              EGM2008 (degree 2190) takes a while to load even for five points.

Output: the CSV file, and a short table on the screen. Units are SI:
  T        disturbing potential [m^2 s^-2]
  zeta_ell height anomaly on the ellipsoid, T/gamma [m]
  Dg_sa    gravity anomaly, spherical approximation [m s^-2]
  dg       gravity disturbance, exact (|g| - gamma) [m s^-2]; undefined at the exact poles
  V_xx, V_zz, V_xz gravitational tensor in the local frame (x north, y west, z up) [s^-2]
N and zeta are not requested: they need a topographic height source, which this example does not use. zeta_ell
(the height anomaly from T) is requested instead.
"""

import argparse
import csv
import sys

import gshtools as gsh

CSV_HEADER = ["id", "lat_deg", "lon_deg", "h_m", "T_m2s2", "zeta_ell_m", "Dg_sa_ms2", "dg_ms2",
              "V_xx_s2", "V_zz_s2", "V_xz_s2"]


def main():
    parser = argparse.ArgumentParser(description="Geopotential functionals at a table of points (gshtools).")
    parser.add_argument("model", help="fully normalised ICGEM .gfc model")
    parser.add_argument("points", help="CSV with header id,lat_deg,lon_deg,h_m")
    parser.add_argument("out", help="output CSV")
    parser.add_argument("nmax", nargs="?", type=int, default=None,
                        help="highest degree, >= 2 (default: the degree of the model)")
    args = parser.parse_args()
    if args.nmax is not None and args.nmax < 2:
        parser.error("nmax must be at least 2")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    ell = gsh.Ellipsoid.GRS80()
    nmax_load = -1 if args.nmax is None else args.nmax
    model = gsh.load_gfc(args.model, nmax_load)

    # Columns: 0 = id, 1 = lat_deg, 2 = lon_deg, 3 = h_m. A line with a non-numeric field is not repaired: stop instead.
    tab = gsh.read_table(args.points, sep=",", header_lines=1, cols=[0, 1, 2, 3])
    if len(tab.skipped_lines):
        raise RuntimeError(f"non-numeric line {tab.skipped_lines[0]} in {args.points}")
    ids, lat, lon, h = tab.columns
    count = len(ids)
    if count == 0:
        raise RuntimeError(f"no points in {args.points}")
    if any(abs(v) >= 90.0 for v in lat):
        raise ValueError("dg is undefined at the exact poles; remove the pole rows from the table")

    pts = gsh.Points.from_geodetic(lat, lon, h, ell.a, ell.e2)

    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.T, gsh.Functional.zeta_ell, gsh.Functional.Dg_sa, gsh.Functional.dg,
                   gsh.Functional.V_xx, gsh.Functional.V_zz, gsh.Functional.V_xz]
    req.nmax = nmax_load
    req.ellipsoid = ell
    res = gsh.compute_functionals(model, pts, req)

    T = res.values[gsh.Functional.T]
    zeta = res.values[gsh.Functional.zeta_ell]
    Dg = res.values[gsh.Functional.Dg_sa]
    dg = res.values[gsh.Functional.dg]
    Vxx = res.values[gsh.Functional.V_xx]
    Vzz = res.values[gsh.Functional.V_zz]
    Vxz = res.values[gsh.Functional.V_xz]

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(CSV_HEADER)
        for i in range(count):
            # 12 significant digits, as the C++ example (ostream precision 12).
            writer.writerow([f"{ids[i]:.12g}", f"{lat[i]:.12g}", f"{lon[i]:.12g}", f"{h[i]:.12g}",
                             f"{T[i]:.12g}", f"{zeta[i]:.12g}", f"{Dg[i]:.12g}", f"{dg[i]:.12g}",
                             f"{Vxx[i]:.12g}", f"{Vzz[i]:.12g}", f"{Vxz[i]:.12g}"])

    print(f"{'id':>4} {'T [m2/s2]':>12} {'zeta [m]':>12} {'Dg [mGal]':>10} {'dg [mGal]':>12} {'Vzz [E]':>10} "
          f"{'Vxz [E]':>12}")
    for i in range(count):
        # Screen units: 1 mGal = 1e-5 m/s^2, 1 Eotvos = 1e-9 s^-2.
        print(f"{ids[i]:4.0f} {T[i]:12.4f} {zeta[i]:10.4f} {Dg[i] * 1e5:10.3f} {dg[i] * 1e5:12.3f} "
              f"{Vzz[i] * 1e9:10.3f} {Vxz[i] * 1e9:12.3f}")
    nmax_text = "model degree" if args.nmax is None else str(args.nmax)
    print(f"points={count}, nmax={nmax_text}, polar points handled by the oracle={res.polar_points}, "
          f"written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
