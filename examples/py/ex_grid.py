"""ex_grid.py: synthesis on latitude-longitude grids (gshtools) of the potential V and the disturbing potential T.

Python equivalent of examples/cpp/ex_grid.cpp. Same arguments, same printed quantities, same regional CSV.

Shows: a global Gauss-Legendre grid (Grid.gauss_legendre), a regional regular grid (Grid.regular), the synthesis on
grids with synth_grid() and its GridTiming report, and the disturbing model (disturbing_model).

The values of V and T on the grids are taken from compute_functionals() applied to the grid nodes
(grid_to_points). gshtools does not expose the conversion V = GM/r S0 on grids, so this example does not write it
out; synth_grid() is called for the timing report and the raw sums, as in the C++ example.

Usage: python examples/py/ex_grid.py <model.gfc> <regional.csv> [nmax=90] [--plot out.png]
  nmax: truncation degree of the model and degree of the global grid, >= 2. The Gauss-Legendre grid has nmax+1
        latitudes and 2(nmax+1) longitudes. The regional grid is fixed: 201 x 201 nodes, latitudes 45 to 35 deg N
        and longitudes 110 to 100 deg W (step 0.05 deg). Its values are written to regional.csv, one node per line.
  --plot out.png: optional. Writes a figure of V and T on the regional grid (needs matplotlib, imported only then).

Output: grid sizes, synthesis timing, sample values of V and T at the radius R of the model (m^2 s^-2), and the
regional CSV with columns lat_sph_deg (spherical latitude), lon_deg, V_m2s2, T_m2s2. Timing lumps and fft are
thread-seconds; total is wall-clock time.
"""

import argparse
import sys

import numpy as np

import gshtools as gsh


def synthesise_fields(model_prep, dist_prep, grid):
    """Raw synthesis of the model and of the disturbing model on a grid, with the timing report of both."""
    _, tv = gsh.synth_grid(model_prep, grid, 0)
    _, tt = gsh.synth_grid(dist_prep, grid, 0)
    return tv, tt


def print_timing(name, tv, tt):
    print(f"{name} synthesis (model + disturbing model): wall {tv.total + tt.total:.3f} s, "
          f"lumps {tv.lumps + tt.lumps:.3f} thread-s, fft {tv.fft + tt.fft:.3f} thread-s")


def values_on_grid(model, grid):
    """V and T at every node of the grid (row by row, latitude slow), with the node points."""
    pts = gsh.grid_to_points(grid)
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.V, gsh.Functional.T]
    req.nmax = model.nmax
    res = gsh.compute_functionals(model, pts, req)
    return pts, res.values[gsh.Functional.V], res.values[gsh.Functional.T]


def print_sample(grid, V, T, row, col):
    idx = row * grid.nlon + col
    print(f"  row {row:4d} (lat_sph {np.degrees(grid.lat[row]):9.4f} deg), col {col:4d}: "
          f"V = {V[idx]:.6e} m^2/s^2, T = {T[idx]:.6e} m^2/s^2")


def write_plot(path, nmax, nlat, nlon, pts, V, T):
    """Two panels, V and T on the regional grid. matplotlib is imported here, only when a plot is requested."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lat = np.degrees(pts.lat).reshape(nlat, nlon)
    lon = np.degrees(pts.lon).reshape(nlat, nlon)
    extent = [lon.min(), lon.max(), lat.min(), lat.max()]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for ax, field, title in ((axes[0], V, "V [m^2/s^2]"), (axes[1], T, "T [m^2/s^2]")):
        im = ax.imshow(field.reshape(nlat, nlon), origin="upper", extent=extent, cmap="viridis", aspect="auto")
        ax.set_title(title)
        ax.set_xlabel("longitude [deg]")
        ax.set_ylabel("spherical latitude [deg]")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle(f"Regional grid, model nmax = {nmax}")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Gauss-Legendre and regional grids of V and T (gshtools).")
    parser.add_argument("model", help="fully normalised ICGEM .gfc model")
    parser.add_argument("regional_csv", help="output CSV of the regional grid")
    parser.add_argument("nmax", nargs="?", type=int, default=90, help="truncation degree, >= 2 (default 90)")
    parser.add_argument("--plot", metavar="OUT.PNG", help="also write a figure of the regional grid (needs matplotlib)")
    args = parser.parse_args()
    if args.nmax < 2:
        parser.error("nmax must be at least 2")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    nmax = args.nmax

    model = gsh.load_gfc(args.model, nmax)
    dist = gsh.disturbing_model(model)
    prep_model = gsh.prepare(model)
    prep_dist = gsh.prepare(dist)
    R = model.R

    # Global grid on the sphere of radius R.
    gl = gsh.Grid.gauss_legendre(nmax, R)
    tv, tt = synthesise_fields(prep_model, prep_dist, gl)
    _, V_gl, T_gl = values_on_grid(model, gl)
    print(f"Gauss-Legendre grid, nmax={nmax}: {gl.nlat} latitudes x {gl.nlon} longitudes = {gl.nlat * gl.nlon} nodes")
    print_timing("global", tv, tt)
    print_sample(gl, V_gl, T_gl, 0, 0)
    print_sample(gl, V_gl, T_gl, gl.nlat // 2, 0)

    # Regional regular grid: direct summation path (it is not a full circle).
    reg = gsh.Grid.regular(45.0, 35.0, 201, -110.0, -100.0, 201, R)
    tv, tt = synthesise_fields(prep_model, prep_dist, reg)
    pts_reg, V_reg, T_reg = values_on_grid(model, reg)
    print(f"Regional grid: {reg.nlat} x {reg.nlon} nodes, full circle = {'yes' if reg.full_circle() else 'no'}")
    print_timing("regional", tv, tt)
    print_sample(reg, V_reg, T_reg, 0, 0)
    print_sample(reg, V_reg, T_reg, reg.nlat - 1, reg.nlon - 1)

    # Node positions from the library (geocentric latitude and longitude in radians, converted for output).
    lat_deg = np.degrees(pts_reg.lat)
    lon_deg = np.degrees(pts_reg.lon)
    with open(args.regional_csv, "w", encoding="utf-8") as f:  # text mode: platform line ends, as the C++ example
        f.write("lat_sph_deg,lon_deg,V_m2s2,T_m2s2\n")
        for i in range(len(V_reg)):
            f.write(f"{lat_deg[i]:.12g},{lon_deg[i]:.12g},{V_reg[i]:.12g},{T_reg[i]:.12g}\n")
    print(f"regional values written to {args.regional_csv}")

    if args.plot:
        write_plot(args.plot, nmax, reg.nlat, reg.nlon, pts_reg, V_reg, T_reg)
        print(f"figure written to {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
