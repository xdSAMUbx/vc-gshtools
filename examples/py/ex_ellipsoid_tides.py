"""ex_ellipsoid_tides.py: GRS80 normal gravity, coordinate conversions and permanent-tide conversions (gshtools).

Python equivalent of examples/cpp/ex_ellipsoid_tides.cpp. Same arguments, same printed quantities.

Shows: Somigliana's closed form against the spherical-harmonic normal gravity of the ellipsoid (Ellipsoid);
geodetic, ECEF and geocentric spherical coordinates (gshtools.coords and Points.from_geodetic); and the conversion of
C20, of the height N, of the ellipsoidal height and of gravity between the tide-free, zero-tide and mean-tide systems
(gshtools.tides). A shift is defined as "value in the target system minus value in the source system".

Usage: python examples/py/ex_ellipsoid_tides.py <model.gfc>
  The model file is read only for its C20 (degree 2). Its header must name the tide system (EGM2008: tide_free).

Output: normal-gravity differences in mGal (1 mGal = 1e-5 m/s^2), round-trip errors of the coordinate conversions, and
the tide shifts. C20 cannot be converted to mean-tide: the library refuses that request, and the example shows the
refusal.
"""

import argparse
import math
import sys

import gshtools as gsh
from gshtools import tides as tide

TF = tide.TideSystem.TideFree
ZT = tide.TideSystem.ZeroTide
MT = tide.TideSystem.MeanTide


def normal_gravity_section(ell):
    print("Normal gravity on GRS80, at the surface (h = 0):")
    for lat in (0.0, 45.0, 80.0):
        lat_c = gsh.coords.geodetic_to_geocentric_latitude(lat, ell)
        r = ell.radius(math.radians(lat_c))
        g_sh = ell.normal_gravity(math.radians(lat_c), r).gamma
        g_closed = ell.somigliana(math.radians(lat))
        print(f"  lat {lat:5.1f} deg: Somigliana {g_closed:.10f} m/s2, spherical harmonics {g_sh:.10f} m/s2, "
              f"difference {(g_sh - g_closed) * 1e5:.3e} mGal")


def coordinates_section(ell):
    lat, lon, h = 45.5, -70.25, 123.4  # geodetic, degrees and metres
    p = gsh.coords.geodetic_to_ecef(lat, lon, h, ell)
    back = gsh.coords.ecef_to_geodetic(p, ell)
    s = gsh.coords.geodetic_to_spherical(lat, lon, h, ell)
    # Points.from_geodetic stores the geocentric latitude (radians) and the radius of the same conversion.
    pts = gsh.Points.from_geodetic([lat], [lon], [h], ell.a, ell.e2)
    print(f"Coordinates of geodetic ({lat:.2f}, {lon:.2f}, {h:.1f} m):")
    print(f"  ECEF: X = {p.x:.4f}, Y = {p.y:.4f}, Z = {p.z:.4f} m")
    print(f"  ECEF back to geodetic: lat {back.lat_deg:.10f}, lon {back.lon_deg:.10f}, h {back.h:.6f} m "
          f"(error in h {back.h - h:.2e} m)")
    print(f"  geocentric spherical: lat {s.lat_deg:.10f} deg, lon {s.lon_deg:.6f} deg, r {s.r:.4f} m")
    print(f"  Points.from_geodetic version (radians): difference in latitude "
          f"{pts.lat[0] - math.radians(s.lat_deg):.2e} rad, in r {pts.r[0] - s.r:.2e} m")


def tides_section(model):
    lat, N, h = 45.0, 30.0, 100.0  # degrees, metres, metres
    zt = tide.convert_model_c20(model, ZT)

    print(f'Permanent-tide systems (file header tide_system = "{model.tide_system}", k20 = {tide.K20_DEFAULT:.2f}):')
    print(f"  C20 tide-free  = {model.c(2, 0):.12e}")
    print(f"  C20 zero-tide  = {zt.c(2, 0):.12e}, shift {zt.c(2, 0) - model.c(2, 0):.6e}")
    print(f"  c20_shift(tide-free -> zero-tide) = {tide.c20_shift(TF, ZT):.6e}")
    try:
        tide.c20_shift(TF, MT)
    except ValueError as e:
        print(f"  c20_shift(tide-free -> mean-tide) refused: {e}")

    print(f"Shifts at lat {lat:.1f} deg (value in target minus value in source):")
    print(f"  height N, tide-free -> zero-tide: {tide.height_shift(lat, TF, ZT):.6f} m, "
          f"tide-free -> mean-tide: {tide.height_shift(lat, TF, MT):.6f} m")
    print(f"  N = {N:.1f} m converted to mean-tide: {tide.convert_geoid_height(N, lat, TF, MT):.6f} m")
    print(f"  ellipsoidal h = {h:.1f} m converted to mean-tide: "
          f"{tide.convert_ellipsoidal_height(h, lat, TF, MT):.6f} m")
    print(f"  gravity: direct tidal term g_T = {tide.gravity_tide_mgal(lat):.6f} mGal, "
          f"tide-free -> mean-tide shift = {tide.gravity_tide_shift_mgal(lat, TF, MT):.6f} mGal")
    print(f'  parse_tide("zero") = {tide.to_string(tide.parse_tide("zero"))}')


def main():
    parser = argparse.ArgumentParser(description="GRS80, coordinates and tide systems (gshtools).")
    parser.add_argument("model", help="ICGEM .gfc model (only its degree-2 terms are read)")
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # the library message contains a combining macron

    ell = gsh.Ellipsoid.GRS80()
    model = gsh.load_gfc(args.model, 2)
    normal_gravity_section(ell)
    coordinates_section(ell)
    tides_section(model)
    return 0


if __name__ == "__main__":
    sys.exit(main())
