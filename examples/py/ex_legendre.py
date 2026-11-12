"""ex_legendre.py: fully normalised associated Legendre functions (gshtools.legendre) at one latitude.

Python equivalent of examples/cpp/ex_legendre.cpp. Same arguments, same printed quantities.

Shows: the three recurrences sfc(), msfc() and era() and how they agree; the dispatcher compute(); the sum-rule
identity sum_{m=0}^{n} Pbar_nm^2 = 2n+1; and the latitude derivative dphi(), checked against a central finite
difference of sfc(). No model file is needed.

Usage: python examples/py/ex_legendre.py [lat_deg=45] [nmax=200]
  lat_deg: GEOCENTRIC (spherical) latitude in degrees, -90..90.
  nmax:    maximum degree, 2..1800. This is the range where sfc() is valid, so all recurrences can be compared.

Output (all values dimensionless, fully normalised with 4 pi, no Condon-Shortley phase):
  max |a - b| is the largest absolute difference and max rel is the largest relative difference over the entries
  where the reference is not below 1e-100 (msfc and era differ from sfc only by rounding).
  The sum rule should give a deviation near 1e-13; it holds at every latitude.
  dP/dphi is the derivative with respect to latitude in radians; the finite difference uses a step of 1e-6 rad.

The packed triangle has tri_size(nmax) entries; element (n, m) is at idx_tri(n, m, nmax). The C++ helpers are not
exposed by gshtools, so the index is written here with the same layout (legendre.hpp): column m starts at
m (nmax + 1) - m (m - 1) / 2 and holds n = m .. nmax.
"""

import argparse
import math
import sys

import numpy as np

import gshtools as gsh

MAX_SFC_DEGREE = 1800  # sfc() is accurate up to this degree (legendre.hpp)


def packed_index(n, m, nmax):
    """Position of (n, m) in the packed triangle (same layout as idx_tri in core.hpp)."""
    return m * (nmax + 1) - m * (m - 1) // 2 + (n - m)


def max_diff(a, ref):
    """Largest absolute and relative difference of a from ref (relative only where |ref| > 1e-100)."""
    err = np.abs(a - ref)
    abs_max = float(err.max())
    mask = np.abs(ref) > 1e-100
    rel_max = float((err[mask] / np.abs(ref[mask])).max())
    return abs_max, rel_max


def main():
    parser = argparse.ArgumentParser(description="Fully normalised Legendre functions at one latitude (gshtools).")
    parser.add_argument("lat_deg", nargs="?", type=float, default=45.0,
                        help="geocentric latitude in degrees, -90..90 (default 45)")
    parser.add_argument("nmax", nargs="?", type=int, default=200,
                        help=f"maximum degree, 2..{MAX_SFC_DEGREE} (default 200)")
    args = parser.parse_args()
    if not -90.0 <= args.lat_deg <= 90.0:
        parser.error("latitude must be in [-90, 90] degrees")
    if not 2 <= args.nmax <= MAX_SFC_DEGREE:
        parser.error(f"nmax must be in 2..{MAX_SFC_DEGREE}")
    lat_deg, nmax = args.lat_deg, args.nmax

    phi = math.radians(lat_deg)
    t, u = math.sin(phi), math.cos(phi)  # t = sin(phi), u = cos(phi)
    leg = gsh.legendre
    coefs = leg.make_coefs(nmax)
    size = (nmax + 1) * (nmax + 2) // 2

    P_sfc = leg.sfc(t, u, coefs)
    P_msfc = leg.msfc(t, u, coefs)
    P_era = leg.era(t, u, coefs)
    P_auto = leg.compute(t, u, coefs)
    assert P_sfc.shape == (size,)

    d_msfc = max_diff(P_msfc, P_sfc)
    d_era = max_diff(P_era, P_sfc)
    d_auto = max_diff(P_auto, P_sfc)

    # Sum rule: for every n, sum_m Pbar_nm^2 = 2n+1 on the unit sphere (t^2 + u^2 = 1).
    worst_sum = 0.0
    for n in range(nmax + 1):
        idx = [packed_index(n, m, nmax) for m in range(n + 1)]
        s = float(np.sum(P_sfc[idx] ** 2))
        worst_sum = max(worst_sum, abs(s / (2 * n + 1) - 1.0))

    # Latitude derivative (analytic, dphi) against a central finite difference of sfc().
    dP = leg.dphi(nmax, P_sfc)
    h = 1e-6
    P_plus = leg.sfc(math.sin(phi + h), math.cos(phi + h), coefs)
    P_minus = leg.sfc(math.sin(phi - h), math.cos(phi - h), coefs)
    fd = (P_plus - P_minus) / (2.0 * h)
    worst_fd = float(np.max(np.abs(dP - fd)))

    n_show, m_show = nmax, 2
    i_show = packed_index(n_show, m_show, nmax)

    print(f"ex_legendre: lat={lat_deg:.6f} deg (geocentric), nmax={nmax}, fully normalised (4 pi)")
    print(f"  sfc vs msfc:      max |a-b| = {d_msfc[0]:.3e}, max rel = {d_msfc[1]:.3e}")
    print(f"  sfc vs era:       max |a-b| = {d_era[0]:.3e}, max rel = {d_era[1]:.3e}")
    print(f"  compute vs sfc:   max |a-b| = {d_auto[0]:.3e} (compute uses sfc up to degree {MAX_SFC_DEGREE})")
    print(f"  sum rule:         max |sum_m Pbar_nm^2/(2n+1) - 1| = {worst_sum:.3e} over n = 0..{nmax}")
    print(f"  (n,m) = ({n_show},{m_show}):  P = {P_sfc[i_show]:.12e}, dP/dphi = {dP[i_show]:.12e}, "
          f"finite difference = {fd[i_show]:.12e}")
    print(f"  dP/dphi vs finite difference over the whole triangle: max |diff| = {worst_fd:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
