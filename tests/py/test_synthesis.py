"""Raw-sum synthesis: the grid path, the point path and the oracle are three independent routes to the same sums."""

import math

import pytest

import numpy as np


def test_grid_fft_path_matches_point_path(gsh, egm2008_nmax60):
    """synth_grid (FFT over longitude for a full circle) equals synth_points at the same nodes."""
    nmax = 30
    prep = gsh.prepare(egm2008_nmax60, 0, nmax)
    grid = gsh.Grid.gauss_legendre(nmax, 6371000.0)
    fields, timing = gsh.synth_grid(prep, grid, 0)
    assert timing.total >= 0.0
    assert grid.full_circle()
    pts = gsh.grid_to_points(grid)
    by_points = gsh.synth_points(prep, pts, 0)
    s0_grid = np.asarray(fields.f[0])
    s0_pts = np.asarray(by_points.f[0])
    assert s0_grid.shape == s0_pts.shape == (grid.nlat * grid.nlon,)
    scale = np.max(np.abs(s0_pts))
    np.testing.assert_allclose(s0_grid, s0_pts, rtol=0.0, atol=1e-10 * scale)


def test_oracle_matches_point_path(gsh, egm2008_nmax60):
    """reference_point (full triangle, MSFC) agrees with synth_points at points away from the poles."""
    nmax = 30
    model = egm2008_nmax60
    coefs = gsh.legendre.make_coefs(nmax)
    prep = gsh.prepare(model, 0, nmax)
    lat_deg = np.array([-41.0, 3.5, 66.0])
    lon_deg = np.array([-70.0, 120.0, 15.0])
    r = np.array([6371000.0, 6378137.0, 6400000.0])
    pts = gsh.Points.from_spherical(lat_deg, lon_deg, r)
    by_points = gsh.synth_points(prep, pts, 0)
    for i in range(len(pts)):
        ref = gsh.reference_point(model, coefs, math.radians(lat_deg[i]), math.radians(lon_deg[i]), r[i], 0, nmax)
        assert math.isclose(ref.raw[0], float(by_points.f[0][i]), rel_tol=1e-10, abs_tol=1e-30)


def test_regular_grid_is_regional_and_not_full_circle(gsh):
    grid = gsh.Grid.regular(45.0, 30.0, 7, 250.0, 260.0, 11, 6378136.3)
    assert (grid.nlat, grid.nlon) == (7, 11)
    assert not grid.full_circle()
    assert math.degrees(grid.lat[1]) == pytest.approx(42.5)


def test_isa_name_is_one_of_the_three(gsh):
    assert gsh.active_isa_name() in ("avx512", "avx2", "scalar")
