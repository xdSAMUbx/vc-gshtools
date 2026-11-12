"""compute_functionals: agreement with the raw-sum path and error propagation.

The model is EGM2008 truncated at degree 60 (fixture egm2008_nmax60). Points are geodetic (degrees, metres) and
stored by Points as geocentric radians.
"""

import numpy as np
import pytest


def make_points(gsh):
    lat = np.array([45.0, -12.5, 89.8, 0.0])
    lon = np.array([-70.0, 77.0, 0.0, 150.0])
    h = np.array([0.0, 500.0, 0.0, 1200.0])
    return gsh.Points.from_geodetic(lat, lon, h)


def test_V_equals_raw_sum_path(gsh, egm2008_nmax60):
    """(b) V = GM/r S0 from compute_functionals equals the same expression built from synth_points directly."""
    model = egm2008_nmax60
    pts = make_points(gsh)
    prep = gsh.prepare(model, 0, 60)
    raw = gsh.synth_points(prep, pts, 0)
    s0 = raw.f[0]
    v_raw = model.GM / pts.r * s0

    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.V]
    req.nmin = 0
    req.nmax = 60
    res = gsh.compute_functionals(model, pts, req)
    np.testing.assert_allclose(res.values[gsh.Functional.V], v_raw, rtol=1e-12, atol=0.0)


def test_laplace_identity_of_the_tensor(gsh, egm2008_nmax60):
    """(b) V_xx + V_yy + V_zz = 0 away from the sources (harmonic potential)."""
    model = egm2008_nmax60
    pts = make_points(gsh)
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.V_xx, gsh.Functional.V_yy, gsh.Functional.V_zz]
    req.nmax = 60
    v = gsh.compute_functionals(model, pts, req).values
    trace = v[gsh.Functional.V_xx] + v[gsh.Functional.V_yy] + v[gsh.Functional.V_zz]
    scale = max(np.max(np.abs(v[k])) for k in v)
    assert np.max(np.abs(trace)) < 1e-12 * scale


def test_result_structure_and_polar_points(gsh, egm2008_nmax60):
    model = egm2008_nmax60
    pts = make_points(gsh)
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.g, gsh.Functional.g_x, gsh.Functional.g_y, gsh.Functional.g_z, gsh.Functional.xi]
    req.nmax = 60
    res = gsh.compute_functionals(model, pts, req)
    assert set(res.values.keys()) == set(req.outputs)
    for arr in res.values.values():
        assert arr.shape == (len(pts),)
        assert arr.dtype == np.float64
    assert res.polar_points == 1  # lat 89.8 deg is below the cos(phi) threshold
    assert res.points_without_height == 0
    assert "model_synthesis" in res.timing_s
    # g is the magnitude of the gravity vector of the LNOF frame (definition of the functional).
    gx, gy, gz = (res.values[gsh.Functional.g_x], res.values[gsh.Functional.g_y], res.values[gsh.Functional.g_z])
    np.testing.assert_allclose(res.values[gsh.Functional.g], np.sqrt(gx**2 + gy**2 + gz**2), rtol=1e-12)


def test_threads_do_not_change_values(gsh, egm2008_nmax60):
    model = egm2008_nmax60
    pts = make_points(gsh)
    out = []
    for threads in (1, 2):
        req = gsh.FunctionalRequest()
        req.outputs = [gsh.Functional.V, gsh.Functional.dV_r]
        req.nmax = 60
        req.synth.threads = threads
        out.append(gsh.compute_functionals(model, pts, req).values)
    for k in out[0]:
        np.testing.assert_allclose(out[1][k], out[0][k], rtol=1e-12, atol=0.0)


def test_nmax_above_model_raises_value_error(gsh, egm2008_nmax60):
    """(c) nmax > model nmax raises ValueError."""
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.V]
    req.nmax = 61
    with pytest.raises(ValueError):
        gsh.compute_functionals(egm2008_nmax60, make_points(gsh), req)


def test_N_without_height_raises_value_error(gsh, egm2008_nmax60):
    """(c) N (or zeta) without a height source raises ValueError."""
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.N]
    req.nmax = 60
    assert req.height is None
    with pytest.raises(ValueError):
        gsh.compute_functionals(egm2008_nmax60, make_points(gsh), req)


def test_eta_at_the_pole_raises_value_error(gsh, egm2008_nmax60):
    """(c) eta contains 1/cos(lat) and is undefined at the exact pole."""
    pole = gsh.Points(np.array([np.pi / 2]), np.array([0.0]), np.array([6356752.314]))  # radians, metres
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.eta]
    req.nmax = 60
    with pytest.raises(ValueError):
        gsh.compute_functionals(egm2008_nmax60, pole, req)


def test_height_grid_marks_points_outside_as_nan(gsh, egm2008_nmax60):
    """N and zeta are NaN where the height grid has no data; the count is reported."""
    grid = gsh.HeightGrid()
    grid.lat_south, grid.lat_north = 40.0, 50.0
    grid.lon_west, grid.lon_east = 285.0, 295.0  # -75 .. -65 deg
    grid.nlat, grid.nlon = 3, 3
    grid.h = np.full(9, 100.0, dtype=np.float32)
    assert grid.dlat() == pytest.approx(5.0)
    assert grid.interpolate(45.0, -70.0) == pytest.approx(100.0, abs=1e-4)
    assert np.isnan(grid.interpolate(0.0, 0.0))

    pts = gsh.Points.from_geodetic(np.array([45.0, 0.0]), np.array([-70.0, 0.0]), np.array([0.0, 0.0]))
    req = gsh.FunctionalRequest()
    req.outputs = [gsh.Functional.N, gsh.Functional.zeta]
    req.nmax = 60
    req.height = gsh.height_from_grid(grid, 0.0)
    res = gsh.compute_functionals(egm2008_nmax60, pts, req)
    assert res.points_without_height == 1
    assert np.isfinite(res.values[gsh.Functional.N][0])
    assert np.isnan(res.values[gsh.Functional.N][1])


def test_functional_names_and_enum(gsh):
    names = gsh.all_functionals()
    assert len(names) == 45
    assert gsh.functional_name(gsh.Functional.V) == "V"
    assert gsh.functional_name(gsh.Functional.zeta_ell) == "zeta_ell"
    assert gsh.functional_name(gsh.Functional.Dg_sa) == "Dg_sa"
