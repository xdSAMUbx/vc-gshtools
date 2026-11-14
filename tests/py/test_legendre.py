"""Legendre functions: agreement of the algorithms and the layout of the packed triangle.

Arguments are t = sin(phi) and u = cos(phi), phi geocentric latitude in radians (see legendre.hpp).
"""

import math

import numpy as np
import pytest


def idx_tri(n, m, nmax):
    """Packed-triangle index (layout of core.hpp: column m holds n = m..nmax)."""
    col_offset = m * (nmax + 1) - m * (m - 1) // 2
    return col_offset + (n - m)


LATS_DEG = [-89.0, -60.0, -23.5, 0.0, 12.25, 45.0, 79.0]


@pytest.mark.parametrize("lat_deg", LATS_DEG)
def test_sfc_and_msfc_agree(gsh, lat_deg):
    """(a) SFC and MSFC give the same values wherever SFC is valid (nmax <= 1800)."""
    nmax = 200
    phi = math.radians(lat_deg)
    t, u = math.sin(phi), math.cos(phi)
    c = gsh.legendre.make_coefs(nmax)
    p_sfc = gsh.legendre.sfc(t, u, c)
    p_msfc = gsh.legendre.msfc(t, u, c)
    assert p_sfc.shape == p_msfc.shape == (idx_tri(nmax, nmax, nmax) + 1,)
    mask = np.abs(p_sfc) > 1e-200  # entries below double range are not compared
    np.testing.assert_allclose(p_msfc[mask], p_sfc[mask], rtol=1e-11, atol=0.0)


@pytest.mark.parametrize("lat_deg", [-45.0, 10.0, 60.0])
def test_era_agrees_with_msfc(gsh, lat_deg):
    nmax = 100
    phi = math.radians(lat_deg)
    t, u = math.sin(phi), math.cos(phi)
    c = gsh.legendre.make_coefs(nmax)
    p_era = gsh.legendre.era(t, u, c)
    p_msfc = gsh.legendre.msfc(t, u, c)
    mask = np.abs(p_era) > 1e-200
    np.testing.assert_allclose(p_msfc[mask], p_era[mask], rtol=1e-10, atol=0.0)


def test_compute_dispatches_to_sfc_below_1800(gsh):
    nmax = 150
    phi = math.radians(33.0)
    t, u = math.sin(phi), math.cos(phi)
    c = gsh.legendre.make_coefs(nmax)
    np.testing.assert_array_equal(gsh.legendre.compute(t, u, c), gsh.legendre.sfc(t, u, c))


def test_first_terms_match_their_closed_forms(gsh):
    """P(0,0) = 1, P(1,0) = sqrt(3) t, P(1,1) = sqrt(3) u (the normalisation of legendre.hpp)."""
    nmax = 10
    phi = math.radians(28.0)
    t, u = math.sin(phi), math.cos(phi)
    p = gsh.legendre.sfc(t, u, gsh.legendre.make_coefs(nmax))
    assert p[idx_tri(0, 0, nmax)] == pytest.approx(1.0, abs=0.0)
    assert p[idx_tri(1, 0, nmax)] == pytest.approx(math.sqrt(3) * t, rel=1e-14)
    assert p[idx_tri(1, 1, nmax)] == pytest.approx(math.sqrt(3) * u, rel=1e-14)


def test_dphi_matches_central_difference(gsh):
    """The non-singular derivative dP/dphi agrees with a central difference of sfc (radians)."""
    nmax = 40
    c = gsh.legendre.make_coefs(nmax)
    phi0 = math.radians(37.0)
    h = 1e-6

    def p_at(phi):
        return gsh.legendre.sfc(math.sin(phi), math.cos(phi), c)

    p = p_at(phi0)
    dp = gsh.legendre.dphi(nmax, p)
    numeric = (p_at(phi0 + h) - p_at(phi0 - h)) / (2 * h)
    np.testing.assert_allclose(dp, numeric, atol=1e-7, rtol=0.0)


def test_pole_is_finite_and_zero_for_m_ge_1(gsh):
    """At a pole u = 0: all m >= 1 entries are exactly zero, no NaN."""
    nmax = 30
    p = gsh.legendre.sfc(1.0, 0.0, gsh.legendre.make_coefs(nmax))
    assert np.all(np.isfinite(p))
    for n in range(nmax + 1):
        for m in range(1, n + 1):
            assert p[idx_tri(n, m, nmax)] == 0.0


def test_coefs_and_dphi_errors(gsh):
    c = gsh.legendre.make_coefs(20)
    assert c.nmax == 20
    with pytest.raises(ValueError):
        gsh.legendre.dphi(20, np.zeros(5))  # shorter than tri_size(20)


def test_legendre_points_shapes_and_default_selection(gsh):
    """orders=None selects every (n, m): count = tri_size(nmax); dP is None without derivative."""
    nmax = 40
    lats = np.radians([-30.0, 0.0, 75.0])
    r = gsh.legendre_points(lats, nmax)
    count = (nmax + 1) * (nmax + 2) // 2
    assert r.P.shape == (3, count)
    assert r.P.dtype == np.float64
    assert r.dP is None
    assert r.n.shape == r.m.shape == (count,)
    assert r.n.dtype == r.m.dtype == np.int32
    assert int(r.n[0]) == 0 and int(r.m[0]) == 0
    assert all(int(b) >= int(a) for a, b in zip(r.m[:-1], r.m[1:]))  # ascending m


def test_legendre_points_matches_compute(gsh):
    """Each stored value equals the entry of gsh.compute for the same point (exact)."""
    nmax = 60
    lats = [0.0, math.pi / 2, -math.pi / 2, 0.4]
    c = gsh.legendre.make_coefs(nmax)
    r = gsh.legendre_points(lats, nmax, coefs=c)
    for i, phi in enumerate(lats):
        full = gsh.legendre.compute(math.sin(phi), math.cos(phi), c)
        expected = [full[idx_tri(int(n), int(m), nmax)] for n, m in zip(r.n, r.m)]
        np.testing.assert_array_equal(r.P[i], np.array(expected))


def test_legendre_points_subset_and_nmin(gsh):
    """orders=[0, 5, 20] with nmin=7: order m stores degrees max(m, nmin)..nmax."""
    nmax = 30
    r = gsh.legendre_points(np.radians([10.0, -55.0]), nmax, orders=[20, 0, 5], nmin=7)
    expected_count = (nmax - 7 + 1) + (nmax - 7 + 1) + (nmax - 20 + 1)
    assert r.P.shape == (2, expected_count)
    assert sorted(set(int(m) for m in r.m)) == [0, 5, 20]
    assert int(r.n[r.m == 20].min()) == 20
    assert int(r.n[r.m == 5].min()) == 7
    assert int(r.n.min()) >= 7


def test_legendre_points_derivative_equals_dphi(gsh):
    nmax = 40
    phi = math.radians(37.0)
    c = gsh.legendre.make_coefs(nmax)
    r = gsh.legendre_points([phi], nmax, derivative=True, coefs=c)
    full = gsh.legendre.sfc(math.sin(phi), math.cos(phi), c)
    dp = gsh.legendre.dphi(nmax, full)
    assert r.dP.shape == r.P.shape
    expected = [dp[idx_tri(int(n), int(m), nmax)] for n, m in zip(r.n, r.m)]
    np.testing.assert_array_equal(r.dP[0], np.array(expected))


def test_legendre_points_bad_order_raises(gsh):
    with pytest.raises(ValueError):
        gsh.legendre_points([0.1], 10, orders=[11])
    with pytest.raises(ValueError):
        gsh.legendre_points([0.1], 10, orders=[3, 3])
    with pytest.raises(ValueError):
        gsh.legendre_points([0.1], 10, orders=[-1])
