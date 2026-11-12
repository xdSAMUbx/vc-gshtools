"""(d) NumPy input handling, error messages, and the coordinate and tide wrappers."""

import numpy as np
import pytest


def test_list_and_array_inputs_agree(gsh):
    lat = [45.0, -12.5]
    lon = [-70.0, 77.0]
    h = [0.0, 500.0]
    from_lists = gsh.Points.from_geodetic(lat, lon, h)
    from_arrays = gsh.Points.from_geodetic(np.array(lat), np.array(lon), np.array(h))
    np.testing.assert_array_equal(from_lists.lat, from_arrays.lat)
    np.testing.assert_array_equal(from_lists.r, from_arrays.r)


def test_float64_contiguous_and_float32_and_strided_agree(gsh):
    """Contiguous float64 is used as it is; float32 and strided views are converted, with the same result."""
    lat = np.linspace(-60, 60, 7)
    lon = np.linspace(-170, 170, 7)
    h = np.zeros(7)
    base = gsh.Points.from_spherical(lat, lon, 6371000.0 + h)
    f32 = gsh.Points.from_spherical(lat.astype(np.float32), lon.astype(np.float32), (6371000.0 + h).astype(np.float32))
    big = np.repeat(lat, 2)  # strided view below picks every second element
    strided = gsh.Points.from_spherical(big[::2], lon, 6371000.0 + h)
    np.testing.assert_allclose(f32.lat, base.lat, rtol=1e-6)
    np.testing.assert_array_equal(strided.lat, base.lat)


def test_non_one_dimensional_input_is_a_clear_value_error(gsh):
    with pytest.raises(ValueError, match="one-dimensional"):
        gsh.Points.from_spherical(np.zeros((2, 2)), np.zeros(4), np.zeros(4))


def test_length_mismatch_is_a_clear_value_error(gsh):
    with pytest.raises(ValueError, match="different sizes"):
        gsh.Points.from_spherical(np.zeros(3), np.zeros(2), np.zeros(3))


def test_points_views_are_read_only(gsh):
    pts = gsh.Points.from_spherical([1.0, 2.0], [3.0, 4.0], [6.4e6, 6.4e6])
    assert pts.lat.flags.writeable is False
    with pytest.raises(ValueError):
        pts.lat[0] = 0.0
    assert pts.size == 2 and len(pts) == 2


def test_model_index_out_of_range_is_index_error(gsh, egm2008_nmax60):
    with pytest.raises(IndexError):
        egm2008_nmax60.c(61, 0)
    with pytest.raises(IndexError):
        egm2008_nmax60.c(2, 3)


def test_read_table_skips_bad_lines(gsh, tmp_path):
    path = tmp_path / "t.csv"
    path.write_text("x,y,z\n1.0,2.0,3.0\n4.0,oops,6.0\n7.0,8.0\n+9.5,10,11\n")
    t = gsh.read_table(str(path), ",", 1, [0, 1, 2])
    assert [c.tolist() for c in t.columns] == [[1.0, 9.5], [2.0, 10.0], [3.0, 11.0]]
    assert list(t.skipped_lines) == [3, 4]


def test_read_table_missing_file_is_runtime_error(gsh, tmp_path):
    with pytest.raises(RuntimeError):
        gsh.read_table(str(tmp_path / "does_not_exist.txt"), None, 0, [0])


def test_coords_round_trip_and_array_form(gsh):
    ell = gsh.Ellipsoid.GRS80()
    p = gsh.coords.geodetic_to_ecef(45.5, -70.25, 123.4, ell)
    back = gsh.coords.ecef_to_geodetic(p, ell)
    assert back.lat_deg == pytest.approx(45.5, abs=1e-10)
    assert back.lon_deg == pytest.approx(-70.25, abs=1e-10)
    assert back.h == pytest.approx(123.4, abs=1e-6)
    x, y, z = gsh.coords.geodetic_to_ecef(np.array([45.5, 0.0]), np.array([-70.25, 0.0]), np.array([123.4, 0.0]), ell)
    assert x[0] == pytest.approx(p.x, rel=1e-15)
    assert z.shape == (2,)


def test_tides_mean_tide_c20_raises(gsh):
    TS = gsh.tides.TideSystem
    with pytest.raises(ValueError):
        gsh.tides.c20_shift(TS.MeanTide, TS.TideFree)
    assert gsh.tides.parse_tide("zero_tide") == TS.ZeroTide


def test_normal_gravity_matches_somigliana(gsh):
    ell = gsh.Ellipsoid.GRS80()
    g_closed = ell.somigliana(np.radians(45.0))
    assert g_closed == pytest.approx(9.80619, abs=1e-4)
