"""Commission and omission errors, and the error report (errors.hpp, error_report.hpp)."""

import json
import math

import numpy as np
import pytest


def test_commission_error_is_positive_and_finite(gsh, egm2008_sigma_nmax60):
    sigma = egm2008_sigma_nmax60
    pts = gsh.Points.from_geodetic(np.array([0.0, 45.0, 80.0]), np.array([0.0, 0.0, 0.0]), np.zeros(3))
    opts = gsh.ErrorOptions()
    opts.nmin = 2
    opts.nmax = 60
    ce = gsh.commission_error(sigma, pts, gsh.ErrorFunctional.ZetaEll, opts)
    assert ce.shape == (3,)
    assert np.all(np.isfinite(ce)) and np.all(ce > 0)


def test_commission_nmax_above_sigma_raises(gsh, egm2008_sigma_nmax60):
    pts = gsh.Points.from_geodetic(np.array([10.0]), np.array([10.0]), np.array([0.0]))
    opts = gsh.ErrorOptions()
    opts.nmax = 61
    with pytest.raises(ValueError):
        gsh.commission_error(egm2008_sigma_nmax60, pts, gsh.ErrorFunctional.T, opts)


def test_kaula_omission_formula(gsh, egm2008_sigma_nmax60):
    """Kaula rule: R 1e-5 / (nmax + 1) with R the model radius."""
    R = egm2008_sigma_nmax60.R
    assert gsh.omission_kaula_undulation(R, 60) == pytest.approx(R * 1e-5 / 61, rel=1e-12)


def test_degree_variances_match_coefficients(gsh, egm2008_nmax60):
    """sum over m of (C^2 + S^2) at degree n equals the degree variance."""
    model = egm2008_nmax60
    dv = gsh.degree_variances(model)
    n = 30
    total = sum(model.c(n, m) ** 2 + model.s(n, m) ** 2 for m in range(n + 1))
    assert dv[n] == pytest.approx(total, rel=1e-12)


def test_error_report_fields(gsh, egm2008_sigma_nmax60, egm2008_nmax60):
    # Run degree 30 against a 60-degree coefficient model, so the Tscherning-Rapp tail above 30 is included.
    report = gsh.make_error_report(
        egm2008_sigma_nmax60,
        30,
        6371000.0,
        9.81,
        [],
        None,
        gsh.ReportOptions(),
        egm2008_nmax60,
    )
    assert report.nmax == 30
    assert report.model_name == "EGM2008"
    assert report.numerical.available is False
    assert report.omission.kaula == pytest.approx(report.R * 1e-5 / 31, rel=1e-12)
    assert report.omission.has_tr_model is True
    assert report.omission.coefficient_nmax == 60
    assert report.commission.global_N > 0
    parsed = json.loads(gsh.error_report_json(report))
    assert parsed["run"]["nmax"] == 30
    assert parsed["omission"]["tr_model_available"] is True
    assert {"numerical", "omission", "commission"} <= set(parsed)
    assert isinstance(gsh.format_error_report(report), str)


def test_relative_error_and_digit_bound(gsh):
    assert gsh.relative_error(1.0 + 1e-9, 1.0) == pytest.approx(1e-9, rel=1e-6)
    assert gsh.relative_error(0.5, 0.0) == pytest.approx(0.5)
    b = gsh.digit_loss_bound(100, math.pi / 2)
    assert math.isfinite(b) and b > 0
