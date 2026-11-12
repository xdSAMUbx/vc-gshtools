"""Shared helpers for the pyshtools and pyharm runners of the SHS benchmark.

This module imports neither gsh_tools nor any spherical-harmonic engine, so it can be loaded
next to either library. Everything the two runners must agree on lives here:

  * GRS80 constants (same values as src/cpp/gsh_utils/ellipsoid.hpp).
  * Geodetic -> geocentric conversion (exact, on the GRS80 ellipsoid).
  * Normal potential U (gravitational part of the GRS80 level ellipsoid): closed form in
    oblate-spheroidal coordinates (u, beta), including the q-function flattening term
    (Heiskanen & Moritz, Physical Geodesy). Only the gravitational part is used; the
    centrifugal part of the normal potential is not part of T = V - U.
  * Normal gravity gamma_SH(r, phi'): the project's spherical-harmonic series of the ellipsoid
    (zonal terms C_{2k,0}, k = 0..10), ported line by line from
    Ellipsoid::normal_gravity in src/cpp/gsh_utils/ellipsoid.cpp.
  * Gravity disturbance dg = |g| - gamma_SH with the same centrifugal terms as
    src/cpp/gsh_expansion/functionals.cpp (function `dg` in the gm >= 1 block).
  * Output writers for the result CSV, its companion .md and results/timings.csv.

Units: lengths m, angles rad inside the functions (geocentric latitude phi' and longitude
lambda), gravity m/s^2 inside, mGal (1e-5 m/s^2) and m^2/s^2 at the output.
"""

import csv
import datetime
import math
import os

import numpy as np

# GRS80 constants. Same values as src/cpp/gsh_utils/ellipsoid.hpp (Ellipsoid defaults).
GRS80_A = 6378137.0                 # semi-major axis [m]
GRS80_E2 = 0.00669438002290         # first eccentricity squared [-]
GRS80_GM = 3.986005e14              # GM of the ellipsoid [m^3/s^2]
GRS80_J2 = 1.08263e-3               # dynamic form factor [-]
GRS80_OMEGA = 7.292115e-5           # angular velocity [rad/s]
MGAL_PER_MS2 = 1.0e5                # 1 m/s^2 = 1e5 mGal

POINT_COLUMNS = ['id', 'lat_deg', 'lon_deg', 'h_m']
RESULT_COLUMNS = ['id', 'lat_deg', 'lon_deg', 'h_m', 'V_m2s2', 'T_m2s2', 'dg_mgal']
FUNCTIONAL_TAG = 'V_T_dg'   # timings notes: all three functionals come from the same run
TIMING_COLUMNS = ['engine', 'recurrence', 'set', 'nmax', 'n_points', 'threads', 'repeat',
                  'wall_s', 'cpu_s', 'notes']


# ----------------------------------------------------------------------------------------
# Input points
# ----------------------------------------------------------------------------------------
def read_points(path):
    """Read a points CSV with header ``id,lat_deg,lon_deg,h_m``.

    lat/lon are geodetic (GRS80), h is the ellipsoidal height in metres.
    Returns a dict of numpy arrays: ids (str), lat_deg, lon_deg, h_m (float64).
    """
    with open(path, 'r', newline='', encoding='utf-8') as fh:
        reader = csv.reader(fh)
        header = [c.strip() for c in next(reader)]
        if header != POINT_COLUMNS:
            raise ValueError(f'{path}: header must be {POINT_COLUMNS}, got {header}')
        rows = [row for row in reader if row and any(c.strip() for c in row)]
    ids = np.array([row[0].strip() for row in rows], dtype=object)
    lat = np.array([float(row[1]) for row in rows], dtype=np.float64)
    lon = np.array([float(row[2]) for row in rows], dtype=np.float64)
    h = np.array([float(row[3]) for row in rows], dtype=np.float64)
    if np.any(np.abs(lat) > 90.0):
        raise ValueError(f'{path}: latitude outside [-90, 90]')
    return {'ids': ids, 'lat_deg': lat, 'lon_deg': lon, 'h_m': h}


def geodetic_to_geocentric(lat_deg, lon_deg, h_m):
    """Exact GRS80 geodetic -> geocentric conversion.

    Returns (r, phi_c, lam) with r the geocentric distance [m], phi_c the geocentric latitude
    [rad] (atan2(Z, p), not the reduced latitude), lam the longitude [rad].
    """
    phi = np.radians(lat_deg)
    lam = np.radians(lon_deg)
    sin_phi = np.sin(phi)
    n = GRS80_A / np.sqrt(1.0 - GRS80_E2 * sin_phi * sin_phi)   # prime vertical radius
    p = (n + h_m) * np.cos(phi)                                  # distance from the Z axis
    z = (n * (1.0 - GRS80_E2) + h_m) * sin_phi
    r = np.hypot(p, z)
    phi_c = np.arctan2(z, p)
    return r, phi_c, lam


# ----------------------------------------------------------------------------------------
# Normal field
# ----------------------------------------------------------------------------------------
def normal_potential_closed(r, phi_c):
    """Gravitational part V_n of the GRS80 normal potential [m^2/s^2], closed form.

    Oblate-spheroidal coordinates (u, beta) of the point, with E = a e and b = a sqrt(1-e^2):
        u^2 = [(r^2 - E^2) + sqrt((r^2 - E^2)^2 + 4 E^2 Z^2)] / 2,   Z = r sin(phi_c),
        sin^2(beta) = Z^2 / u^2.
    Level-ellipsoid normal potential (Heiskanen & Moritz, Physical Geodesy):
        V_n = (GM/E) arctan(E/u) + 1/2 omega^2 a^2 (q(u)/q(b)) (sin^2 beta - 1/3),
        q(x) = 1/2 [ (1 + 3 x^2/E^2) arctan(E/x) - 3 x/E ].
    The second term is the flattening correction of the level ellipsoid. Without it the
    first term alone is the potential of a homogeneous spheroid (J2 = e^2/5), not GRS80.
    Verified against the project's zonal series (spec section 6; ellipsoid.cpp): agreement
    to ~1e-6 m^2/s^2 on a 6e7 m^2/s^2 level for geocentric latitudes 0..90 deg.
    The centrifugal term 1/2 omega^2 (X^2 + Y^2) is NOT included: T = V - V_n.
    """
    e_lin = GRS80_A * math.sqrt(GRS80_E2)
    b_semi = GRS80_A * math.sqrt(1.0 - GRS80_E2)
    z = r * np.sin(phi_c)
    d = r * r - e_lin * e_lin
    u2 = 0.5 * (d + np.sqrt(d * d + 4.0 * e_lin * e_lin * z * z))
    u = np.sqrt(u2)
    sin2_beta = z * z / u2

    def q_func(x):
        return 0.5 * ((1.0 + 3.0 * x * x / (e_lin * e_lin)) * np.arctan(e_lin / x)
                      - 3.0 * x / e_lin)

    q0 = q_func(b_semi)
    flattening = (0.5 * GRS80_OMEGA ** 2 * GRS80_A ** 2 * (q_func(u) / q0)
                  * (sin2_beta - 1.0 / 3.0))
    return GRS80_GM / e_lin * np.arctan(e_lin / u) + flattening


def ellipsoid_zonal(k):
    """Fully normalised zonal coefficient C^Ell_{2k,0}, k = 0..10 (project formula, sec. 6)."""
    if k == 0:
        return 1.0
    c20 = -GRS80_J2 / math.sqrt(5.0)
    if k == 1:
        return c20
    e2 = GRS80_E2
    sign = 1.0 if k % 2 == 0 else -1.0
    return (sign * 3.0 * e2 ** k / ((2 * k + 1) * (2 * k + 3) * math.sqrt(4 * k + 1))
            * (1.0 - k - 5.0 * math.sqrt(5.0) * k * c20 / e2))


def normal_gravity_gamma(r, phi_c):
    """Normal gravity gamma_SH(r, phi') [m/s^2], ported from Ellipsoid::normal_gravity.

    Zonal spherical-harmonic series of the ellipsoid up to degree 20, with the centrifugal
    term of the rotating frame added as in the C++ code. dP_n0/dphi = sqrt(n(n+1)/2) P_n1
    is used to avoid the pole singularity. Legendre functions carry no Condon-Shortley phase.
    """
    r = np.asarray(r, dtype=np.float64)
    phi_c = np.asarray(phi_c, dtype=np.float64)
    t = np.sin(phi_c)
    u = np.cos(phi_c)
    pole = np.abs(phi_c) == np.pi / 2
    t = np.where(pole, np.sign(phi_c), t)
    u = np.where(pole, 0.0, u)
    q = GRS80_A / r

    p0 = [np.ones_like(t), math.sqrt(3.0) * t]
    p1 = [np.zeros_like(t), math.sqrt(3.0) * u]
    for n in range(2, 21):
        nn = float(n)
        p0.append(math.sqrt(4 * nn * nn - 1) / nn * t * p0[n - 1]
                  - (nn - 1) / nn * math.sqrt((2 * nn + 1) / (2 * nn - 3)) * p0[n - 2])
        a1 = math.sqrt((2 * nn + 1) * (2 * nn - 1) / ((nn - 1) * (nn + 1)))
        b1 = (math.sqrt((2 * nn + 1) * nn * (nn - 2) / ((nn - 1) * (nn + 1) * (2 * nn - 3)))
              if n >= 3 else 0.0)
        p1.append(a1 * t * p1[n - 1] - b1 * p1[n - 2])

    sr = np.zeros_like(t)
    sp = np.zeros_like(t)
    qn = np.ones_like(t)
    for k in range(0, 11):
        n = 2 * k
        c = ellipsoid_zonal(k)
        sr += qn * (n + 1) * c * p0[n]
        sp += qn * c * math.sqrt(n * (n + 1) / 2.0) * p1[n]
        qn = qn * q * q
    dvr = -GRS80_GM / (r * r) * sr
    dvp = GRS80_GM / r * sp
    dvcr = GRS80_OMEGA ** 2 * r * u * u
    dvcp = -GRS80_OMEGA ** 2 * r * r * u * t
    g_x = (dvp + dvcp) / r
    g_z = dvr + dvcr
    return np.hypot(g_x, g_z)


def gravity_disturbance_mgal(v_r, v_lat, v_lon, r, phi_c):
    """Gravity disturbance dg = |g| - gamma_SH in mGal, with the project's definition.

    Inputs are the gradient of the gravitational potential V of the model (no centrifugal
    term) in the geocentric spherical frame:
        v_r   = dV/dr,   v_lat = dV/dphi' (geocentric latitude),   v_lon = dV/dlambda.
    The centrifugal gradient of W = V + 1/2 omega^2 r^2 cos^2 phi' is added here, then
    g = |(g_x, g_y, g_z)| with g_x = (V_phi + Vc_phi)/r, g_y = -V_lambda/(r cos phi'),
    g_z = V_r + Vc_r, exactly as functionals.cpp does.
    """
    cos_p = np.cos(phi_c)
    sin_p = np.sin(phi_c)
    w2 = GRS80_OMEGA ** 2
    vc_r = w2 * r * cos_p * cos_p
    vc_lat = -w2 * r * r * cos_p * sin_p
    g_x = (v_lat + vc_lat) / r
    g_y = -v_lon / (r * cos_p)
    g_z = v_r + vc_r
    g = np.sqrt(g_x * g_x + g_y * g_y + g_z * g_z)
    return (g - normal_gravity_gamma(r, phi_c)) * MGAL_PER_MS2


def disturbing_potential(v, r, phi_c):
    """T = V - U [m^2/s^2] with U from :func:`normal_potential_closed`."""
    return v - normal_potential_closed(r, phi_c)


# ----------------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------------
def write_result_csv(path, points, v_vals, t_vals, dg_vals):
    """Write the result ``id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal`` with 15 significant digits.

    V and T in m^2/s^2, dg in mGal. One file per engine / recurrence / point set; all three
    functionals come from the same run.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        fh.write(','.join(RESULT_COLUMNS) + '\n')
        for i in range(len(t_vals)):
            fh.write(f"{points['ids'][i]},{points['lat_deg'][i]:.15g},{points['lon_deg'][i]:.15g},"
                     f"{points['h_m'][i]:.15g},{v_vals[i]:.15g},{t_vals[i]:.15g},{dg_vals[i]:.15g}\n")


def write_result_md(md_path, csv_name, engine_label, command, params, model_info, n_points, notes):
    """Write the companion .md of a result CSV: column table plus the exact command/params.

    The table covers the three columns V_m2s2, T_m2s2 and dg_mgal of the one result CSV.
    """
    os.makedirs(os.path.dirname(os.path.abspath(md_path)), exist_ok=True)
    lines = [
        f'# {csv_name}',
        '',
        f'Engine: **{engine_label}**. Generated by the runner script on '
        f'{datetime.datetime.now().isoformat(timespec="seconds")}. V, T and dg are computed in the same run.',
        '',
        '## Columns',
        '',
        '| column | unit | meaning |',
        '|---|---|---|',
        '| id | - | point identifier from the input CSV |',
        '| lat_deg | deg | geodetic latitude on GRS80 (input, echoed) |',
        '| lon_deg | deg | longitude (input, echoed) |',
        '| h_m | m | ellipsoidal height on GRS80 (input, echoed) |',
        '| V_m2s2 | m^2/s^2 | gravitational potential V of the model at the point (no normal potential '
        'subtracted; degree 0 included, C00 as in the file) |',
        '| T_m2s2 | m^2/s^2 | disturbing potential T = V - U at the point (V: model, '
        'U: GRS80 normal potential, closed form) |',
        '| dg_mgal | mGal | gravity disturbance dg = \\|g\\| - gamma_SH at the point '
        '(g includes the centrifugal term; gamma_SH: project series) |',
        '',
        '## Command and parameters',
        '',
        '```',
        command,
        '```',
        '',
    ]
    for key, value in params:
        lines.append(f'- {key}: `{value}`')
    lines += ['', '## Model', '']
    for key, value in model_info:
        lines.append(f'- {key}: `{value}`')
    lines += ['', f'Points: {n_points}.', '']
    if notes:
        lines += ['## Notes', ''] + [f'- {n}' for n in notes] + ['']
    with open(md_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines))


def append_timing(path, row):
    """Append one repetition to results/timings.csv (header written if the file is new)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    new = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, 'a', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        if new:
            writer.writerow(TIMING_COLUMNS)
        writer.writerow([row[c] for c in TIMING_COLUMNS])
