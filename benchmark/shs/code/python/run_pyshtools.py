"""pyshtools runner for the SHS benchmark (definitions: benchmark/shs/results/csv.md, engine 'pyshtools').

One run evaluates at each point the three functionals defined in results/csv.md from the tide-free model, in the same
synthesis pass, and writes one result CSV:
  * V_m2s2 = GM/r sum_n (R/r)^n S_n    (gravitational potential of the model, degree 0 included, C00 from the file)
  * T_m2s2 = V - U                     (disturbing potential, m^2/s^2; V is the value above, not a second synthesis)
  * dg_mgal = |g| - gamma_SH           (gravity disturbance, mGal; from the gradient of the same run)

Usage:
  python code/python/run_pyshtools.py --gfc <tide_free.gfc> --points <csv> --nmax N
      --out results/pyshtools_native_<set>.csv --timings results/timings.csv
      --set colorado|latsweep --repeat R --threads K

Library calls used (pyshtools 4.14.1, all public names):
  * pyshtools.SHGravCoeffs.from_file(path, format='icgem', lmax=N, normalization='4pi',
    csphase=1)              -> model reader (pyshtools' own GFC/ICGEM reader)
  * SHGravCoeffs.to_array(normalization='4pi', csphase=1, lmax=N) -> cilm[2, N+1, N+1]
  * pyshtools.expand.MakeGridPoint(cilm, lat, lon, lmax, norm=1, csphase=1)
                            -> sum_{n,m} Cbar cos + Sbar sin Pbar(sin phi') at (phi', lambda)
  * pyshtools.gravmag.MakeGravGridPoint(cilm, gm, r0, r, lat, lon, lmax, omega=0)
                            -> (B_r, B_theta, B_phi) = grad V in geocentric spherical components

Conventions and how each difference was resolved:
  * Normalisation: fully normalised 4pi (geodesy), the GFC convention. Checked against an
    independent numpy recursion on a synthetic model: exact agreement.
  * Condon-Shortley phase: none (csphase=1 in pyshtools means "do NOT apply (-1)^m").
    This is the GFC convention. Checked as above.
  * Radius: the library needs r0 and GM from the model. The model's GM and R are read by
    pyshtools; the script checks them against the header (gfc_reader.py) to 1e-9 relative
    (GM) and 1e-6 m (R).
  * Point radius: pyshtools' MakeGridPoint has NO radius argument, it returns the
    unit-radius sum. The potential at radius r is V = GM/r * sum_n (R/r)^n S_n. The script
    scales the degree-n coefficients by (R/r)^n and groups the points by identical r
    (np.unique). Cost: one scaled copy of cilm per distinct radius. This copy is inside the
    timed block. For points at h = 0 the number of distinct radii equals the number of distinct
    latitudes.
  * Gradient: MakeGravGridPoint is called with omega = 0, so no centrifugal term. Its
    components are converted to partial derivatives of V: dV/dr = B_r,
    dV/dphi' = -r B_theta (theta = colatitude), dV/dlambda = r cos(phi') B_phi. The
    centrifugal term and the normal gravity are added in shs_common.gravity_disturbance_mgal.
  * Geometry: geodetic (GRS80) -> geocentric (phi', r) is done exactly in shs_common.
    pyshtools receives the geocentric latitude in degrees (its functions take degrees).
  * Ellipsoidal handling: V, T and dg are evaluated at the actual point (r includes h), not on
    the reference ellipsoid. This matches the project (functionals.cpp, dg uses the point's
    R(i); N and zeta use the foot point, not used here).
  * U: closed form of the GRS80 level-ellipsoid gravitational potential (shs_common).
  * gamma_SH: project's spherical-harmonic series of the ellipsoid (shs_common), not
    pyshtools' Somigliana NormalGravity (which is evaluated on the geodetic latitude and
    would differ from the project's definition).
  * Threads: --threads sets OMP_NUM_THREADS. pyshtools' point evaluation is not expected to
    use it (the Fortran point routines are believed to be serial; not verified). The value is
    recorded in timings.csv regardless.
  * Timing: wall (time.perf_counter) and cpu (time.process_time) over the synthesis calls of
    V and grad V (one pass). Model load time is in notes as load_s=... .

All text in English. No gsh_tools import (libomp conflict with pyharm in the same process).
"""

import argparse
import os
import sys
import time


def parse_args(argv):
    p = argparse.ArgumentParser(description='pyshtools SHS benchmark runner')
    p.add_argument('--gfc', required=True, help='tide-free .gfc model')
    p.add_argument('--points', required=True, help='points CSV: id,lat_deg,lon_deg,h_m')
    p.add_argument('--nmax', required=True, type=int, help='maximum degree of the synthesis')
    p.add_argument('--out', required=True, help='result CSV, e.g. results/pyshtools_native_colorado.csv')
    p.add_argument('--timings', required=True, help='timings CSV (rows appended)')
    p.add_argument('--set', required=True, choices=['colorado', 'latsweep'], dest='point_set')
    p.add_argument('--repeat', type=int, default=1, help='number of timed repetitions')
    p.add_argument('--threads', type=int, default=1, help='OMP_NUM_THREADS for the process')
    return p.parse_args(argv)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    if args.repeat < 1 or args.threads < 1 or args.nmax < 0:
        raise SystemExit('repeat, threads must be >= 1 and nmax >= 0')

    # Must be set before the engine's native code is loaded.
    os.environ['OMP_NUM_THREADS'] = str(args.threads)

    import platform
    import numpy as np
    import pyshtools

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gfc_reader import read_gfc_header
    from shs_common import (FUNCTIONAL_TAG, read_points, geodetic_to_geocentric, disturbing_potential,
                            gravity_disturbance_mgal, write_result_csv, write_result_md, append_timing)
    from pyshtools.expand import MakeGridPoint
    from pyshtools.gravmag import MakeGravGridPoint

    nmax = args.nmax
    hdr = read_gfc_header(args.gfc)
    if nmax > hdr.max_degree:
        raise SystemExit(f'nmax={nmax} exceeds max_degree={hdr.max_degree} of the model')

    # ---- model load (timed separately, reported in notes) ----
    t0 = time.perf_counter()
    model = pyshtools.SHGravCoeffs.from_file(args.gfc, format='icgem', lmax=nmax,
                                             normalization='4pi', csphase=1)
    cilm = model.to_array(normalization='4pi', csphase=1, lmax=nmax)
    load_s = time.perf_counter() - t0
    gm, r0 = float(model.gm), float(model.r0)
    if abs(gm - hdr.gm) > 1e-9 * hdr.gm or abs(r0 - hdr.radius) > 1e-6:
        raise SystemExit(f'GM/R mismatch: pyshtools ({gm}, {r0}) vs header ({hdr.gm}, {hdr.radius})')
    if cilm.shape != (2, nmax + 1, nmax + 1):
        raise SystemExit(f'unexpected coefficient shape {cilm.shape}')

    # ---- points ----
    pts = read_points(args.points)
    r, phi_c, lam = geodetic_to_geocentric(pts['lat_deg'], pts['lon_deg'], pts['h_m'])
    lat_c_deg = np.degrees(phi_c)
    lon_deg = pts['lon_deg']
    n = len(r)

    # Group points by identical radius: the potential needs a radius-scaled copy of cilm.
    uniq_r, inverse = np.unique(r, return_inverse=True)
    order = np.argsort(inverse, kind='stable')
    counts = np.bincount(inverse, minlength=len(uniq_r))
    groups = np.split(order, np.cumsum(counts)[:-1])

    def synth():
        """Synthesis calls only: V of the model and its gradient (dV/dr, dV/dphi', dV/dlambda), one pass."""
        v = np.empty(n)
        for rk, idx in zip(uniq_r, groups):
            scale = (r0 / rk) ** np.arange(nmax + 1)
            scaled = cilm * scale[None, :, None]
            vals = MakeGridPoint(scaled, lat=lat_c_deg[idx], lon=lon_deg[idx], lmax=nmax,
                                 norm=1, csphase=1)
            v[idx] = gm / rk * np.asarray(vals, dtype=np.float64)
        b = np.empty((n, 3))
        for i in range(n):
            b[i] = MakeGravGridPoint(cilm, gm=gm, r0=r0, r=float(r[i]), lat=float(lat_c_deg[i]),
                                     lon=float(lon_deg[i]), lmax=nmax, omega=0.0)
        v_r = b[:, 0]
        v_lat = -r * b[:, 1]
        v_lon = r * np.cos(phi_c) * b[:, 2]
        return v, v_r, v_lat, v_lon

    # ---- timed repetitions ----
    last = None
    for rep in range(1, args.repeat + 1):
        w0 = time.perf_counter()
        c0 = time.process_time()
        last = synth()
        cpu = time.process_time() - c0
        wall = time.perf_counter() - w0
        notes = (f'load_s={load_s:.3f}; distinct_radii={len(uniq_r)}; '
                 f'threads_effect_unverified; platform={platform.system()}; functional={FUNCTIONAL_TAG}')
        append_timing(args.timings, {
            'engine': 'pyshtools', 'recurrence': 'native', 'set': args.point_set, 'nmax': nmax,
            'n_points': n, 'threads': args.threads, 'repeat': rep,
            'wall_s': f'{wall:.6f}', 'cpu_s': f'{cpu:.6f}',
            'notes': notes,
        })
        print(f'repeat {rep}: wall={wall:.3f}s cpu={cpu:.3f}s')

    command = 'python ' + ' '.join(['code/python/run_pyshtools.py'] + argv)
    common_params = [('nmax', nmax), ('set', args.point_set), ('repeat', args.repeat),
                     ('threads', args.threads), ('normalization', '4pi (fully normalised)'),
                     ('csphase', '1 (no Condon-Shortley phase)'),
                     ('functionals', 'V, T, dg from the same synthesis pass'),
                     ('gravity disturbance', '|g| - gamma_SH with centrifugal term, omega=7.292115e-5')]
    model_info = [('file', os.path.basename(args.gfc)), ('modelname', hdr.modelname),
                  ('GM [m^3/s^2]', f'{gm:.10e}'), ('R [m]', f'{r0:.6f}'),
                  ('tide_system', hdr.tide_system), ('load_s', f'{load_s:.3f}')]
    radii_note = ('pyshtools has no radius argument for the potential: coefficients are scaled per '
                  f'distinct radius ({len(uniq_r)} radii); this cost is inside the timed block.')

    v, v_r, v_lat, v_lon = last
    t_vals = disturbing_potential(v, r, phi_c)
    dg_vals = gravity_disturbance_mgal(v_r, v_lat, v_lon, r, phi_c)
    write_result_csv(args.out, pts, v, t_vals, dg_vals)

    write_result_md(
        os.path.splitext(args.out)[0] + '.md',
        csv_name=os.path.basename(args.out),
        engine_label=f'pyshtools {pyshtools.__version__} (native, single algorithm)',
        command=command,
        params=common_params,
        model_info=model_info,
        n_points=n,
        notes=[
            'V = GM/r * sum_n (R/r)^n S_n from MakeGridPoint on the radius-scaled coefficients: the model potential '
            'itself, no normal potential subtracted, degree 0 included (C00 as in the file).',
            'T = V - U, U = GRS80 level-ellipsoid gravitational potential (closed form), with the V of this run.',
            'Functionals evaluated at the actual point (geodetic h), not on the ellipsoid.',
            radii_note,
            'Gradient from MakeGravGridPoint (omega=0); centrifugal term added in code.',
        ],
    )

    print('first rows:')
    for i in range(min(3, n)):
        print(f"{pts['ids'][i]},{pts['lat_deg'][i]},{pts['lon_deg'][i]},{pts['h_m'][i]},"
              f"{v[i]:.15g},{t_vals[i]:.15g},{dg_vals[i]:.15g}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
