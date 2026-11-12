"""pyharm runner for the SHS benchmark (definitions: benchmark/shs/results/csv.md, engine 'pyharm').

One run evaluates at each point the three functionals defined in results/csv.md from the tide-free model, in the same
synthesis pass, and writes one result CSV:
  * V_m2s2 = pyharm.shs.point(...)     (gravitational potential of the model, mu/r sum (R/r)^n ...; degree 0 included,
                                        C00 as in the file)
  * T_m2s2 = V - U                     (disturbing potential, m^2/s^2; V is the value above, not a second synthesis)
  * dg_mgal = |g| - gamma_SH           (gravity disturbance, mGal; from shs.point_grad1 of the same run)

Usage:
  python code/python/run_pyharm.py --gfc <tide_free.gfc> --points <csv> --nmax N
      --out results/pyharm_native_<set>.csv --timings results/timings.csv
      --set colorado|latsweep --repeat R --threads K

Library calls used (pyharm 0.4.11, all public names):
  * pyharm.shc.Shc.from_file('gfc', path, nmax=N)  -> model reader (CHarm's GFC reader).
    The Shc carries mu (= GM) and r (= R) from the file; the script checks them against the
    header (gfc_reader.py).
  * pyharm.crd.PointSctr.from_arrays(lat_geocentric_rad, lon_rad, r_m) -> scattered points.
  * pyharm.shs.point(pnt, shcs, nmax)          -> potential V = mu/r sum (R/r)^n ... at each point.
  * pyharm.shs.point_grad1(pnt, shcs, nmax)    -> three LNOF components of grad V, returned as
    (f[P], f[L], f[R]) = (V_phi'/r, -V_lambda/(r cos phi'), V_r). Order and signs were verified
    against finite differences of V on a synthetic model (scratchpad test): all three agree
    with the finite differences within the finite-difference accuracy (about 1e-5 relative
    for the two tangential components, with step 1e-7 rad).

Conventions and how each difference was resolved:
  * Normalisation: fully normalised 4pi (geodesy). Checked against an independent numpy
    recursion on a synthetic model: exact agreement (~1e-16 relative).
  * Condon-Shortley phase: none. CHarm's shs uses the GFC convention. Verified as above.
  * Units: CHarm works in the units of the file (m, m^3/s^2). mu and r are taken from the
    file. Output V is m^2/s^2 directly, no rescaling.
  * Point radius: pyharm evaluates at a per-point radius natively (solid harmonic
    synthesis), so no coefficient scaling is needed.
  * Gradient: point_grad1 returns the gradient in LNOF with the 1/r factor on the latitude
    component (see above). The script converts back to partial derivatives of V:
    dV/dphi' = r f[P], dV/dlambda = -r cos(phi') f[L], dV/dr = f[R]. The centrifugal term
    and the normal gravity are added in shs_common.gravity_disturbance_mgal.
  * Geometry: geodetic (GRS80) -> geocentric (phi', r) is done exactly in shs_common.
    pyharm receives the geocentric latitude and longitude in radians.
  * Ellipsoidal handling: V, T and dg at the actual point (r includes h), as in the project's
    functionals.cpp (dg uses the point's R(i)).
  * U: closed form of the GRS80 level-ellipsoid gravitational potential (shs_common).
  * gamma_SH: project's spherical-harmonic series of the ellipsoid (shs_common).
  * Legendre algorithm: pyharm's shs point synthesis exposes one algorithm (no recurrence
    option), so there is a single run, labelled 'native'. The leg module only changes the
    storage ordering of Fourier coefficients (PMNJ/PMJN) and is not used by shs.point.
  * Threads: OMP_NUM_THREADS is set from --threads before pyharm is imported. Whether
    CHarm's shs point synthesis actually uses OpenMP was NOT verified: at nmax 360 with 50
    points, 1 and 4 threads gave the same wall time within noise.
  * Timing: wall (time.perf_counter) and cpu (time.process_time) over the synthesis calls
    shs.point and shs.point_grad1 (one pass). Model load time is in notes as load_s=... .

All text in English. No gsh_tools import (libomp conflict with pyharm in the same process).
"""

import argparse
import os
import sys
import time


def parse_args(argv):
    p = argparse.ArgumentParser(description='pyharm SHS benchmark runner')
    p.add_argument('--gfc', required=True, help='tide-free .gfc model')
    p.add_argument('--points', required=True, help='points CSV: id,lat_deg,lon_deg,h_m')
    p.add_argument('--nmax', required=True, type=int, help='maximum degree of the synthesis')
    p.add_argument('--out', required=True, help='result CSV, e.g. results/pyharm_native_colorado.csv')
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

    # Must be set before CHarm's OpenMP runtime is loaded by "import pyharm".
    os.environ['OMP_NUM_THREADS'] = str(args.threads)

    import platform
    import numpy as np

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gfc_reader import read_gfc_header
    from shs_common import (FUNCTIONAL_TAG, read_points, geodetic_to_geocentric, disturbing_potential,
                            gravity_disturbance_mgal, write_result_csv, write_result_md, append_timing)
    import pyharm as ph

    nmax = args.nmax
    hdr = read_gfc_header(args.gfc)
    if nmax > hdr.max_degree:
        raise SystemExit(f'nmax={nmax} exceeds max_degree={hdr.max_degree} of the model')

    # ---- model load (timed separately, reported in notes) ----
    t0 = time.perf_counter()
    shcs = ph.shc.Shc.from_file('gfc', args.gfc, nmax=nmax)
    load_s = time.perf_counter() - t0
    gm, r0 = float(shcs.mu), float(shcs.r)
    if abs(gm - hdr.gm) > 1e-9 * hdr.gm or abs(r0 - hdr.radius) > 1e-6:
        raise SystemExit(f'GM/R mismatch: pyharm ({gm}, {r0}) vs header ({hdr.gm}, {hdr.radius})')

    # ---- points ----
    pts = read_points(args.points)
    r, phi_c, lam = geodetic_to_geocentric(pts['lat_deg'], pts['lon_deg'], pts['h_m'])
    n = len(r)
    pnt = ph.crd.PointSctr.from_arrays(np.ascontiguousarray(phi_c),
                                       np.ascontiguousarray(lam),
                                       np.ascontiguousarray(r))

    def synth():
        """Synthesis calls only: potential V of the model and the LNOF gradient, one pass."""
        v = ph.shs.point(pnt, shcs, nmax)
        g_p, g_l, g_r = ph.shs.point_grad1(pnt, shcs, nmax)
        return v, g_p, g_l, g_r

    # ---- timed repetitions ----
    last = None
    for rep in range(1, args.repeat + 1):
        w0 = time.perf_counter()
        c0 = time.process_time()
        last = synth()
        cpu = time.process_time() - c0
        wall = time.perf_counter() - w0
        notes = (f'load_s={load_s:.3f}; shs.point+point_grad1; platform={platform.system()}; '
                 f'functional={FUNCTIONAL_TAG}')
        append_timing(args.timings, {
            'engine': 'pyharm', 'recurrence': 'native', 'set': args.point_set, 'nmax': nmax,
            'n_points': n, 'threads': args.threads, 'repeat': rep,
            'wall_s': f'{wall:.6f}', 'cpu_s': f'{cpu:.6f}',
            'notes': notes,
        })
        print(f'repeat {rep}/{args.repeat}: wall={wall:.3f}s cpu={cpu:.3f}s')

    command = 'python ' + ' '.join(['code/python/run_pyharm.py'] + argv)
    v, g_p, g_l, g_r = last
    v = np.asarray(v, dtype=np.float64)
    v_r = np.asarray(g_r, dtype=np.float64)
    v_lat = np.asarray(g_p, dtype=np.float64) * r
    v_lon = -np.asarray(g_l, dtype=np.float64) * r * np.cos(phi_c)
    t_vals = disturbing_potential(v, r, phi_c)
    dg_vals = gravity_disturbance_mgal(v_r, v_lat, v_lon, r, phi_c)
    write_result_csv(args.out, pts, v, t_vals, dg_vals)

    write_result_md(
        os.path.splitext(args.out)[0] + '.md',
        csv_name=os.path.basename(args.out),
        engine_label=f'pyharm {ph.misc.get_version()} (native, single algorithm)',
        command=command,
        params=[('nmax', nmax), ('set', args.point_set), ('repeat', args.repeat),
                ('threads', args.threads), ('normalization', '4pi (fully normalised)'),
                ('csphase', 'none (GFC convention)'),
                ('functionals', 'V, T, dg from the same synthesis pass'),
                ('gravity disturbance', '|g| - gamma_SH with centrifugal term, omega=7.292115e-5')],
        model_info=[('file', os.path.basename(args.gfc)), ('modelname', hdr.modelname),
                    ('GM [m^3/s^2]', f'{gm:.10e}'), ('R [m]', f'{r0:.6f}'),
                    ('tide_system', hdr.tide_system), ('load_s', f'{load_s:.3f}')],
        n_points=n,
        notes=[
            'V = pyharm.shs.point: the model potential itself, no normal potential subtracted, degree 0 included '
            '(C00 as in the file).',
            'T = V - U, U = GRS80 level-ellipsoid gravitational potential (closed form), with the V of this run.',
            'Functionals evaluated at the actual point (geodetic h), not on the ellipsoid.',
            'Gradient from shs.point_grad1 (LNOF, latitude component divided by r); '
            'converted back to partial derivatives of V.',
        ],
    )

    print('first rows:')
    for i in range(min(3, n)):
        print(f"{pts['ids'][i]},{pts['lat_deg'][i]},{pts['lon_deg'][i]},{pts['h_m'][i]},"
              f"{v[i]:.15g},{t_vals[i]:.15g},{dg_vals[i]:.15g}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
