"""pyharm grid runner for the grid synthesis benchmark (see benchmark/grid/code/README.md, engine 'pyharm').

Synthesises V and T on the grids of the benchmark with pyharm 0.4.11 (a Python wrapper of CHarm 0.4.11):
  gl, dh1, dh2   pyharm.crd.PointGridGL / PointGridDH1 / PointGridDH2 (nmax, r = R), the CHarm constructors
  regional       pyharm.crd.PointGrid.from_arrays(lat, lon, r) on the 201 x 201 grid (lat 45 -> 35 deg, lon -110 ->
                 -100 deg, step 0.05 deg), with constant r = R
  pyharm.shs.point(pnt, shcs, nmax) returns V on the grid as an array of shape (nlat, nlon).

Definitions (benchmark/shs/CONTRACT.md):
  V_m2s2 = GM/R sum_n sum_m (...) with r = R, degree 0 included (C00 as in the file); no normal potential subtracted.
  T_m2s2 = V - U, U = shs_common.normal_potential_closed(R, geocentric latitude) (closed form; depends on latitude only).

Timing: the wall time of pyharm.shs.point() only (time.perf_counter). U and T are computed after the timed runs.
Peak memory: resource.getrusage(RUSAGE_SELF).ru_maxrss (KB on Linux/WSL), in MB; -1 where the resource module is
not available (Windows). The peak is cumulative over the process.

Output (never a whole grid):
  <out-dir>/pyharm_<case>_rows.csv   header row,lat_deg,lon_deg,V_m2s2,T_m2s2 ; values with %.17g ; row is 0-based,
                                     0 = northernmost
  <out-dir>/timings.csv              appended: engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb

Row selection: global grids use the rule of the benchmark (rows closest to 0, 30, 60, 85, 89.5 deg north and to
-30, -60, -85, -89.5 deg, plus rows 0 and nlat-1); regional default rows 0, 50, 100, 150, 200; --rows overrides.

Usage:
  python code/python/run_pyharm_grid.py --gfc <tide_free.gfc> --nmax N --case gl|dh1|dh2|regional --threads K
      --repeat R --out-dir <dir> [--rows i,j,k]

Library calls used (pyharm 0.4.11, public names):
  pyharm.shc.Shc.from_file('gfc', path, nmax=N), shcs.mu, shcs.r
  pyharm.crd.PointGridGL(nmax, r=R), pyharm.crd.PointGridDH1(nmax, r=R), pyharm.crd.PointGridDH2(nmax, r=R)
  pyharm.crd.PointGrid.from_arrays(lat, lon, r), pnt.lat, pnt.lon (radians)
  pyharm.shs.point(pnt, shcs, nmax)

All text in English. Shared helpers are imported from benchmark/shs/code/python (gfc_reader.py, shs_common.py).
"""
import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SHS_PY = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'shs', 'code', 'python'))
sys.path.insert(0, SHS_PY)

REG_N = 201
REG_LAT_NORTH = 45.0
REG_LAT_SOUTH = 35.0
REG_LON_WEST = -110.0
REG_LON_EAST = -100.0
CASES = ('gl', 'dh1', 'dh2', 'regional')
TIMINGS_HEADER = 'engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb'


def parse_args(argv):
    p = argparse.ArgumentParser(description='pyharm grid synthesis runner')
    p.add_argument('--gfc', required=True, help='tide-free .gfc model')
    p.add_argument('--nmax', required=True, type=int, help='maximum degree of the synthesis')
    p.add_argument('--case', required=True, choices=CASES)
    p.add_argument('--threads', type=int, default=1, help='OMP_NUM_THREADS for the process')
    p.add_argument('--repeat', type=int, default=1, help='number of timed repetitions')
    p.add_argument('--out-dir', required=True, help='directory for the rows CSV and timings.csv (must exist)')
    p.add_argument('--rows', default=None, help='comma-separated 0-based row indices (overrides the default)')
    return p.parse_args(argv)


def peak_mem_mb():
    try:
        import resource
    except ImportError:
        return -1.0
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def recurrence_label(lon_rad, nmax):
    """Expectation of the CHarm path, from the conditions documented in charm_shs.h (not read back from CHarm):
    'fft' when nlon > 1, (nlon-1)//2 >= nmax, constant dlon, lon[0] == 0 and lon[-1] + dlon == 2 pi; else 'chebyshev'.
    """
    thr = 1e-9
    nlon = len(lon_rad)
    if nlon < 2:
        return 'chebyshev'
    dlon = lon_rad[1] - lon_rad[0]
    ok = ((nlon - 1) // 2 >= nmax and abs(lon_rad[0]) < thr
          and abs(lon_rad[-1] + dlon - 2.0 * math.pi) < thr)
    if ok:
        ok = all(abs((lon_rad[j] - lon_rad[j - 1]) - dlon) <= thr for j in range(2, nlon))
    return 'fft' if ok else 'chebyshev'


def nearest_row(lat_deg, target):
    return min(range(len(lat_deg)), key=lambda i: abs(lat_deg[i] - target))


def select_rows(args, lat_deg):
    nlat = len(lat_deg)
    if args.rows is not None:
        rows = set()
        for tok in args.rows.split(','):
            v = int(tok)
            if not 0 <= v < nlat:
                raise SystemExit('--rows: index out of range')
            rows.add(v)
        return sorted(rows)
    if args.case == 'regional':
        return sorted({0, (nlat - 1) // 4, (nlat - 1) // 2, 3 * (nlat - 1) // 4, nlat - 1})
    rows = {nearest_row(lat_deg, t) for t in (0.0, 30.0, 60.0, 85.0, 89.5, -30.0, -60.0, -85.0, -89.5)}
    rows.update({0, nlat - 1})
    return sorted(rows)


def append_timing(path, row):
    need_header = not os.path.exists(path)
    with open(path, 'a', encoding='utf-8', newline='') as f:
        if need_header:
            f.write(TIMINGS_HEADER + '\n')
        f.write(','.join(str(x) for x in row) + '\n')


def write_rows(path, rows, lat_deg, lon_deg, v, t):
    """Writes the selected rows only; v and t are numpy arrays of shape (nlat, nlon)."""
    with open(path, 'w', encoding='utf-8', newline='') as f:
        f.write('row,lat_deg,lon_deg,V_m2s2,T_m2s2\n')
        for i in rows:
            la = f'{lat_deg[i]:.17g}'
            for j in range(len(lon_deg)):
                f.write(f'{i},{la},{lon_deg[j]:.17g},{v[i, j]:.17g},{t[i, j]:.17g}\n')


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.repeat < 1 or args.threads < 1 or args.nmax < 0:
        raise SystemExit('repeat and threads must be >= 1 and nmax >= 0')
    # Must be set before CHarm's OpenMP runtime is loaded by "import pyharm".
    os.environ['OMP_NUM_THREADS'] = str(args.threads)

    import numpy as np
    from gfc_reader import read_gfc_header
    from shs_common import normal_potential_closed
    import pyharm as ph

    hdr = read_gfc_header(args.gfc)
    nmax = args.nmax
    if nmax > hdr.max_degree:
        raise SystemExit(f'nmax={nmax} exceeds max_degree={hdr.max_degree} of the model')
    R = float(hdr.radius)

    shcs = ph.shc.Shc.from_file('gfc', args.gfc, nmax=nmax)
    gm, r0 = float(shcs.mu), float(shcs.r)
    if abs(gm - hdr.gm) > 1e-9 * hdr.gm or abs(r0 - hdr.radius) > 1e-6:
        raise SystemExit(f'GM/R mismatch: pyharm ({gm}, {r0}) vs header ({hdr.gm}, {hdr.radius})')

    if args.case == 'gl':
        pnt = ph.crd.PointGridGL(nmax, r=R)
    elif args.case == 'dh1':
        pnt = ph.crd.PointGridDH1(nmax, r=R)
    elif args.case == 'dh2':
        pnt = ph.crd.PointGridDH2(nmax, r=R)
    else:
        k = np.arange(REG_N)
        lat = np.radians(REG_LAT_NORTH - (REG_LAT_NORTH - REG_LAT_SOUTH) * k / (REG_N - 1))
        lon = np.radians(REG_LON_WEST + (REG_LON_EAST - REG_LON_WEST) * k / (REG_N - 1))
        rr = np.full(REG_N, R)
        pnt = ph.crd.PointGrid.from_arrays(np.ascontiguousarray(lat), np.ascontiguousarray(lon),
                                           np.ascontiguousarray(rr))

    lat_rad = np.asarray(pnt.lat, dtype=float)
    lon_rad = np.asarray(pnt.lon, dtype=float)
    nlat, nlon = lat_rad.size, lon_rad.size
    if not lat_rad[0] > lat_rad[-1]:
        raise SystemExit('grid is not ordered north to south (row 0 must be the northernmost)')
    recurrence = recurrence_label(lon_rad, nmax)
    lat_deg = np.degrees(lat_rad).tolist()
    lon_deg = np.degrees(lon_rad)
    u_row = np.asarray(normal_potential_closed(R, lat_rad), dtype=float)

    timings_path = os.path.join(args.out_dir, 'timings.csv')
    v = None
    for rep in range(1, args.repeat + 1):
        t0 = time.perf_counter()
        v = ph.shs.point(pnt, shcs, nmax)
        wall = time.perf_counter() - t0
        append_timing(timings_path, ['pyharm', args.case, recurrence, args.threads, nmax, nlat, nlon, rep,
                                     f'{wall:.9g}', f'{peak_mem_mb():.3f}'])
        print(f'repeat {rep}: wall={wall:.6f} s, recurrence={recurrence}')

    v = np.asarray(v, dtype=float)
    if v.shape != (nlat, nlon) or not np.all(np.isfinite(v)):
        raise SystemExit(f'unexpected V shape {v.shape} or non-finite values')
    t = v - u_row[:, None]
    rows = select_rows(args, lat_deg)
    rows_path = os.path.join(args.out_dir, f'pyharm_{args.case}_rows.csv')
    write_rows(rows_path, rows, lat_deg, lon_deg, v, t)
    print(f'run_pyharm_grid: {args.case} nmax={nmax} nlat={nlat} nlon={nlon} recurrence={recurrence} '
          f'pyharm {ph.misc.get_version()}, rows written: {len(rows)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
