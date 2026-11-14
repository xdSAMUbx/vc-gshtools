"""pyshtools grid runner for the grid synthesis benchmark (see benchmark/grid/code/README.md, engine 'pyshtools').

Synthesises V and T on the Driscoll-Healy grids of the benchmark with pyshtools 4.14.1, using the SHTOOLS Fortran
backend explicitly:
  dh1   pyshtools' MakeGridDH(cilm, lmax=N, sampling=1), n x n with n = 2N+2
  dh2   pyshtools' MakeGridDH(cilm, lmax=N, sampling=2), n x 2n
Not run here (exit codes, same as shtools_grid_bench):
  regional   exit 3: pyshtools has no local grid synthesis
  gl         exit 4: MakeGridGLQ places 2*lmax+1 longitudes, not the contract grid with 2N+2 longitudes

Backend: pyshtools 4.14 prefers the ducc0 backend when it is installed. This runner calls the 'shtools' backend
module directly (pyshtools.backends.backend_module('shtools')), so the timed library is the Fortran SHTOOLS one.

Definitions (benchmark/shs/CONTRACT.md):
  V_m2s2 = GM/R * MakeGridDH(cilm) with fully normalised 4 pi coefficients, no Condon-Shortley phase, r = R (no
           radius scaling), degree 0 included (C00 as in the file).
  T_m2s2 = V - U, U = shs_common.normal_potential_closed(R, geocentric latitude).

Grid: row i = 0 is 90 deg, lat_i = 90 - i*180/n; lon_j = j*360/nlong (computed here from the grid dimensions, the
library's documented sampling of MakeGridDH).

Timing: time.perf_counter() around the MakeGridDH call only. U and T are computed after the timed runs.
Peak memory: resource.getrusage(RUSAGE_SELF).ru_maxrss in MB (KB on Linux/WSL); -1 where unavailable.

Output (never a whole grid):
  <out-dir>/pyshtools_<case>_rows.csv   header row,lat_deg,lon_deg,V_m2s2,T_m2s2 ; %.17g ; row 0-based
  <out-dir>/timings.csv                 appended: engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb
Row selection: as run_pyharm_grid.py (global rule; --rows overrides).

Usage:
  python code/python/run_pyshtools_grid.py --gfc <tide_free.gfc> --nmax N --case dh1|dh2|gl|regional --threads K
      --repeat R --out-dir <dir> [--rows i,j,k]

Library calls used (pyshtools 4.14.1):
  pyshtools.SHGravCoeffs.from_file(path, format='icgem', lmax=N, normalization='4pi', csphase=1)
  model.to_array(normalization='4pi', csphase=1, lmax=N), model.gm, model.r0
  pyshtools.backends.backend_module('shtools').MakeGridDH(cilm, lmax=N, norm=1, sampling=1|2, csphase=1)
    (keyword names as in the pyshtools backend documentation; the compiled module was not inspected, see README)

All text in English. Shared helpers come from benchmark/shs/code/python (gfc_reader.py, shs_common.py).
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SHS_PY = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'shs', 'code', 'python'))
sys.path.insert(0, SHS_PY)

CASES = ('gl', 'dh1', 'dh2', 'regional')
TIMINGS_HEADER = 'engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb'


def parse_args(argv):
    p = argparse.ArgumentParser(description='pyshtools grid synthesis runner')
    p.add_argument('--gfc', required=True, help='tide-free .gfc model')
    p.add_argument('--nmax', required=True, type=int, help='maximum degree of the synthesis')
    p.add_argument('--case', required=True, choices=CASES)
    p.add_argument('--threads', type=int, default=1, help='OMP_NUM_THREADS for the process (informational)')
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

    if args.case == 'regional':
        print('run_pyshtools_grid: case regional is not run: pyshtools has no local (regional) grid synthesis.',
              file=sys.stderr)
        return 3
    if args.case == 'gl':
        print('run_pyshtools_grid: case gl is not run: MakeGridGLQ gives 2*lmax+1 longitudes (spacing '
              '360/(2*lmax+1)), not the contract grid with 2*nmax+2 longitudes.', file=sys.stderr)
        return 4
    sampling = 1 if args.case == 'dh1' else 2

    os.environ['OMP_NUM_THREADS'] = str(args.threads)

    import numpy as np
    import pyshtools
    from gfc_reader import read_gfc_header
    from shs_common import normal_potential_closed

    hdr = read_gfc_header(args.gfc)
    nmax = args.nmax
    if nmax > hdr.max_degree:
        raise SystemExit(f'nmax={nmax} exceeds max_degree={hdr.max_degree} of the model')
    R = float(hdr.radius)

    model = pyshtools.SHGravCoeffs.from_file(args.gfc, format='icgem', lmax=nmax,
                                             normalization='4pi', csphase=1)
    cilm = model.to_array(normalization='4pi', csphase=1, lmax=nmax)
    gm, r0 = float(model.gm), float(model.r0)
    if abs(gm - hdr.gm) > 1e-9 * hdr.gm or abs(r0 - hdr.radius) > 1e-6:
        raise SystemExit(f'GM/R mismatch: pyshtools ({gm}, {r0}) vs header ({hdr.gm}, {hdr.radius})')
    if cilm.shape != (2, nmax + 1, nmax + 1):
        raise SystemExit(f'unexpected cilm shape {cilm.shape}')

    # The 'shtools' (Fortran) backend, not the preferred one (ducc0 when installed).
    make_grid_dh = pyshtools.backends.backend_module('shtools').MakeGridDH
    preferred = pyshtools.backends.preferred_backend()

    n = 2 * nmax + 2
    nlong = sampling * n
    lat_deg = [90.0 - i * 180.0 / n for i in range(n)]
    lon_deg = np.arange(nlong) * 360.0 / nlong
    phi_c = np.radians(lat_deg)
    u_row = np.asarray(normal_potential_closed(R, phi_c), dtype=float)

    timings_path = os.path.join(args.out_dir, 'timings.csv')
    grid = None
    for rep in range(1, args.repeat + 1):
        t0 = time.perf_counter()
        grid = make_grid_dh(cilm, lmax=nmax, norm=1, sampling=sampling, csphase=1)
        wall = time.perf_counter() - t0
        grid = np.asarray(grid, dtype=float)
        if grid.shape != (n, nlong):
            raise SystemExit(f'unexpected MakeGridDH shape {grid.shape}, expected {(n, nlong)}')
        append_timing(timings_path, ['pyshtools', args.case, 'native', args.threads, nmax, n, nlong, rep,
                                     f'{wall:.9g}', f'{peak_mem_mb():.3f}'])
        print(f'repeat {rep}: wall={wall:.6f} s')

    if not np.all(np.isfinite(grid)):
        raise SystemExit('non-finite values in the grid')
    v = (gm / R) * grid
    t = v - u_row[:, None]
    rows = select_rows(args, lat_deg)
    rows_path = os.path.join(args.out_dir, f'pyshtools_{args.case}_rows.csv')
    write_rows(rows_path, rows, lat_deg, lon_deg, v, t)
    print(f'run_pyshtools_grid: {args.case} nmax={nmax} nlat={n} nlon={nlong} backend={preferred} '
          f'pyshtools {pyshtools.__version__} (MakeGridDH from the shtools backend), rows written: {len(rows)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
