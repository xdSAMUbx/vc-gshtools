function out = run_graflab_grid(gfc_path, nmax, grid_case, fnALFs, out_dir, nrep, rows, graflab_dir)
%RUN_GRAFLAB_GRID GrafLab 2.1.4 grid mode (point_type = 0) for the regional case of the grid benchmark.
%
%   out = run_graflab_grid(gfc_path, nmax, grid_case, fnALFs, out_dir, nrep, rows, graflab_dir)
%
%   gfc_path     tide-free .gfc model (the same file as the SHS benchmark)
%   nmax         maximum degree of the expansion
%   grid_case    must be 'regional'. The global cases (gl, dh1, dh2) are refused: 9.6 million nodes at degree
%                2190 need prohibitive memory and hours in GrafLab.
%   fnALFs       1 standard forward column (sfc), 2 modified forward column (msfc), 3 extended range (era)
%   out_dir      existing directory for the rows CSV and timings.csv
%   nrep         number of timed repetitions (default 1)
%   rows         0-based row indices to write (default [0 50 100 150 200]); row 0 is 45 deg north
%   graflab_dir  folder that contains GrafLab.m (not shipped with this repository; see README.md)
%
%   Grid (benchmark definition): latitude 45 -> 35 deg (north to south), longitude -110 -> -100 deg, step 0.05 deg,
%   201 x 201 nodes, spherical (geocentric) latitude at the constant radius R of the model header (coord = 1, h = 0).
%
%   Functionals: V (functional 11, GrafLab takes the file's coefficients directly, C00 included) and T
%   (functional 5, GrafLab's disturbing potential with its own normal field). Note that T here is GrafLab's
%   functional 5, not the closed-form U of shs_common.py used by the other engines; the SHS benchmark already
%   compares the two at the 1e-6 m^2/s^2 level (see benchmark/shs/code/python/shs_common.py).
%
%   Output: GrafLab returns, for nargout == 1 in grid mode, the matrix [lat, lon, V, T] with one row per node
%   (GrafLab.m, grid branch: output = [fi; lambda; P] transposed). The script places each row by its own
%   latitude and longitude, so the node order of GrafLab is not assumed.
%
%   Files written (never a whole grid):
%     <out_dir>/graflab_regional_<recurrence>_rows.csv   header row,lat_deg,lon_deg,V_m2s2,T_m2s2 (%.17g)
%     <out_dir>/timings.csv   appended: engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb
%
%   Peak memory is not measured here: peak_mem_mb is written as -1 (see README.md).
%   Written but NOT run by its author. GrafLab was not executed for this file.

if nargin < 6 || isempty(nrep), nrep = 1; end
if nargin < 7 || isempty(rows), rows = [0 50 100 150 200]; end
if nargin < 8 || isempty(graflab_dir)
    error('run_graflab_grid:noGrafLab', 'Give the folder that contains GrafLab.m as the eighth argument.');
end

if ~strcmp(grid_case, 'regional')
    error('run_graflab_grid:globalGrid', ['GrafLab is not run on the global grids (case "%s"): 9.6 million nodes ' ...
        'at degree 2190 need prohibitive memory and hours. Use grid_case = ''regional''.'], grid_case);
end
if ~any(fnALFs == [1 2 3])
    error('run_graflab_grid:fnALFs', 'fnALFs must be 1 (sfc), 2 (msfc) or 3 (era).');
end
rec_names = {'sfc', 'msfc', 'era'};
recurrence = rec_names{fnALFs};
if nrep < 1, error('run_graflab_grid:nrep', 'nrep must be >= 1.'); end

here = fileparts(mfilename('fullpath'));
addpath(graflab_dir);
addpath(fullfile(here, '..', '..', '..', 'shs', 'code', 'matlab'));   % read_gfc_header.m (shared with the SHS runs)

hdr = read_gfc_header(gfc_path);
if nmax > hdr.max_degree
    error('run_graflab_grid:nmax', 'nmax=%d exceeds the degree of the model (%d).', nmax, hdr.max_degree);
end

% ---- GrafLab('OK', ...) argument list: same order as graflab_evaluate.m (GrafLab.m, function header) ----
GM = hdr.GM;               % GrafLab reads GM and R from the .gfc header anyway
R = hdr.R;
nmin = 0;
ellipsoid = 1;             % GRS80, as in the SHS benchmark
coord = 1;                 % spherical input coordinates: latitude is the geocentric one, radius R + h
point_type = 0;            % grid mode
lat_min = 35; lat_step = 0.05; lat_max = 45;     % scalars: GrafLab builds lat_min:lat_step:lat_max
lon_min = -110; lon_step = 0.05; lon_max = -100;
h = 0;                     % height above the sphere of radius R (m)
Input_data_path = [];
lat = []; lon = []; h2 = [];
Output_path = fullfile(tempdir, 'graflab_grid_unused_prefix');   % not written: exports are disabled
Functional_or_Commission = 0;
Functional = [11; 5];      % V, T (output columns in this order)
DTM_path = [];
Export_data_txt = 0;
Export_report = 0;
Export_data_mat = 0;
Display_data = 0;
Graphic_format = [];
Colormap = [];
Number_of_colors = [];
DPI = [];
Status_bar = 0;

nlat = round((lat_max - lat_min) / lat_step) + 1;   % 201
nlon = round((lon_max - lon_min) / lon_step) + 1;   % 201

timings_path = fullfile(out_dir, 'timings.csv');
maxthreads = maxNumCompThreads;
res = [];
for rep = 1:nrep
    tic;
    res = GrafLab('OK', GM, R, nmin, nmax, ellipsoid, gfc_path, coord, point_type, ...
        lat_min, lat_step, lat_max, lon_min, lon_step, lon_max, h, Input_data_path, ...
        lat, lon, h2, Output_path, Functional_or_Commission, Functional, fnALFs, DTM_path, ...
        Export_data_txt, Export_report, Export_data_mat, ...
        Display_data, Graphic_format, Colormap, Number_of_colors, DPI, Status_bar);
    wall = toc;
    append_timing(timings_path, recurrence, maxthreads, nmax, nlat, nlon, rep, wall);
    fprintf('run_graflab_grid: repeat %d, recurrence %s, wall %.3f s\n', rep, recurrence, wall);
end

% ---- place the nodes by their own coordinates (row 1 = 45 deg north, column 1 = -110 deg) ----
if size(res, 2) ~= 4 || size(res, 1) ~= nlat * nlon
    error('run_graflab_grid:shape', 'unexpected GrafLab output size %d x %d (expected %d x 4).', ...
        size(res, 1), size(res, 2), nlat * nlon);
end
i_idx = round((lat_max - res(:, 1)) / lat_step) + 1;
j_idx = round((res(:, 2) - lon_min) / lon_step) + 1;
if any(i_idx < 1 | i_idx > nlat | j_idx < 1 | j_idx > nlon)
    error('run_graflab_grid:index', 'GrafLab node outside the requested grid.');
end
Vg = nan(nlat, nlon);
Tg = nan(nlat, nlon);
lin = sub2ind([nlat nlon], i_idx, j_idx);
Vg(lin) = res(:, 3);
Tg(lin) = res(:, 4);
if any(isnan(Vg(:))) || any(isnan(Tg(:)))
    error('run_graflab_grid:missing', 'some grid nodes were not returned by GrafLab.');
end

rows = sort(unique(rows(:)'));
if any(rows < 0 | rows > nlat - 1)
    error('run_graflab_grid:rows', 'rows must be in [0, %d].', nlat - 1);
end
rows_path = fullfile(out_dir, sprintf('graflab_regional_%s_rows.csv', recurrence));
fid = fopen(rows_path, 'w');
if fid < 0, error('run_graflab_grid:io', 'cannot write %s (does out_dir exist?)', rows_path); end
fprintf(fid, 'row,lat_deg,lon_deg,V_m2s2,T_m2s2\n');
for k = 1:numel(rows)
    i = rows(k);                                  % 0-based row
    la = lat_max - i * lat_step;
    for j = 0:nlon-1
        lo = lon_min + j * lon_step;
        fprintf(fid, '%d,%.17g,%.17g,%.17g,%.17g\n', i, la, lo, Vg(i + 1, j + 1), Tg(i + 1, j + 1));
    end
end
fclose(fid);

out = struct('V', Vg, 'T', Tg, 'recurrence', recurrence, 'rows_file', rows_path);
end

function append_timing(path, recurrence, threads, nmax, nlat, nlon, rep, wall)
need_header = ~isfile(path);
fid = fopen(path, 'a');
if fid < 0, error('run_graflab_grid:timings', 'cannot append to %s (does out_dir exist?)', path); end
if need_header
    fprintf(fid, 'engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb\n');
end
fprintf(fid, 'graflab,regional,%s,%d,%d,%d,%d,%d,%.9g,-1\n', recurrence, threads, nmax, nlat, nlon, rep, wall);
fclose(fid);
end
