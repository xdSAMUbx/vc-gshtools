function run_graflab(gfc_path, points_csv, nmax, recurrence, set_name, out_dir, repeat_count, graflab_dir)
%RUN_GRAFLAB Synthesise V (gravitational potential), T (disturbing potential) and dg (gravity disturbance) in one
%   GrafLab 2.1.4 call per repetition (headless, point-wise mode). One result CSV per recurrence and point set.
%
%   run_graflab(gfc_path, points_csv, nmax, recurrence, set_name, out_dir, repeat_count, graflab_dir)
%
%   gfc_path      Tide-free ICGEM .gfc file (its name must end in ".gfc").
%   points_csv    CSV with header id,lat_deg,lon_deg,h_m (GRS80 geodetic latitude,
%                 ellipsoidal height in metres).
%   nmax          Maximum degree of the expansion (integer, <= degree of the .gfc).
%   recurrence    'sfc'  standard forward column   (GrafLab fnALFs = 1)
%                 'msfc' modified forward column   (GrafLab fnALFs = 2)
%                 'era'  extended-range arithmetic (GrafLab fnALFs = 3)
%   set_name      Point-set label used in file names, e.g. 'colorado' or 'latsweep'.
%   out_dir       Results folder (benchmark/shs/results). Receives
%                 graflab_<recurrence>_<set_name>.csv (id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal),
%                 its .md companion, and appends rows to timings.csv.
%   repeat_count  Number of timed repetitions. Each repetition is a full GrafLab
%                 call (model load + synthesis); see graflab_evaluate.m.
%   graflab_dir   Folder that contains GrafLab.m. GrafLab is NOT shipped in this
%                 repository (see code/matlab/README.md). Optional: if omitted,
%                 the environment variable GRAFLAB_DIR is used.
%
%   One call computes the three functionals: GrafLab Functional = [11; 5; 20] (V, T, dg), the output columns
%   are V, T and dg in that order (see graflab_evaluate.m).
%
%   Example (from PowerShell):
%     matlab -batch "cd('C:\...\benchmark\shs\code\matlab'); run_graflab('C:\...\XGM2019e_2190_tide_free.gfc', 'C:\...\points_colorado.csv', 360, 'sfc', 'colorado', 'C:\...\benchmark\shs\results', 1, 'C:\...\basura')"

if nargin < 8 || isempty(graflab_dir)
    graflab_dir = getenv('GRAFLAB_DIR');
end
n_functionals = 3;                    % V, T, dg
value_cols = 4:(3 + n_functionals);   % columns of GrafLab output holding the functional values

if isempty(graflab_dir) || ~isfile(fullfile(graflab_dir, 'GrafLab.m'))
    error('run_graflab:noGrafLab', ...
        'GrafLab.m not found in "%s". Pass its folder as the 8th argument or set GRAFLAB_DIR. See code/matlab/README.md.', ...
        graflab_dir);
end
graflab_dir = to_absolute(graflab_dir);

gfc_path = to_absolute(gfc_path);
points_csv = to_absolute(points_csv);
out_dir = to_absolute(out_dir);
if ~isfile(gfc_path)
    error('run_graflab:noGfc', 'Geopotential model not found: %s', gfc_path);
end
if ~isfile(points_csv)
    error('run_graflab:noPoints', 'Points file not found: %s', points_csv);
end

switch recurrence
    case 'sfc'
        fnALFs = 1;
    case 'msfc'
        fnALFs = 2;
    case 'era'
        fnALFs = 3;
    otherwise
        error('run_graflab:badRecurrence', 'recurrence must be sfc, msfc or era (got "%s").', recurrence);
end
validateattributes(nmax, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(repeat_count, {'numeric'}, {'scalar', 'integer', '>=', 1});
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

addpath(graflab_dir);

hdr = read_gfc_header(gfc_path);
if nmax > hdr.max_degree
    error('run_graflab:nmaxTooLarge', 'nmax=%d exceeds the degree of the model (%d).', nmax, hdr.max_degree);
end

pts = read_points_csv(points_csv);
n_points = numel(pts.lat);

p = struct('gfc_path', gfc_path, 'GM', hdr.GM, 'R', hdr.R, 'nmax', nmax, ...
    'fnALFs', fnALFs, 'lat', pts.lat, 'lon', pts.lon, 'h', pts.h);

% Model load time, measured with a separate 1-point call through the same
% GrafLab pipeline (header parse + textscan of the whole .gfc + 1 point).
% An untimed warm-up call comes first, so the probe is not a cold-start call
% (a cold probe made the first subtraction negative in the smoke test).
p1 = p;
p1.lat = pts.lat(1); p1.lon = pts.lon(1); p1.h = pts.h(1);
graflab_evaluate(p1);
c0 = cputime; t0 = tic;
graflab_evaluate(p1);
load_s = toc(t0);
load_cpu = cputime - c0;
fprintf('run_graflab: model load (warm, 1-point call) = %.3f s wall, %.3f s cpu\n', load_s, load_cpu);

threads = maxNumCompThreads;
timings_csv = fullfile(out_dir, 'timings.csv');
notes_base = sprintf('load_s=%.3f load_cpu_s=%.3f; wall_s = total_call_wall_s - load_s (load from separate 1-point call); tide_system=%s; GrafLab 2.1.4; MATLAB %s; fnALFs=%d; functional=V_T_dg', ...
    load_s, load_cpu, hdr.tide_system, version('-release'), fnALFs);

out_first = [];
max_diff = zeros(repeat_count, n_functionals);
for r = 1:repeat_count
    c0 = cputime; t0 = tic;
    out = graflab_evaluate(p);
    total_wall = toc(t0);
    total_cpu = cputime - c0;

    if size(out, 1) ~= n_points || size(out, 2) ~= 3 + n_functionals
        error('run_graflab:badOutput', 'Unexpected GrafLab output size [%d %d], expected [%d %d].', ...
            size(out, 1), size(out, 2), n_points, 3 + n_functionals);
    end
    if any(~isfinite(out(:, value_cols)))
        error('run_graflab:nonFinite', 'GrafLab returned non-finite V, T or dg values.');
    end

    if r == 1
        out_first = out;
    else
        max_diff(r, :) = max(abs(out(:, value_cols) - out_first(:, value_cols)), [], 1);
    end

    wall_s = total_wall - load_s;
    cpu_s = total_cpu - load_cpu;
    notes = sprintf('%s; total_call_wall_s=%.3f', notes_base, total_wall);
    append_timing_row(timings_csv, 'graflab', recurrence, set_name, nmax, n_points, threads, r, wall_s, cpu_s, notes);
    fprintf('run_graflab: repeat %d/%d  total=%.3f s  synthesis(wall)=%.3f s\n', r, repeat_count, total_wall, wall_s);
end

meta = struct();
meta.gfc_name = repo_rel(gfc_path);      % repository-relative, no local paths in outputs
meta.points_name = repo_rel(points_csv);
meta.tide_system = hdr.tide_system;
meta.product_type = hdr.product_type;
meta.GM = hdr.GM;
meta.R = hdr.R;
meta.nmax = nmax;
meta.fnALFs = fnALFs;
meta.recurrence = recurrence;
meta.set_name = set_name;
meta.load_s = load_s;
meta.load_cpu_s = load_cpu;
meta.threads = threads;
meta.repeat_count = repeat_count;
meta.max_diff = max_diff;
meta.graflab_dir = '<GrafLab folder: user-supplied, not redistributed>';
meta.matlab_release = version('-release');
meta.command = sprintf('run_graflab(''%s'', ''%s'', %d, ''%s'', ''%s'', ''%s'', %d, ''<GrafLab folder>'')', ...
    meta.gfc_name, meta.points_name, nmax, recurrence, set_name, repo_rel(out_dir), repeat_count);

write_result_files(out_dir, recurrence, set_name, pts, out_first, meta);
fprintf('run_graflab: wrote %s\n', fullfile(out_dir, sprintf('graflab_%s_%s.csv', recurrence, set_name)));
end

function p = to_absolute(p)
f = java.io.File(p);
p = char(f.getAbsolutePath());
end

function r = repo_rel(p)
%REPO_REL Path relative to benchmark/shs (for outputs). Other paths are reduced to their file name.
p = strrep(p, char(92), '/');   % Windows separators to '/'
k = strfind(p, 'benchmark/shs/');
if ~isempty(k)
    r = p(k(end) + numel('benchmark/shs/'):end);
else
    [~, nm, ext] = fileparts(p);
    r = [nm ext];
end
end
