function run_graflab_meridian(gfc_path, out_csv, graflab_dir, fnALFs)
%RUN_GRAFLAB_MERIDIAN GrafLab 2.1.4 reference values for the meridional benchmark (see ../README.md).
%
%   run_graflab_meridian(gfc_path, out_csv, graflab_dir, fnALFs)
%
%   gfc_path     Tide-free ICGEM .gfc, full degree (EGM2008 in this benchmark).
%   out_csv      Result CSV, written with %.17g (full double precision). Columns and units: see README.md.
%   graflab_dir  Folder that contains GrafLab.m. GrafLab is NOT shipped in this repository. Optional: when
%                omitted, the environment variable GRAFLAB_DIR is used.
%   fnALFs       Legendre recurrence of GrafLab: 1 standard forward column, 2 modified forward column,
%                3 extended-range arithmetic (default 3: the reference recurrence of this benchmark).
%
%   Points: 181 points, geodetic latitude -90 .. +90 deg in steps of 1 deg, longitude 0 deg, height 0 m,
%   GRS80 (GrafLab ellipsoid = 1). Point-wise mode (point_type = 2), no grid, no interpolation.
%
%   GrafLab accepts at most 4 functionals per call, so two calls are made:
%     call 1: Functional = [5; 22; 19; 20]   T, zeta_Ell, Dg_sa (gravity anomaly, spherical approx.), dg
%     call 2: Functional = [14; 15]          V_xx,V_yy,V_zz (code 14) and V_xy,V_xz,V_yz (code 15)
%   GrafLab code 14 exports the columns [Vxx, -Vyy, Vzz] (see GrafLab.m, the "Pg=[Pg -Vyy(:)]" line):
%   the sign of the Vyy column is written as GrafLab exports it and is named Vyy_as_exported_E in the CSV.
%
%   Example (PowerShell, from the repository root):
%     matlab -batch "cd('benchmark\meridian\matlab'); run_graflab_meridian('data\cache\modelos\EGM2008\EGM2008.gfc', 'benchmark\meridian\results\graflab_era_egm2008_2190.csv', 'C:\...\GrafLab', 3)"

if nargin < 4 || isempty(fnALFs)
    fnALFs = 3;
end
if nargin < 3 || isempty(graflab_dir)
    graflab_dir = getenv('GRAFLAB_DIR');
end
if isempty(graflab_dir) || ~isfile(fullfile(graflab_dir, 'GrafLab.m'))
    error('run_graflab_meridian:noGrafLab', ...
        'GrafLab.m not found in "%s". Pass its folder as the 3rd argument or set GRAFLAB_DIR.', graflab_dir);
end
if ~isfile(gfc_path)
    error('run_graflab_meridian:noGfc', 'Geopotential model not found: %s', gfc_path);
end
if ~any(fnALFs == [1 2 3])
    error('run_graflab_meridian:badRecurrence', 'fnALFs must be 1, 2 or 3.');
end
addpath(graflab_dir);

[GM, R, nmax] = read_icgem_header(gfc_path);

lat = (-90:1:90)';
lon = zeros(size(lat));
h = zeros(size(lat));
n_points = numel(lat);

t0 = tic;
c1 = graflab_call(gfc_path, GM, R, nmax, fnALFs, lat, lon, h, [5; 22; 19; 20]);
t_call1 = toc(t0);
fprintf('run_graflab_meridian: call 1 (T, zeta_Ell, Dg_sa, dg) done in %.1f s\n', t_call1);

t0 = tic;
c2 = graflab_call(gfc_path, GM, R, nmax, fnALFs, lat, lon, h, [14; 15]);
t_call2 = toc(t0);
fprintf('run_graflab_meridian: call 2 (V tensor, codes 14 and 15) done in %.1f s\n', t_call2);

if size(c1, 1) ~= n_points || size(c1, 2) ~= 7 || size(c2, 1) ~= n_points || size(c2, 2) ~= 9
    error('run_graflab_meridian:badOutput', 'Unexpected GrafLab output sizes [%d %d] and [%d %d].', ...
        size(c1, 1), size(c1, 2), size(c2, 1), size(c2, 2));
end

% Column layout of GrafLab output: 1..3 = lat, lon, h; then one column per functional entry (see header).
T_m2s2 = c1(:, 4);
zeta_ell_m = c1(:, 5);
Dg_sa_mgal = c1(:, 6);
dg_mgal = c1(:, 7);
Vxx_E = c2(:, 4);
Vyy_as_exported_E = c2(:, 5);
Vzz_E = c2(:, 6);
Vxy_E = c2(:, 7);
Vxz_E = c2(:, 8);
Vyz_E = c2(:, 9);

% Only the gravity disturbance at the exact poles may be non-finite (GrafLab 1/cos(phi) terms); report, do not fix.
bad = ~isfinite([T_m2s2 zeta_ell_m Dg_sa_mgal Vxx_E Vyy_as_exported_E Vzz_E Vxy_E Vxz_E Vyz_E]);
if any(bad(:))
    fprintf('run_graflab_meridian: %d non-finite GrafLab values written as NaN\n', nnz(bad));
end
dg_nonfinite = ~isfinite(dg_mgal);

out_dir = fileparts(out_csv);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
fid = fopen(out_csv, 'w');
if fid < 0
    error('run_graflab_meridian:open', 'Cannot write %s', out_csv);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, 'lat_deg,lon_deg,h_m,T_m2s2,zeta_ell_m,Dg_sa_mgal,dg_mgal,Vxx_E,Vyy_as_exported_E,Vzz_E,Vxy_E,Vxz_E,Vyz_E\n');
for i = 1:n_points
    fprintf(fid, '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
        lat(i), lon(i), h(i), T_m2s2(i), zeta_ell_m(i), Dg_sa_mgal(i), dg_mgal(i), ...
        Vxx_E(i), Vyy_as_exported_E(i), Vzz_E(i), Vxy_E(i), Vxz_E(i), Vyz_E(i));
end
fprintf('run_graflab_meridian: wrote %s (%d points, fnALFs=%d, nmax=%d, dg non-finite at %d points)\n', ...
    out_csv, n_points, fnALFs, nmax, nnz(dg_nonfinite));
end

function out = graflab_call(gfc_path, GM, R, nmax, fnALFs, lat, lon, h, functional)
%GRAFLAB_CALL One headless GrafLab 2.1.4 call (argument list as in code/matlab/graflab_evaluate.m of the SHS benchmark).
nmin = 0;                 % required by dg (Functional 20); the tensors of code 14 and 15 accept nmin = 0
ellipsoid = 1;            % GRS80 inside GrafLab
coord = 0;                % ellipsoidal input (geodetic latitude, ellipsoidal height)
point_type = 2;           % point-wise mode
lat_min = []; lat_step = []; lat_max = [];
lon_min = []; lon_step = []; lon_max = [];
hgrid = [];
Input_data_path = [];
Output_path = fullfile(tempdir, 'graflab_unused_prefix');   % not written: exports are disabled
Functional_or_Commission = 0;
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

out = GrafLab('OK', GM, R, nmin, nmax, ellipsoid, gfc_path, coord, point_type, ...
    lat_min, lat_step, lat_max, lon_min, lon_step, lon_max, hgrid, Input_data_path, ...
    lat, lon, h, Output_path, Functional_or_Commission, functional, fnALFs, DTM_path, ...
    Export_data_txt, Export_report, Export_data_mat, ...
    Display_data, Graphic_format, Colormap, Number_of_colors, DPI, Status_bar);
end

function [GM, R, nmax] = read_icgem_header(path)
%READ_ICGEM_HEADER earth_gravity_constant, radius and max_degree of an ICGEM .gfc header.
fid = fopen(path, 'r');
if fid < 0
    error('run_graflab_meridian:open', 'Cannot open %s', path);
end
cleaner = onCleanup(@() fclose(fid));
GM = NaN; R = NaN; nmax = NaN; norm_ok = false;
while true
    s = fgetl(fid);
    if ~ischar(s)
        error('run_graflab_meridian:noEndOfHead', 'No end_of_head line in %s', path);
    end
    [key, rest] = strtok(strtrim(s));
    arg = strtok(strtrim(rest));
    switch lower(key)
        case 'earth_gravity_constant'
            GM = str2double(arg);
        case 'radius'
            R = str2double(arg);
        case 'max_degree'
            nmax = str2double(arg);
        case 'norm'
            norm_ok = strcmpi(arg, 'fully_normalized');
        case 'end_of_head'
            break;
    end
end
if ~norm_ok
    error('run_graflab_meridian:norm', 'GrafLab needs fully normalized coefficients (header norm line).');
end
end
