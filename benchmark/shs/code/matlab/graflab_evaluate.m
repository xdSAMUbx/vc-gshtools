function out = graflab_evaluate(p)
%GRAFLAB_EVALUATE One headless GrafLab 2.1.4 call: V, T and dg at the given points, in one pass.
%
%   p fields: gfc_path, GM, R, nmax, fnALFs, lat (deg), lon (deg), h (m).
%
%   Mapping to GrafLab('OK', ...) argument list (see GrafLab.m, function header):
%     GM, R              values of the .gfc header (GrafLab overwrites them with the file's values)
%     nmin = 0           required: GrafLab refuses dg (Functional 20) when nmin > 0
%     nmax               degree of the expansion
%     ellipsoid = 1      GRS80 (GM 3.986005e14, a 6378137, e^2 0.006694380022903416, omega 7292115e-11)
%     GGM_path           the tide-free .gfc (GrafLab reads GM and R from its header)
%     coord = 0          ellipsoidal input coordinates (geodetic latitude, ellipsoidal height)
%     point_type = 2     point-wise mode: lat/lon/h2 are arrays, no grid, no interpolation
%     Functional = [11; 5; 20]
%                        11 = gravitational potential V (m^2/s^2). GrafLab takes the coefficients of the file
%                             directly for this functional (no delta C00 term, see GrafLab.m, grav==1 branch),
%                             so degree 0 is C00 of the file.
%                         5 = disturbing potential T (m^2/s^2)
%                        20 = gravity disturbance, full (not spherical approximation), mGal
%                        The three functionals are computed in one call. Reading GrafLab.m, the coefficient sets
%                        they need (C for V and dg, delta C for T, normal-field C for dg) are each built from the
%                        same file as in the former separate calls; this has not been checked by running GrafLab.
%     fnALFs             1 standard forward column, 2 modified forward column, 3 extended-range arithmetic
%     Export_* = 0, Status_bar = 0
%                        nargout == 1 makes GrafLab return the matrix in memory with full
%                        double precision: rows = points, columns = lat, lon, h, then one column per
%                        functional in the order of the Functional vector, i.e. V, T, dg. GrafLab does not sort
%                        the Functional vector: the functional loop runs over Functional(:) in the given order and
%                        appends each result as a column (GrafLab.m, the "P=[P Pg]" assembly). No file is written.
%
%   The GrafLab 'OK' branch is entered only when nargin > 1 (GUI is used when nargin == 1).

GM = p.GM;
R = p.R;
nmin = 0;
nmax = p.nmax;
ellipsoid = 1;
GGM_path = p.gfc_path;
coord = 0;
point_type = 2;
lat_min = []; lat_step = []; lat_max = [];
lon_min = []; lon_step = []; lon_max = [];
h = [];
Input_data_path = [];
lat = p.lat(:);
lon = p.lon(:);
h2 = p.h(:);
Output_path = fullfile(tempdir, 'graflab_unused_prefix'); % not written: exports are disabled
Functional_or_Commission = 0;
Functional = [11; 5; 20];     % V, T, dg (output columns in this order)
fnALFs = p.fnALFs;
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

out = GrafLab('OK', GM, R, nmin, nmax, ellipsoid, GGM_path, coord, point_type, ...
    lat_min, lat_step, lat_max, lon_min, lon_step, lon_max, h, Input_data_path, ...
    lat, lon, h2, Output_path, Functional_or_Commission, Functional, fnALFs, DTM_path, ...
    Export_data_txt, Export_report, Export_data_mat, ...
    Display_data, Graphic_format, Colormap, Number_of_colors, DPI, Status_bar);
end
