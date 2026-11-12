function hdr = read_gfc_header(path)
%READ_GFC_HEADER Read the ICGEM header of a .gfc file (up to end_of_head).
%   Returns GM, R, max_degree, norm, tide_system, product_type. Errors if the
%   coefficients are not fully normalized (GrafLab requires that) or if the first
%   data line does not have the 7 columns GrafLab reads (gfc deg ord C S sigC sigS).

fid = fopen(path, 'r');
if fid < 0
    error('read_gfc_header:open', 'Cannot open %s', path);
end
cleaner = onCleanup(@() fclose(fid));

hdr = struct('GM', NaN, 'R', NaN, 'max_degree', NaN, 'norm', '', ...
    'tide_system', '', 'product_type', '', 'n_columns_first_data', NaN);

found_end = false;
while true
    s = fgetl(fid);
    if ~ischar(s)
        break;
    end
    [key, rest] = strtok(strtrim(s));
    val = strtrim(rest);
    arg = strtok(val);
    switch lower(key)
        case 'earth_gravity_constant'
            hdr.GM = str2double(arg);
        case 'radius'
            hdr.R = str2double(arg);
        case 'max_degree'
            hdr.max_degree = str2double(arg);
        case 'norm'
            hdr.norm = val;
        case 'tide_system'
            hdr.tide_system = val;
        case 'product_type'
            hdr.product_type = val;
        case 'end_of_head'
            found_end = true;
    end
    if found_end
        break;
    end
end
if ~found_end
    error('read_gfc_header:noEndOfHead', 'No end_of_head line in %s', path);
end

if ~strcmpi(hdr.norm, 'fully_normalized')
    error('read_gfc_header:norm', 'GrafLab needs fully normalized coefficients; header says "%s".', hdr.norm);
end

s = fgetl(fid);
if ischar(s)
    hdr.n_columns_first_data = numel(strsplit(strtrim(s)));
end
if hdr.n_columns_first_data < 7
    warning('read_gfc_header:columns', ...
        'First data line of %s has %d columns; GrafLab expects 7 (gfc deg ord C S sigC sigS).', ...
        path, hdr.n_columns_first_data);
end
end
