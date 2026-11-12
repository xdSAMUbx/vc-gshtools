function pts = read_points_csv(path)
%READ_POINTS_CSV Read a benchmark points file with header id,lat_deg,lon_deg,h_m.
%   Ids are kept as text so that any id style is preserved verbatim.

fid = fopen(path, 'r');
if fid < 0
    error('read_points_csv:open', 'Cannot open %s', path);
end
cleaner = onCleanup(@() fclose(fid));

header = strtrim(fgetl(fid));
if ~strcmp(header, 'id,lat_deg,lon_deg,h_m')
    error('read_points_csv:header', 'Unexpected header "%s" in %s', header, path);
end

ids = {};
lat = [];
lon = [];
h = [];
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    line = strtrim(line);
    if isempty(line)
        continue;
    end
    c = strsplit(line, ',');
    if numel(c) ~= 4
        error('read_points_csv:columns', 'Expected 4 columns, got %d in line: %s', numel(c), line);
    end
    ids{end+1, 1} = strtrim(c{1}); %#ok<AGROW>
    lat(end+1, 1) = str2double(c{2}); %#ok<AGROW>
    lon(end+1, 1) = str2double(c{3}); %#ok<AGROW>
    h(end+1, 1) = str2double(c{4}); %#ok<AGROW>
end

if isempty(lat)
    error('read_points_csv:empty', 'No points in %s', path);
end
if any(~isfinite([lat; lon; h])) || any(lat < -90 | lat > 90)
    error('read_points_csv:values', 'Non-finite values or latitude outside [-90, 90] in %s', path);
end

pts = struct('id', {ids}, 'lat', lat, 'lon', lon, 'h', h);
end
