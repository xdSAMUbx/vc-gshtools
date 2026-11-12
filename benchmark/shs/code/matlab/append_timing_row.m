function append_timing_row(path, engine, recurrence, set_name, nmax, n_points, threads, repeat_index, wall_s, cpu_s, notes)
%APPEND_TIMING_ROW Append one repetition to results/timings.csv (format of results/timings.csv).
%   Header: engine,recurrence,set,nmax,n_points,threads,repeat,wall_s,cpu_s,notes
%   The notes field is quoted, so commas are allowed in it.

if ~isfile(path)
    fid = fopen(path, 'w');
    fprintf(fid, 'engine,recurrence,set,nmax,n_points,threads,repeat,wall_s,cpu_s,notes\n');
else
    fid = fopen(path, 'a');
end
if fid < 0
    error('append_timing_row:open', 'Cannot open %s', path);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s,%s,%s,%d,%d,%d,%d,%.6f,%.6f,"%s"\n', engine, recurrence, set_name, ...
    nmax, n_points, threads, repeat_index, wall_s, cpu_s, strrep(notes, '"', ''''));
end
