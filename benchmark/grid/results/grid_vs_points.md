# Grid path vs point path (gsh)

Criterion: max|grid - points| / max|points| <= 1e-12 for each field (fixed before the run).

Model: `EGM2008.gfc`, truncated to the degree of each case. Spherical latitudes, constant radius R per row.
Grid side: `synth_grid` on the full grid (V: grad 0 of the model; T and dg: grad 1 of the disturbing model), assembled
with the formulas of `functionals.cpp`. Point side: `grid_to_points` on the selected rows, then `compute_functionals`
with outputs V, T, dg. Selected rows: auto (global) or the 5 regional rows; exact-pole rows excluded.

| case | nmax | rows | points | field | max abs diff | max abs value | relative diff | row of max diff (lat deg) | verdict |
|---|---|---|---|---|---|---|---|---|---|
| gl | 2190 | 11 | 48202 | V | 7.302e-07 | 6.253e+07 | 1.168e-14 | row 1095 (0.000000) | PASS |
| gl | 2190 | 11 | 48202 | T | 9.436e-12 | 1.011e+03 | 9.335e-15 | row 1095 (0.000000) | PASS |
| gl | 2190 | 11 | 48202 | dg | 9.634e-16 | 3.512e-03 | 2.743e-13 | row 1095 (0.000000) | PASS |
| regional | 2190 | 5 | 1005 | V | 2.086e-07 | 6.250e+07 | 3.338e-15 | row 200 (35.000000) | PASS |
| regional | 2190 | 5 | 1005 | T | 3.382e-12 | 2.716e+02 | 1.245e-14 | row 150 (37.500000) | PASS |
| regional | 2190 | 5 | 1005 | dg | 2.662e-17 | 1.222e-03 | 2.178e-14 | row 150 (37.500000) | PASS |
| dh1 | 720 | 10 | 14420 | V | 4.545e-07 | 6.253e+07 | 7.268e-15 | row 721 (0.000000) | PASS |
| dh1 | 720 | 10 | 14420 | T | 7.276e-12 | 1.011e+03 | 7.199e-15 | row 961 (-29.958391) | PASS |
| dh1 | 720 | 10 | 14420 | dg | 2.069e-16 | 2.449e-03 | 8.448e-14 | row 721 (0.000000) | PASS |

Wall time of each side in this single accuracy run (not a timing benchmark):

| case | grid side (s) | point side (s) |
|---|---|---|
| gl | 1.5 | 5344.2 |
| regional | 0.5 | 0.5 |
| dh1 | 0.2 | 107.3 |
