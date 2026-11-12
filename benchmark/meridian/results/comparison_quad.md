# Meridian benchmark: engines against the binary128 reference

Reference: `ref_quad/ref_quad.cpp` (IEEE binary128, ~1e-34), EGM2008 n = 2190, GRS80, lon 0, h = 0, 181 points. Coefficients, GM, R and GRS80 constants are rounded to double first (as the engines read them).

Cell = max|engine - reference| / max|reference|, verdict against 1e-12, latitude of the maximum, and for dg the maximum in ULP of g (1 ULP = 1.78e-15 m/s^2). dg excludes the two poles (gsh rejects 1/cos(phi) there); V_xx excludes the poles (no d2P closed form there in the reference).

### Position: exact decimal latitude

| functional | gsh | GrafLab era | GrafLab msfc | GrafLab sfc |
|---|---|---|---|---|
| dg | 3.98e-13 (yes) @+75, 0.4 ULP(g) | 5.56e-11 (no) @+35, 52.0 ULP(g) | 2.67e-11 (no) @+42, 25.0 ULP(g) | 4.00e-03 (no) @+62, 3741760830.0 ULP(g) |
| V_xz | 5.89e-13 (yes) @-89 | 4.68e-12 (no) @+90 | 1.42e-12 (no) @+90 | 1.05e-02 (no) @+69 |
| V_xx (control) | 1.64e-14 (yes) @+89 | 7.45e-14 (yes) @+43 | 7.40e-14 (yes) @+43 | 2.80e-04 (no) @+62 |
| V_zz (control) | 7.53e-15 (yes) @+15 | 9.67e-14 (yes) @+90 | 5.25e-14 (yes) @+43 | 7.07e-04 (no) @+62 |

### Position: latitude as the double radian used by the engines

| functional | gsh | GrafLab era | GrafLab msfc | GrafLab sfc |
|---|---|---|---|---|
| dg | 3.99e-13 (yes) @+75, 0.4 ULP(g) | 5.56e-11 (no) @+35, 52.0 ULP(g) | 2.67e-11 (no) @+42, 25.0 ULP(g) | 4.00e-03 (no) @+62, 3741760830.0 ULP(g) |
| V_xz | 5.91e-13 (yes) @-89 | 4.68e-12 (no) @+90 | 1.42e-12 (no) @+90 | 1.05e-02 (no) @+69 |
| V_xx (control) | 1.71e-14 (yes) @+89 | 7.53e-14 (yes) @+43 | 7.48e-14 (yes) @+43 | 2.80e-04 (no) @+62 |
| V_zz (control) | 7.67e-15 (yes) @+15 | 9.67e-14 (yes) @+90 | 5.27e-14 (yes) @+43 | 7.07e-04 (no) @+62 |
