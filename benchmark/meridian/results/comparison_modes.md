# Meridian benchmark: all GrafLab Legendre modes

EGM2008, degree and order 2190, GRS80, longitude 0 deg, 181 points, h = 0 m. gsh is the same single run (`compute_functionals`) in every column. GrafLab 2.1.4 is run with `fnALFs` = 3 (era), 2 (msfc) and 1 (sfc).

Metric: max |d| / max |reference|, criterion <= 1e-12; the verdict is in brackets. The last columns compare two GrafLab runs with each other (spread of GrafLab itself).

| functional | gsh vs GrafLab era | gsh vs GrafLab msfc | gsh vs GrafLab sfc | GrafLab msfc vs era | GrafLab sfc vs era |
|---|---|---|---|---|---|
| T | 8.99e-14 (yes) | 9.26e-14 (yes) | 3.85e-05 (no) | 9.66e-15 (yes) | 3.85e-05 (no) |
| zeta_ell | 1.06e-13 (yes) | 1.16e-13 (yes) | 3.85e-05 (no) | 9.67e-15 (yes) | 3.85e-05 (no) |
| Dg_sa (gravity anomaly, sa) | 7.74e-13 (yes) | 5.69e-13 (yes) | 4.17e-03 (no) | 5.12e-13 (yes) | 4.17e-03 (no) |
| dg (gravity disturbance) | 5.52e-11 (no) | 2.71e-11 (no) | 4.00e-03 (no) | 5.02e-11 (no) | 4.00e-03 (no) |
| V_xx | 1.98e-13 (yes) | 7.67e-14 (yes) | 2.80e-04 (no) | 1.36e-13 (yes) | 2.80e-04 (no) |
| V_yy | 3.90e-13 (yes) | 1.34e-13 (yes) | 1.14e-03 (no) | 2.57e-13 (yes) | 1.14e-03 (no) |
| V_zz | 9.72e-14 (yes) | 5.30e-14 (yes) | 7.07e-04 (no) | 5.86e-14 (yes) | 7.07e-04 (no) |
| V_xz | 4.67e-12 (no) | 1.41e-12 (no) | 1.05e-02 (no) | 3.26e-12 (no) | 1.05e-02 (no) |

dg excludes the two exact poles for the gsh columns (gsh rejects `1/cos(phi)` there); in the GrafLab-against-GrafLab columns all 181 points are used.
