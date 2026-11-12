# Meridian benchmark: gsh against GrafLab 2.1.4 (era)

EGM2008 (`data/cache/modelos/EGM2008/EGM2008.gfc`), degree and order 2190, GRS80, one meridian at longitude 0 deg, 181 points of geodetic latitude -90..+90 deg, h = 0 m.

Success criterion: max |d| / max |GrafLab| <= 1e-12 (rounding level, stated by the user). Each row is reported as it is; no tolerance was changed.

| functional | unit | n | max abs(d) | max abs(GrafLab) | max abs(d)/max abs(GrafLab) | max rel(pt) | rms(d) | lat of max abs(d) [deg] | meets 1e-12 |
|---|---|---|---|---|---|---|---|---|---|
| T | m^2 s^-2 | 181 | 4.630e-11 | 5.151e+02 | 8.989e-14 | 1.749e-12 | 2.062e-11 | 85.0 | yes |
| zeta_ell | m | 181 | 5.588e-12 | 5.254e+01 | 1.064e-13 | 1.698e-12 | 2.185e-12 | -87.0 | yes |
| Dg_sa (gravity anomaly, sa) | m s^-2 | 181 | 1.232e-15 | 1.592e-03 | 7.739e-13 | 5.102e-11 | 1.727e-16 | 90.0 | yes |
| dg (gravity disturbance) | m s^-2 | 179 | 9.176e-14 | 1.663e-03 | 5.518e-11 | 2.295e-08 | 2.089e-14 | 35.0 | no |
| V_xx | s^-2 | 181 | 3.138e-19 | 1.588e-06 | 1.976e-13 | 2.028e-13 | 2.987e-20 | 90.0 | yes |
| V_yy | s^-2 | 181 | 6.118e-19 | 1.567e-06 | 3.904e-13 | 3.998e-13 | 4.733e-20 | 90.0 | yes |
| V_zz | s^-2 | 181 | 3.062e-19 | 3.149e-06 | 9.724e-14 | 9.948e-14 | 3.613e-20 | 90.0 | yes |
| V_xz | s^-2 | 181 | 2.786e-19 | 5.971e-08 | 4.666e-12 | 1.463e-10 | 2.812e-20 | 90.0 | no |

## Conventions and mapping

- gsh: `compute_functionals` (`Functional::T, zeta_ell, Dg_sa, dg, V_xx, V_yy, V_zz, V_xz`), SI units.
- GrafLab: `Functional = [5; 22; 19; 20]` and `[14; 15]`, `fnALFs = 3` (extended-range arithmetic). Units as GrafLab returns them: mGal and Eotvos, converted here to SI.
- dg at the exact poles: gsh does not evaluate it (1/cos(phi)); those two rows are excluded (n = 179).
- zeta: GrafLab functional 23 needs a DTM file, which this benchmark does not have. The comparison uses zeta_Ell (GrafLab 22) against gsh `zeta_ell` = T/gamma. The full zeta is not compared.
- Dg (gravity anomaly): GrafLab has only the spherical-approximation anomaly (functional 19); the gsh counterpart is `Dg_sa`. The exact gravity anomaly is not a GrafLab output.
- V_yy: GrafLab code 14 writes `[Vxx, -Vyy, Vzz]` (GrafLab.m, `Pg=[Pg -Vyy(:)]`). The comparison is made with the column as exported; the sign-flipped comparison is reported below as a diagnostic.
- max rel(pt) is the largest |d|/|GrafLab| at points with a non-zero GrafLab value. It is not meaningful where the functional passes through zero (e.g. near the equator), so the normalised column is the criterion.

## Diagnostics

- Laplace residual max|Vxx+Vyy+Vzz| / max|Vxx|: gsh 2.000e-14, GrafLab as exported 1.653e-14, GrafLab with the Vyy sign flipped 1.973e+00.
- gsh V_yy against -GrafLab Vyy (sign flipped): max abs(d) 3.134e-06, normalised 2.000e+00.

Raw numbers: `comparison.csv`. Inputs: `gsh_meridian.csv`, `graflab_era_meridian.csv`.
