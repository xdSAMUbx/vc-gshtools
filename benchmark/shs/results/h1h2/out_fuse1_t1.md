# out_fuse1_t1.csv

Engine: **gsh** (C++ driver `gsh_shs_bench`), recurrence `fast`, point set `colorado`, functionals V, T and dg computed in one pass.

## Columns

| column | unit | meaning |
|---|---|---|
| id | - | point identifier from the input CSV |
| lat_deg | deg | geodetic latitude on GRS80 (input, echoed) |
| lon_deg | deg | longitude (input, echoed) |
| h_m | m | ellipsoidal height on GRS80 (input, echoed) |
| V_m2s2 | m^2/s^2 | gravitational potential V = GM/r sum_n (R/r)^n sum_m (C cos m lambda + S sin m lambda) P_nm(sin phi') over n = 0..nmax, with the coefficients of the model as read from the file (C00 included); geocentric r and phi' from the GRS80 geodetic input; no normal potential subtracted, so V is not T |
| T_m2s2 | m^2/s^2 | disturbing potential T = V - U at the point (V: model, U: GRS80 normal potential, closed form; T = GM/r sum_n (R/r)^n sum_m (dC cos m lambda + dS sin m lambda) P_nm, with the disturbing coefficients dC, dS of the project) |
| dg_mgal | mGal | gravity disturbance dg = \|g\| - gamma_SH at the point (g: gravity vector of the model in the LNOF frame, including the centrifugal term; gamma_SH: normal gravity of the GRS80 zonal series of the project) |

## Command and parameters

```
gsh_shs_bench --gfc benchmark\shs\data\XGM2019e_2190_tide_free.gfc --points benchmark\shs\data\points_colorado.csv --recurrence fast --nmax 2190 --out <repo>/benchmark/shs/results/h1h2/out_fuse1_t1.csv --timings <repo>/benchmark/shs/results/h1h2/timings_fuse1.csv --set colorado --threads 1 --repeat 8 --fuse 1
```

- recurrence: `fast`
- nmax: `2190`, nmin: `0`
- threads: `1` (0 on the command line means the OpenMP default)
- repeat: `8` (the result is written from the last repetition)
- points: `509`
- method: `gsh::compute_functionals` with the outputs `{V, T, dg}` in one call: the model is synthesised once with gradient order 1 (its S0 gives V, its gradient sums give dg) and the disturbing model once with order 0 (its S0 gives T). SIMD kernels, OpenMP over points.

## Model

- file: `benchmark\shs\data\XGM2019e_2190_tide_free.gfc`
- modelname: `XGM2019e_2159`
- GM [m^3/s^2]: `3.986e+14`, R [m]: `6.37814e+06`
- tide_system (header): `tide_free` (no tide conversion applied by the driver)
- load_s: `0.519628`

## Notes

- V is the full gravitational potential of the model. It is NOT the disturbing potential T = V - U, and no normal potential is subtracted. Compare it with GrafLab Functional 11 with nmin = 0.
- Degree 0 is included in V: the term GM/r * C00 uses C00 from the file (1.0 in the tide-free model).
- Geodetic to geocentric conversion: `Points::from_geodetic` (GRS80), the same for all recurrences.
