# H1 and H2: parallel preparation and fused synthesis

Colorado, 509 points, V + T + dg, XGM2019e nmax 2190, `gsh_shs_bench --recurrence fast`, 8 repetitions, machine idle
(AMD Ryzen 7 7445HS, 6 cores / 12 threads, AVX-512). The reference column is the round measured before these changes
(`output` of that round, same machine and points).

- **H1:** serial setup steps parallelised (recurrence tables moved instead of copied, packing loops in parallel).
- **H2:** when the model and the disturbing set are both needed, one synthesis of the disturbing set at the highest
  gradient order plus a zonal series (degree at most 20) gives the model sums. `--fuse 0` switches it off.

| Threads | Reference median (s) | H1 only (`--fuse 0`) | H1 + H2 (`--fuse 1`) | Gain of H2 |
|---:|---:|---:|---:|---:|
| 1 | 1.432 | 1.448 | 1.066 | -26 % |
| 6 | 0.429 | 0.398 | 0.280 | -30 % |
| 12 | 0.316 | 0.306 | 0.232 | -24 % |

H1 alone is within 0-7 % (inside the noise at 1 thread). H2 is the real gain.

Accuracy, fused vs two-pass on the same points: T and dg identical bit for bit, V differs by 9.5e-15 (relative to
max|V|); criterion 1e-12. The 45-functional case (37.5 % less time) is in [`../../../aero5540/README.md`](../../../aero5540/README.md).

Files: `timings_fuse1.csv`, `timings_fuse0.csv` (one row per repetition), `out_fuse{1,0}_t{1,6,12}.csv` (values),
`log_*.txt`. Reproduce: `gsh_shs_bench --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc --points
benchmark/shs/data/points_colorado.csv --recurrence fast --nmax 2190 --set colorado --threads K --repeat 8 --fuse 0|1
--out <csv> --timings <csv>`.
