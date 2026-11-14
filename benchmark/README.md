# Benchmarks

| folder | what it checks | external reference |
|---|---|---|
| [`shs/`](shs/README.md) | synthesis accuracy and speed, XGM2019e n=2190, 4 engines ([Summary](shs/results/Summary.md)) | GrafLab, CHarm, SHTOOLS |
| [`meridian/`](meridian/README.md) | T, zeta_ell, Dg, dg and tensor along a meridian, EGM2008 n=2190 | GrafLab 2.1.4 |
| [`alf/`](alf/README.md) | Legendre sum-rule identities; fast kernel vs oracle | analytic (none external) |
| [`grid/`](grid/README.md) | grid synthesis (Gauss-Legendre, Driscoll-Healy, regional) vs the same values computed point by point; timings | the library's own point synthesis |
| [`python_overhead/`](python_overhead/README.md) | cost of calling the library from Python (`gshtools`) vs the C++ call, 30 cases, bitwise equality | C++ reference program |
| [`aero5540/`](aero5540/README.md) | official run at degree 5540: 45 functionals on a 283 716-point GRAV-D flight; per-functional time, fused-pass effect, accuracy | sequential oracle |

Speed round after the first shs results: preparation in parallel (H1) and one fused synthesis for the model and disturbing
sets (H2), measured in [`shs/results/h1h2/`](shs/results/h1h2/README.md) and [`aero5540/results/h2_fuse/`](aero5540/README.md).

Rules used by every benchmark here: timings with the machine idle and 8 repetitions (the one-per-functional full-flight
run of `aero5540` is the stated exception); accuracy is computed once, as max|d| / max|ref| with a 1e-12 criterion.
Raw results live in each folder's `results/`.

Error-model validations (commission/omission, null ellipsoid, Laplace trace) are unit tests:
`tests/cpp/test_errors.cpp`, `tests/cpp/test_error_report.cpp`, `tests/cpp/test_sanity_identities.cpp`; see
[`docs/modules/errors.md`](../docs/modules/errors.md).
