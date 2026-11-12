# Benchmarks

| folder | what it checks | external reference |
|---|---|---|
| [`shs/`](shs/README.md) | synthesis accuracy and speed, XGM2019e n=2190, 4 engines ([Summary](shs/results/Summary.md)) | GrafLab, CHarm, SHTOOLS |
| [`meridian/`](meridian/README.md) | T, zeta_ell, Dg, dg and tensor along a meridian, EGM2008 n=2190 | GrafLab 2.1.4 |
| [`alf/`](alf/README.md) | Legendre sum-rule identities; fast kernel vs oracle | analytic (none external) |

Error-model validations (commission/omission, null ellipsoid, Laplace trace) are unit tests:
`tests/test_errors.cpp`, `tests/test_error_report.cpp`, `tests/test_sanity_identities.cpp`; see
[`docs/modules/errors.md`](../docs/modules/errors.md).
