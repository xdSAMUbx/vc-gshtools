# SHS benchmark: input data

Shared inputs of every engine of the spherical-harmonic synthesis benchmark. Definitions of the functionals,
units and file columns: `../results/csv.md`.

## Files

| File | In the repository | Content |
|---|---|---|
| `points_colorado.csv` | yes | 509 points in Colorado (`id,lat_deg,lon_deg,h_m`; geodetic latitude and longitude on GRS80, ellipsoidal height in m) |
| `points_latsweep.csv` | yes | 341 points on the meridian 10 deg E, latitude -85 to +85 deg every 0.5 deg, height 0 m |
| `XGM2019e_2190_tide_free.gfc.xz` | no (public repository: GitHub Release asset `data-v1`) | lossless xz of the tide-free model (65.9 MB) |
| `XGM2019e_2190_tide_free.gfc` | no (ignored by git) | uncompressed tide-free model: `xz -dk` of the file above |
| `SHA256SUMS.txt` | yes | SHA-256 of the original zero-tide model and of the derived files |

## Origin of the gravity model

* Model: **XGM2019e_2159**, from ICGEM (GFZ Potsdam, https://icgem.gfz-potsdam.de). Header: `modelname XGM2019e_2159`,
  `max_degree 2190`, `norm fully_normalized`, `earth_gravity_constant 3.986004415E+14`, `radius 6378136.3`,
  `errors formal`. The header says this is the spherical-harmonic (d/o 2190) version of the spheroidal model d/o 2159.
* Citation, taken from the header of the `.gfc` file: Zingerle, Philipp; Pail, Roland; Gruber, Thomas; Oikonomidou,
  Xanthi (2019): The experimental gravity field model XGM2019e. GFZ Data Services. http://doi.org/10.5880/ICGEM.2019.007
* Licence, from the header: **Creative Commons Attribution 4.0 International (CC BY 4.0)**. Any use or redistribution,
  including the derived file below, must credit the authors above and state that the model was modified
  (tide-free conversion of C20).
* The original zero-tide file `XGM2019e_2159.gfc` is not in this repository. Obtain it from ICGEM and check its
  SHA-256 against `SHA256SUMS.txt`.

## Tide system: derivation of the tide-free model

The original header declares `tide_system zero_tide`. The benchmark uses one common tide-free model, so every engine
reads the same file and no engine applies its own tide conversion.

The derivation is a streaming, byte-level copy of the original file in which **only** these changes are made:

1. Header line `tide_system                 zero_tide` becomes `tide_system                 tide_free` (same spacing).
2. One comment line is inserted directly after that header line, starting with `# Derived from XGM2019e_2159.gfc`.
   It is a header line that the gsh parser ignores (it reads only the keys it knows, up to `end_of_head`).
3. The C20 coefficient (n = 2, m = 0) changes:

   C20_TF = C20_ZT + 3.1108e-8 * 0.3 / sqrt(5) = C20_ZT + 4.173576e-9

| Quantity | Value |
|---|---|
| C20 zero-tide (original, `gfc 2 0` line) | `-4.841694947475625E-04` |
| C20 tide-free (derived) | `-4.841653211714038E-04` |
| Change | `+4.173576E-09` |

The formula and its constant are the ones used for the benchmark (see `../results/csv.md`). The magnitude
(4.1736e-9) is the usual zero-tide/tide-free C20 shift; the sign follows the benchmark definition. The conversion
was not re-derived here.

Verification: `diff` between the original and the derived file reports exactly three differences: the `tide_system`
line (replaced), the inserted comment line, and the C20 `gfc 2 0` line (value replaced, the other fields of that line
byte-identical). Every other line is copied byte for byte. The output format of the new C20 value is `%.15E`.

## Points

* **Colorado** (`points_colorado.csv`): 509 GPSBM benchmark points in Colorado, with the IGS08 geodetic (ellipsoidal)
  height of GRS80 (the `h_m` column). The NAVD88 orthometric height that the source also gives is not used, and no
  geoid height was added. Longitudes are in the -180..180 convention. Provenance: the GPSBM benchmark data of
  Colorado, extracted with a provenance script (`code/python/prepare_data.py`) that is not part of the public flow.
  No licence statement is made for the source data here; check the terms of the original provider before any reuse.
* **Latitude sweep** (`points_latsweep.csv`): made by `code/python/make_latitude_sweep.py`, a deterministic list
  (`id` 1..341, latitude -85.0 + 0.5 i, longitude 10.000000 deg, height 0.0 m). It samples the meridian from pole to
  pole, where the Legendre recurrences are stressed.

## Obtaining and checking the derived files

In the public repository the compressed model is not in the git tree. Download `XGM2019e_2190_tide_free.gfc.xz`
from the Releases page of the repository (tag `data-v1`) into this folder, then check it against the `.xz` line of
`SHA256SUMS.txt` before decompressing:

```
xz -dk benchmark/shs/data/XGM2019e_2190_tide_free.gfc.xz    # keeps the .xz
sha256sum -c benchmark/shs/data/SHA256SUMS.txt               # the line of the original zero-tide file fails unless it is copied here
```

The derived model is reproduced by `code/python/prepare_data.py` (provenance script, see its header: it needs the
private project package and is not runnable from the public tree). Compression of the derived file:
`xz -9 -T0 -k -f -c XGM2019e_2190_tide_free.gfc > XGM2019e_2190_tide_free.gfc.xz`.
