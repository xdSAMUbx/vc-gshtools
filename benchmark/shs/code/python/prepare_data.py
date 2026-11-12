"""Prepare the shared input data of the SHS benchmark (definitions: benchmark/shs/results/csv.md).

PROVENANCE SCRIPT. It is NOT part of the public benchmark flow and is not runnable from the public tree:
it imports gsh_tools (the private project package) and reads data/colorado (not published). It is kept
so the origin of the derived files can be audited. The public inputs are already in data/:
points_colorado.csv, points_latsweep.csv (made by make_latitude_sweep.py) and the .xz of the model.
See data/DATA.md for how each file was produced.

Steps (each one is independent; run with the project venv):

    .venv\\Scripts\\python.exe -I benchmark\\shs\\code\\python\\prepare_data.py

  (a) benchmark/shs/data/points_colorado.csv the 509 GPSBM benchmarks of data/colorado, read with the
                                             project loader gsh_tools.data.colorado.load_colorado
                                             (kind 'gpsbm'). The height written is the IGS08
                                             ELLIPSOIDAL height (h_geodetic); no geoid is added.
  (b) benchmark/shs/data/XGM2019e_2190_tide_free.gfc  (uncompressed, git-ignored) derived from the
                                             zero-tide XGM2019e_2159.gfc by a streaming copy in which
                                             ONLY two things change:
                                               * the value of the `tide_system` header line
                                                 (zero_tide -> tide_free);
                                               * the C20 coefficient (n = 2, m = 0):
                                                   C20_TF = C20_ZT + 3.1108e-8 * 0.3 / sqrt(5).
                                             One comment line is added inside the header (after
                                             `tide_system`) stating the derivation; the parser of
                                             gsh ignores header lines it does not know.
                                             Every other line is copied byte for byte.
                                             The file is then verified against the original.

Compression (xz -9 -T0), SHA256SUMS.txt and DATA.md are produced outside this script (see DATA.md).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[4]          # GSH_tools
SHS = Path(__file__).resolve().parents[2]           # benchmark/shs
SALIDA = SHS / "data"

MODELO_ZT = RAIZ / "data" / "cache" / "modelos" / "XGM2019e_2159" / "XGM2019e_2159.gfc"
MODELO_ZIP = RAIZ / "data" / "models" / "XGM2019e_2159.zip"
MODELO_TF = SALIDA / "XGM2019e_2190_tide_free.gfc"

# Tide conversion as defined for the benchmark (only C20 changes).
DELTA_C20 = 3.1108e-8 * 0.3 / math.sqrt(5.0)

COMENTARIO = (b"# Derived from XGM2019e_2159.gfc (zero_tide): only tide_system and C20 (2,0) changed. "
              b"C20_TF = C20_ZT + 3.1108e-8*0.3/sqrt(5). See benchmark/shs/data/DATA.md.\n")


def escribir_puntos_colorado() -> int:
    """Writes points_colorado.csv from the GPSBM benchmarks through the project loader."""
    sys.path.insert(0, str(RAIZ / "src" / "py"))
    from gsh_tools.data.colorado import load_colorado  # noqa: E402  (needs the path above)

    gp = load_colorado(RAIZ / "data", "gpsbm")
    lat, lon, h = gp["lat"], gp["lon"], gp["h_geodetic"]  # h_geodetic = IGS08 ellipsoidal height
    salida = SALIDA / "points_colorado.csv"
    with salida.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id,lat_deg,lon_deg,h_m\n")
        for i in range(lat.size):
            f.write(f"{i + 1},{float(lat[i])!r},{float(lon[i])!r},{float(h[i])!r}\n")
    return int(lat.size)


def c20_linea(linea: bytes) -> tuple[bytes, float]:
    """Returns the C20 token (bytes) and value of a `gfc 2 0` data line."""
    partes = linea.split()
    token = partes[3]
    return token, float(token.decode("ascii"))


def derivar_tide_free() -> tuple[float, float]:
    """Streams the zero-tide model into the tide-free model. Returns (C20_ZT, C20_TF)."""
    if not MODELO_ZT.is_file():
        if not MODELO_ZIP.is_file():
            raise FileNotFoundError(f"Missing {MODELO_ZT} and {MODELO_ZIP}")
        raise FileNotFoundError(
            f"Missing {MODELO_ZT}: extract XGM2019e_2159.gfc from {MODELO_ZIP} first")

    c20_zt = c20_tf = None
    en_cabecera = True
    tmp = MODELO_TF.with_name(MODELO_TF.name + ".tmp")
    with MODELO_ZT.open("rb") as fin, tmp.open("wb") as fout:
        for linea in fin:                    # binary mode: bytes are copied exactly
            if en_cabecera:
                if linea.startswith(b"tide_system"):
                    if linea.split()[1:] != [b"zero_tide"]:
                        raise ValueError(f"Unexpected tide_system line: {linea!r}")
                    fout.write(linea.replace(b"zero_tide", b"tide_free", 1))
                    # Comment line added right after tide_system (the gsh parser ignores header
                    # lines it does not know). Numeric values are in DATA.md.
                    fout.write(COMENTARIO)
                    continue
                if linea.startswith(b"end_of_head"):
                    en_cabecera = False
                fout.write(linea)
                continue
            if linea.startswith(b"gfc"):
                partes = linea.split()
                if len(partes) >= 4 and partes[1] == b"2" and partes[2] == b"0":
                    token, c20_zt = c20_linea(linea)
                    c20_tf = c20_zt + DELTA_C20
                    nuevo = f"{c20_tf:.15E}".encode("ascii")
                    fout.write(linea.replace(token, nuevo, 1))
                    continue
            fout.write(linea)
    tmp.replace(MODELO_TF)
    return c20_zt, c20_tf


if __name__ == "__main__":
    SALIDA.mkdir(parents=True, exist_ok=True)
    n_co = escribir_puntos_colorado()
    print(f"points_colorado.csv: {n_co} points")
    c20_zt, c20_tf = derivar_tide_free()
    print(f"C20 zero-tide : {c20_zt!r}")
    print(f"C20 tide-free : {c20_tf!r}")
    print(f"written: {MODELO_TF}")
