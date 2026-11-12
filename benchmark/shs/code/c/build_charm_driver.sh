#!/usr/bin/env bash
# Builds the CHarm shared library and the charm_shs_bench driver on native Windows (MSYS2 UCRT64).
#
# Run from the MSYS2 UCRT64 shell (C:\msys64\ucrt64.exe), with PATH starting in /ucrt64/bin:
#     bash benchmark/shs/code/c/build_charm_driver.sh
#
# Prerequisites: CHarm v0.4.11 configured and built in build/external/charm (see README.md in this folder),
# i.e. build/external/charm/.libs/libcharm.a exists. Nothing is installed system-wide and nothing is written
# outside build/external/charm-driver.
#
# Why a DLL: CHarm's public headers declare its functions __declspec(dllimport) on Windows, so the driver
# must link to an import library. libtool in CHarm's autotools build produces only the static archive here,
# so this script links the archive into libcharm.dll (all symbols exported) and creates libcharm.dll.a.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
CHARM="$REPO/build/external/charm"
OUT="$REPO/build/external/charm-driver"
FLAGS="-O3 -march=native -fopenmp"

if [ ! -f "$CHARM/.libs/libcharm.a" ]; then
  echo "error: $CHARM/.libs/libcharm.a not found; build CHarm first (see README.md)" >&2
  exit 1
fi
mkdir -p "$OUT"

# 1) CHarm as a DLL: the whole static archive, FFTW (OpenMP build) from MSYS2 as the FFT backend.
gcc $FLAGS -shared -o "$OUT/libcharm.dll" \
  -Wl,--whole-archive "$CHARM/.libs/libcharm.a" -Wl,--no-whole-archive \
  -lfftw3_omp -lfftw3 -lm \
  -Wl,--export-all-symbols -Wl,--out-implib,"$OUT/libcharm.dll.a"

# 2) The driver. The macro strings are recorded in the notes of the result and timings files.
gcc $FLAGS -O3 \
  -DCHARM_DRIVER_FLAGS='"-O3 -march=native -fopenmp (CHarm and driver), FFTW 3.3.11 OpenMP from MSYS2"' \
  -DCHARM_VERSION_TAG='"0.4.11 (commit 1ead5828)"' \
  -I"$CHARM" "$HERE/charm_shs_bench.c" \
  -L"$OUT" -lcharm -lm -o "$OUT/charm_shs_bench.exe"

echo "built: $OUT/libcharm.dll and $OUT/charm_shs_bench.exe"
echo "run with: export PATH=\"$OUT:\$PATH\" before calling charm_shs_bench.exe"
