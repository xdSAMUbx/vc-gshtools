#!/usr/bin/env bash
# Builds the charm_grid_bench driver on Linux or WSL (Ubuntu or Debian packages), against the CHarm 0.4.11 static
# library already built in build/external/charm (see benchmark/shs/code/c/README.md for how CHarm is built).
#
# Run from the repository root, in WSL:
#     bash benchmark/grid/code/c/build_charm_grid_driver.sh
#
# Prerequisites (not installed by this script):
#   gcc with OpenMP (libgomp), libfftw3-dev with the OpenMP FFTW library (libfftw3_omp), libm.
#   build/external/charm/.libs/libcharm.a must exist.
#
# Nothing is written outside build/external/charm-grid (git-ignored under build/).
# This script was written but NOT run by its author; see benchmark/grid/code/README.md, section "Not verified".
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
CHARM="$REPO/build/external/charm"
OUT="$REPO/build/external/charm-grid"
FLAGS="-O3 -march=native -fopenmp"

if [ ! -f "$CHARM/.libs/libcharm.a" ]; then
  echo "error: $CHARM/.libs/libcharm.a not found; build CHarm first (benchmark/shs/code/c/README.md)" >&2
  exit 1
fi
mkdir -p "$OUT"

# The driver. The version tag is recorded in its stdout line. Link order: the driver, the static CHarm
# archive, then the FFT backend (CHarm's grid FFT path needs FFTW with OpenMP) and libm.
gcc $FLAGS -Wall -Wextra -std=gnu11 \
  -DCHARM_VERSION_TAG='"0.4.11 (commit 1ead5828)"' \
  -I"$CHARM" "$HERE/charm_grid_bench.c" \
  "$CHARM/.libs/libcharm.a" -lfftw3_omp -lfftw3 -lm \
  -o "$OUT/charm_grid_bench"

echo "built: $OUT/charm_grid_bench"
echo "run example (from the repository root; --out-dir must exist):"
echo "  mkdir -p build/grid_out && $OUT/charm_grid_bench --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \\"
echo "      --nmax 2190 --case dh1 --threads 8 --repeat 3 --out-dir build/grid_out"
