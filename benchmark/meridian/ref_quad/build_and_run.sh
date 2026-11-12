#!/usr/bin/env bash
# Builds and runs the binary128 reference (see ref_quad.cpp). Run from WSL (Ubuntu 24.04, gcc >= 13, libquadmath).
#   wsl -d Ubuntu-24.04 -e bash benchmark/meridian/ref_quad/build_and_run.sh [extra ref_quad options]
# Run it from the repository root (a /mnt/c/... path). Output: benchmark/meridian/results/ref_quad_meridian.csv
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
bin="${TMPDIR:-/tmp}/ref_quad"
g++ -O2 -std=gnu++17 -fext-numeric-literals -fopenmp -o "$bin" "$here/ref_quad.cpp" -lquadmath
"$bin" "${GFC:-$root/data/cache/modelos/EGM2008/EGM2008.gfc}" "${OUT:-$root/benchmark/meridian/results/ref_quad_meridian.csv}" "$@"
