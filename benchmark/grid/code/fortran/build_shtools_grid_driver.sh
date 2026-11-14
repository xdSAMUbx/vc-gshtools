#!/usr/bin/env bash
# Builds SHTOOLS 4.14.1 (Fortran 95, OpenMP library) and the shtools_grid_bench driver on Linux or WSL.
#
# Run from the repository root, in WSL (Ubuntu or Debian):
#     bash benchmark/grid/code/fortran/build_shtools_grid_driver.sh
#
# Prerequisites (not installed by this script): gfortran, git, make-free build (this script compiles the object
# list read from SHTOOLS src/Makefile itself), liblapack-dev, libblas-dev, libfftw3-dev.
#
# SHTOOLS is cloned anonymously into build/external/shtools-wsl (git-ignored). The Windows/MSYS2 tree in
# build/external/shtools is NOT used, because its objects are for MSYS2.
# Nothing is written outside build/external. This script was written but NOT run by its author.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
EXT="$REPO/build/external"
SRC="$EXT/shtools-wsl"
BUILD="$EXT/shtools-grid-wsl"
SHTOOLS_TAG=v4.14.1
FFLAGS_LIB="-fopenmp -fPIC -O3 -march=native -std=gnu"
FFLAGS_DRV="-fopenmp -O3 -march=native -std=f2008"

# 1. Sources: SHTOOLS v4.14.1 (anonymous clone; never prompts for credentials).
if [ ! -f "$SRC/src/Makefile" ]; then
    mkdir -p "$EXT"
    GIT_TERMINAL_PROMPT=0 GCM_INTERACTIVE=never \
        git clone --depth 1 --branch "$SHTOOLS_TAG" https://github.com/SHTOOLS/SHTOOLS.git "$SRC"
fi
echo "SHTOOLS source: $(git -C "$SRC" log -1 --format='%H %s')"

# 2. Library: the object list and order of SHTOOLS src/Makefile (target fortran-mp), compiled here directly.
grab() { sed -n "/^$1 = /,/^\$/p" "$SRC/src/Makefile" | sed "s/^$1 = //" | tr -d '\r\\' | tr '\n\t' '  '; }
S0="$(grab SRCS0)"; SL="$(grab SRCSLAPACK)"; SF="$(grab SRCSFFTW)"
MODS="ftypes.f95 SHTOOLS.f95 FFTW3.f95 PlanetsConstants.f95"

mkdir -p "$SRC/include" "$SRC/lib"
cd "$SRC/src"
objs=""
for f in $MODS; do
    gfortran -c $FFLAGS_LIB -J "$SRC/include" "$f"; objs="$objs ${f%.*}.o"
done
for f in $S0 $SF; do
    gfortran -c $FFLAGS_LIB -J "$SRC/include" "$f"; objs="$objs ${f%.*}.o"
done
for f in $SL; do
    gfortran -c $FFLAGS_LIB -J "$SRC/include" -cpp "$f"; objs="$objs ${f%.*}.o"
done
rm -f "$SRC/lib/libSHTOOLS-mp.a"
ar rcs "$SRC/lib/libSHTOOLS-mp.a" $objs
ranlib "$SRC/lib/libSHTOOLS-mp.a"
echo "library: $SRC/lib/libSHTOOLS-mp.a ($(echo $objs | wc -w) objects)"

# 3. Driver: compile and link against the library.
mkdir -p "$BUILD"
gfortran $FFLAGS_DRV -I "$SRC/include" -c "$HERE/shtools_grid_bench.f90" -o "$BUILD/shtools_grid_bench.o"
gfortran -fopenmp "$BUILD/shtools_grid_bench.o" -o "$BUILD/shtools_grid_bench" \
    -L"$SRC/lib" -lSHTOOLS-mp -llapack -lblas -lfftw3 -lm
echo "driver: $BUILD/shtools_grid_bench"
echo "run example (from the repository root; --out-dir must exist):"
echo "  mkdir -p build/grid_out && $BUILD/shtools_grid_bench --gfc benchmark/shs/data/XGM2019e_2190_tide_free.gfc \\"
echo "      --nmax 2190 --case dh1 --threads 1 --repeat 3 --out-dir build/grid_out"
