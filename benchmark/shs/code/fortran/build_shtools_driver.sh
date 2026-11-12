#!/bin/bash
# Build SHTOOLS 4.14.1 (Fortran 95, OpenMP) and the shtools_shs_bench driver on MSYS2 UCRT64.
#
# Run from the MSYS2 UCRT64 shell (C:\msys64\ucrt64.exe), with the packages listed in README.md:
#   bash benchmark/shs/code/fortran/build_shtools_driver.sh
#
# Nothing is installed system-wide and nothing is written outside build/external/ (git-ignored).
# SHTOOLS is fetched anonymously from its public GitHub repository (no credentials are used).
# FFTW (GPL) is only linked here as an external build dependency; nothing is redistributed.
set -euo pipefail

UCRT=/c/msys64/ucrt64; [ -x "$UCRT/bin/gfortran.exe" ] || UCRT=/ucrt64   # MSYS2 UCRT64 prefix
export PATH="$UCRT/bin:/usr/bin:$PATH"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
EXT="$REPO/build/external"
SRC="$EXT/shtools"
BUILD="$EXT/shtools_bench"
SHTOOLS_TAG=v4.14.1
FFLAGS_LIB="-fopenmp -fPIC -O3 -march=native -std=gnu"
FFLAGS_DRV="-fopenmp -O3 -march=native -std=f2008"
LIB_FLAGS_TEXT="SHTOOLS 4.14.1 $FFLAGS_LIB (no fast-math)"

# Work around a gfortran temp-file failure seen when gfortran runs under make in this MSYS2
# install ("Cannot create temporary file in C:\WINDOWS\"). Windows-style paths fix it.
mkdir -p "$BUILD/tmp"
WINTMP="$(cygpath -w "$BUILD/tmp")"
export TMP="$WINTMP" TEMP="$WINTMP" TMPDIR="$WINTMP"

# 1. Sources: SHTOOLS $SHTOOLS_TAG (anonymous clone; never prompts for credentials).
if [ ! -f "$SRC/src/Makefile" ]; then
    mkdir -p "$EXT"
    GIT_TERMINAL_PROMPT=0 GCM_INTERACTIVE=never \
        git clone --depth 1 --branch "$SHTOOLS_TAG" https://github.com/SHTOOLS/SHTOOLS.git "$SRC"
fi
echo "SHTOOLS source: $(git -C "$SRC" log -1 --format='%H %s')"

# 2. Library: the object list and order of SHTOOLS src/Makefile (target fortran-mp), compiled here
#    directly because 'make fortran-mp' fails in this setup (see the workaround above).
#    The list is read from the Makefile itself so that it cannot drift from the release.
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

# 3. Driver: copy the source, record the library flags in the copy, compile and link.
cd "$BUILD"
sed "s|not recorded (compiled without build_shtools_driver.sh)|$LIB_FLAGS_TEXT|" \
    "$HERE/shtools_shs_bench.f90" > "$BUILD/shtools_shs_bench.f90"
gfortran $FFLAGS_DRV -I "$SRC/include" -c "$BUILD/shtools_shs_bench.f90" -o "$BUILD/shtools_shs_bench.o"
gfortran -fopenmp "$BUILD/shtools_shs_bench.o" -o "$BUILD/shtools_shs_bench.exe" \
    -L"$SRC/lib" -lSHTOOLS-mp -llapack -lblas -lfftw3 -lm
echo "driver: $BUILD/shtools_shs_bench.exe"
