#!/bin/bash
# Build libGrid2Grid.so for OpenFOAM v2412 with the gcc 12.4 toolchain and the
# module FFTW (Fortran interface), and install it into the v2412 user library.
#
# Usage:  bash doc/swense_port/grid2grid_build.sh <scratch-dir>
#
# The ThirdParty/Grid2Grid sources are copied to <scratch-dir>; the repository
# copy is not modified. bspline and fymc (static) are built first, then the
# Grid2Grid Fortran modules and the shared library, in the order used by
# ThirdParty/Grid2Grid/src/libGrid2Grid/CMakeLists.txt.
set -e

SCRATCH=${1:?usage: grid2grid_build.sh <scratch-dir>}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
SRC="$REPO/ThirdParty/Grid2Grid"

source /etc/profile.d/modules.sh
module load OpenFoam/v2412/gcc-12.4.0-binutils-2.42-openmpi-5.0.7

# FFTW of the cluster module (gcc 12.4). Its lib directory is also the RPATH.
FFTW_LIB=/appl9/fftw/3.3.10-openmpi-5.0.6-gcc-12.4.0/EPYC9354/gcc-12.4.0/lib
INSTALL=$FOAM_USER_LIBBIN

rm -rf "$SCRATCH" && mkdir -p "$SCRATCH/mod" "$SCRATCH/lib"
(cd "$SRC" && tar --exclude=./build --exclude=./lib -cf - .) | tar -xf - -C "$SCRATCH"

F="-cpp -fPIC -O2 -J$SCRATCH/mod -I$SCRATCH/mod"
cd "$SCRATCH"

# 1. bspline (static): dependency order sub -> oo -> module
for f in bspline_sub_module bspline_oo_module bspline_module; do
    gfortran $F -c src/libBspline/src/$f.f90 -o $f.o
done
ar rcs lib/libbspline.a bspline_sub_module.o bspline_oo_module.o bspline_module.o

# 2. fymc (static)
for f in src/libFyMc/src/mfpGlobal.f90 src/libFyMc/src/general/mfpGeneral.f90 src/libFyMc/src/geometry/mfpGeometry.f90; do
    o=$(basename $f .f90).o
    gfortran $F -c $f -o $o
    objs="$objs $o"
done
ar rcs lib/libfymc.a $objs

# 3. Grid2Grid modules and shared library
G=src/libGrid2Grid
cd $G
for f in modGrid2GridGlobal.f90 fftw3/modFourier_r2c_FFTW3.f90 HDF5/modHDF5Interface.f90 \
         surf2Vol/modSurf2Vol.f90 vol2Vol/modVol2Vol.f90 postGrid2Grid/modPostGrid2Grid.f90 modGrid2Grid.f90; do
    gfortran $F -I$SCRATCH/lib -c $f -o $SCRATCH/g_$(basename $f .f90).o
done
cd "$SCRATCH"
gfortran -shared -o libGrid2Grid.so g_*.o -L"$SCRATCH/lib" -lbspline -lfymc \
    -L"$FFTW_LIB" -lfftw3 -Wl,--disable-new-dtags -Wl,-rpath,"$FFTW_LIB"

cp libGrid2Grid.so "$INSTALL/libGrid2Grid.so"
echo "installed $INSTALL/libGrid2Grid.so"
