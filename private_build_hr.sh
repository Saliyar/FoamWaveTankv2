#!/bin/bash -l
# Private build of branch swense-zone-fix with the height ramp: libwaves2Foam,
# libwaves2FoamProcessing, libswenseWaveTank_dyn and the solver go to
# ./private_hr/{lib,bin}; nothing installed is touched.
set -e
WT=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/fwt_zone_fix
OUT=${OUT:-$WT/private_hr}
export PAGER=cat MODULES_PAGER=cat
module use /services/tools/modulefiles
module load tools gcc/14.2.0 gsl/2.7 fftw/3.3.10 ucx/1.20.0 openmpi/5.0.2 cuda/toolkit/11.0.3 openfoam.com/2506
set +e; source /services/tools/openfoam.com/2506/OpenFOAM-v2506/etc/bashrc; set -e
mkdir -p $OUT/lib $OUT/bin
export FOAM_USER_LIBBIN=$OUT/lib FOAM_USER_APPBIN=$OUT/bin
export WAVES_DIR=$WT WAVES_SRC=$WT/src WAVES_LIBBIN=$OUT/lib WAVES_APPBIN=$OUT/bin
export WAVES_GSL_INCLUDE=/services/tools/gsl/2.7/include WAVES_GSL_LIB=/services/tools/gsl/2.7/lib
export WAVES_XVERSION=0 EXTBRANCH=0 OFPLUSBRANCH=1 WM_PROJECT_VERSION_NUMBER=2506
for d in waves2Foam waves2FoamProcessing swenseWaveTank; do
    wmakeLnInclude -u $WT/src/$d > /dev/null
done
cd $WT/src/waves2Foam && wmake libso > $OUT/build_w2f.log 2>&1 || { grep -E "error" $OUT/build_w2f.log | head; exit 1; }
cd $WT/src/waves2FoamProcessing && wmake libso > $OUT/build_w2fp.log 2>&1 || { grep -E "error" $OUT/build_w2fp.log | head; exit 1; }
cd $WT/src/swenseWaveTank && wmake libso > $OUT/build_lib.log 2>&1 || { grep -E "error" -A3 $OUT/build_lib.log | head -40; exit 1; }
cd $WT/applications/solvers/waveIsoFoamSWENSE && wmake > $OUT/build_solver.log 2>&1 || { grep -E "error" -A3 $OUT/build_solver.log | head -40; exit 1; }
grep -E " error|undefined" $OUT/build_*.log | head
md5sum $OUT/lib/*.so $OUT/bin/*
