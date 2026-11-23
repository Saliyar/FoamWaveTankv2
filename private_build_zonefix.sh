#!/bin/bash -l
# Private build of branch swense-zone-fix. Library and solver go to
# ./private/{lib,bin}; nothing installed is touched.
set -e
WT=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/fwt_zone_fix
E5=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/fwt_rev1_floating/private_e5
export PAGER=cat MODULES_PAGER=cat
module use /services/tools/modulefiles
module load tools gcc/14.2.0 gsl/2.7 fftw/3.3.10 ucx/1.20.0 openmpi/5.0.2 cuda/toolkit/11.0.3 openfoam.com/2506
set +e; source /services/tools/openfoam.com/2506/OpenFOAM-v2506/etc/bashrc; set -e
mkdir -p $WT/private/lib $WT/private/bin
export FOAM_USER_LIBBIN=$WT/private/lib FOAM_USER_APPBIN=$WT/private/bin
export WAVES_DIR=$WT WAVES_SRC=$WT/src WAVES_LIBBIN=$WT/private/lib WAVES_APPBIN=$WT/private/bin
[ -f $WT/private/lib/libwaves2Foam.so ] || cp -p $E5/lib/libwaves2Foam.so $WT/private/lib/
[ -d $WT/src/waves2Foam/lnInclude ] || wmakeLnInclude $WT/src/waves2Foam > /dev/null
cd $WT/src/swenseWaveTank && wmake libso > $WT/private/build_lib.log 2>&1 || { grep -E "error" $WT/private/build_lib.log | head; exit 1; }
cd $WT/applications/solvers/waveIsoFoamSWENSE && wmake > $WT/private/build_solver.log 2>&1 || { grep -E "error" $WT/private/build_solver.log | head; exit 1; }
grep -E " error" $WT/private/build_*.log | head
md5sum $WT/private/lib/*.so $WT/private/bin/*
