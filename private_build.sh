#!/bin/bash -l
# Private build of the momentum-blend SWENSE (branch swense-momentum-rev1).
# Library and solver go to ./private/{lib,bin}; nothing installed is touched.
set -e
WT=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/fwt_momentum_rev1
export PAGER=cat MODULES_PAGER=cat
module use /services/tools/modulefiles
module load tools gcc/14.2.0 gsl/2.7 fftw/3.3.10 ucx/1.20.0 openmpi/5.0.2 cuda/toolkit/11.0.3 openfoam.com/2506
set +e; source /services/tools/openfoam.com/2506/OpenFOAM-v2506/etc/bashrc; set -e
mkdir -p $WT/private/lib $WT/private/bin
export FOAM_USER_LIBBIN=$WT/private/lib
export FOAM_USER_APPBIN=$WT/private/bin
export WAVES_DIR=$WT WAVES_SRC=$WT/src WAVES_LIBBIN=$WT/private/lib WAVES_APPBIN=$WT/private/bin
echo "FOAM_USER_LIBBIN=$FOAM_USER_LIBBIN"; echo "FOAM_USER_APPBIN=$FOAM_USER_APPBIN"
# waves2Foam: library reused from the previous private build (source identical), headers linked here
[ -f $WT/private/lib/libwaves2Foam.so ] || cp -p /gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/fwt_absolute_pressure/private/lib/libwaves2Foam.so $WT/private/lib/
[ -d $WT/src/waves2Foam/lnInclude ] || wmakeLnInclude $WT/src/waves2Foam
cd $WT/src/swenseWaveTank && wmake libso 2>&1 | tee $WT/build_rev1_lib.log | grep -E "error|Error|undefined|^wmake|\.so$" | tail -40
cd $WT/applications/solvers/waveIsoFoamSWENSE && wmake 2>&1 | tee $WT/build_rev1_solver.log | grep -E "error|Error|undefined|^wmake|_dyn$" | tail -60
ls -la $WT/private/lib $WT/private/bin
md5sum $WT/private/lib/*.so $WT/private/bin/*
