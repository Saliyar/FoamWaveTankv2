#!/bin/bash
R=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/wt_SWENSE_dynamic/diag/dynamic_mesh/scripts/run_case.sh
S=/tmp/claude-7201142/-gpfs-hpc-home-projects-dtu-00075-apps-FoamWaveTank-FoamWaveTankv2/001c5316-aea6-4d14-b3d8-47c9f34f36fc/scratchpad/stage3
B=/home/people/sitali/OpenFOAM/sitali-v2506/platforms/linux64GccDPInt64Opt/bin
T=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/wt_SWENSE_dynamic/tutorials
L=/home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/box_fixed/coarsen3/L1_sw
export MODES=/home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/new10/constant/modes_HOS_SWENSE_clean.dat
export FVSOL=/home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/box_fixed/coarsen3/L1_sw/system/fvSolution
export FVSCHEMES=/home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/box_fixed/coarsen3/L1_sw/system/fvSchemes
$R s3_wfs_prod $B/waveFoamSWENSE $T/waveFoamSWENSE_coarse prod 20 $S
$R s3_wfs_static $B/waveFoamSWENSE_dyn $T/waveFoamSWENSE_coarse static 20 $S
$R s3_wfs_zero $B/waveFoamSWENSE_dyn $T/waveFoamSWENSE_coarse zero 20 $S
$R s3_wis_prod $B/waveIsoFoamSWENSE $T/waveIsoFoamSWENSE_coarse prod 20 $S
$R s3_wis_static $B/waveIsoFoamSWENSE_dyn $T/waveIsoFoamSWENSE_coarse static 20 $S
$R s3_wis_zero $B/waveIsoFoamSWENSE_dyn $T/waveIsoFoamSWENSE_coarse zero 20 $S
$R s3_box_prod $B/waveFoamSWENSE $L prod 20 $S
$R s3_box_static $B/waveFoamSWENSE_dyn $L static 20 $S
$R s3_box_zero $B/waveFoamSWENSE_dyn $L zero 20 $S
echo STAGE3_DONE
