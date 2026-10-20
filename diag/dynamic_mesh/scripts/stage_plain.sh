#!/bin/bash
R=/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/wt_SWENSE_dynamic/diag/dynamic_mesh/scripts/run_case.sh
S=/tmp/claude-7201142/-gpfs-hpc-home-projects-dtu-00075-apps-FoamWaveTank-FoamWaveTankv2/001c5316-aea6-4d14-b3d8-47c9f34f36fc/scratchpad/stage3
P=/home/people/sitali/OpenFOAM/sitali-v2506/platforms/linux64GccDPInt64Opt/bin
PS=/tmp/claude-7201142/-gpfs-hpc-home-projects-dtu-00075-apps-FoamWaveTank-FoamWaveTankv2/001c5316-aea6-4d14-b3d8-47c9f34f36fc/scratchpad/stage3/plain_src
$R p_wf_static $P/waveFoam $PS static 20 $S
$R p_wf_zero $P/waveFoam $PS zero 20 $S
$R p_wis_static $P/waveIsoFoam $PS static 20 $S
$R p_wis_zero $P/waveIsoFoam $PS zero 20 $S
echo PLAIN_DONE
