#!/bin/sh
###General options
### -- set the job Name --
#SBATCH --job-name=Mesh25Co1
#SBATCH --time=48:00:00
#SBATCH --partition=windfatq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH -o log.out
#SBATCH -e log.err

#(
  #blockMesh &&
  #cp -pr 0.org/* 0/ &&
  #waveGaugesNProbes &&
  #setWaveParameters &&
  #setWaveField 
#  relaxationZoneLayout 
#  decomposePar
#) 2>&1 | tee log.init

mpirun -np 20 waveFoam -parallel > log.run
# IsoFoam > log.run
