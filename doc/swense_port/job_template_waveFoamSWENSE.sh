#!/bin/sh
### General options (same style as the reference case)
#BSUB -J fwt_swense_waveFoamSWENSE
#BSUB -n 16
#BSUB -R "span[ptile=4]"
#BSUB -W 4:00
#BSUB -o Output_%J.out
#BSUB -e Error_%J.err

# Load modules: v2412 OpenFOAM (gcc 12.4 toolchain, matches the FFTW build below)
module purge
module load OpenFoam/v2412/gcc-12.4.0-binutils-2.42-openmpi-5.0.7

# FFTW: module fftw3/3.3.10-openmpi-5.0.6-gcc-12.4.0 cannot be loaded together
# with this OpenFOAM module (MPI 5.0.6 vs 5.0.7). Its library directory is
# built into libGrid2Grid.so and waveFoamSWENSE as RPATH, so no LD_LIBRARY_PATH change is needed.

# Decompose and run on 16 ranks
decomposePar -force > log.decompose 2>&1
mpirun -np 16 waveFoamSWENSE -parallel > log.foam 2>&1
