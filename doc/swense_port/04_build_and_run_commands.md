# Phase D: build and run commands (exact)

Status: the library and both solvers are built with wmake; the top-level
Allwmake is not run (it updates git submodules). Run on the cluster login node
gbarlogin1; all paths are absolute.

## 1. Environment (every build and run shell)

```bash
source /etc/profile.d/modules.sh
module load OpenFoam/v2412/gcc-12.4.0-binutils-2.42-openmpi-5.0.7
source /zhome/19/7/191420/OpenFOAM/sitali-v2412/FoamWaveTankv2/bin/bashrc x
```

FFTW: the cluster module `fftw3/3.3.10-openmpi-5.0.6-gcc-12.4.0` (the only
gcc 12.4 FFTW) cannot be loaded together with the OpenFOAM module above, because
it requires mpi/5.0.6 and OpenFOAM loads mpi/5.0.7. The module's FFTW library has
no MPI dependency (ldd shows only libc and libm), and it exports the Fortran
interface (dfftw_*). So its library directory
`/appl9/fftw/3.3.10-openmpi-5.0.6-gcc-12.4.0/EPYC9354/gcc-12.4.0/lib` is used
directly, as RPATH (`-Wl,--disable-new-dtags -Wl,-rpath,<dir>`) in libGrid2Grid.so
and in each solver. No `LD_LIBRARY_PATH` change is needed. The OpenFOAM module
puts the ThirdParty FFTW (without Fortran symbols) on LD_LIBRARY_PATH; the RPATH
(DT_RPATH) takes precedence over it.

## 2. Grid2Grid for v2412 (done once; re-run only if the library is lost)

```bash
bash /zhome/19/7/191420/OpenFOAM/sitali-v2412/FoamWaveTankv2/doc/swense_port/grid2grid_build.sh \
     /tmp/claude-291420/-zhome-19-7-191420-Softwares/2527d8ca-d112-4b82-9e7f-320845e881a7/scratchpad/g2g_build
```

The script loads the modules itself, compiles bspline and fymc (static), then
Grid2Grid, and installs libGrid2Grid.so into the v2412 user library. Verified:
the rebuilt file is byte-identical to the installed one.

## 3. SWENSE library and solvers (wmake, not Allwmake)

```bash
cd /zhome/19/7/191420/OpenFOAM/sitali-v2412/FoamWaveTankv2/src/swenseWaveTank
wmake libso

cd /zhome/19/7/191420/OpenFOAM/sitali-v2412/FoamWaveTankv2/applications/solvers/waveFoamSWENSE
wmake

cd /zhome/19/7/191420/OpenFOAM/sitali-v2412/FoamWaveTankv2/applications/solvers/waveIsoFoamSWENSE
wmake
```

Installed to `$FOAM_USER_LIBBIN` (libswenseWaveTank.so) and `$FOAM_USER_APPBIN`
(waveFoamSWENSE, waveIsoFoamSWENSE). Compiler: gcc 12.4.0 (module).

Library dependencies, from ldd (module environment, no LD_LIBRARY_PATH change):

- libGrid2Grid.so -> libfftw3.so.3 (module lib dir), libgfortran.so.5 (gcc 12.4)
- waveFoamSWENSE -> libwaves2Foam.so, libswenseWaveTank.so, libGrid2Grid.so,
  libfftw3.so.3 (module lib dir)
- waveIsoFoamSWENSE -> libgeometricVoF.so, libswenseWaveTank.so,
  libGrid2Grid.so, libfftw3.so.3 (module lib dir)

## 4. Run (per case directory, job script = bash.sh)

Each run directory holds the same inputs. Submit with LSF:

```bash
cd /work3/sitali/SWENSEComparison/RegularWaveTesting/FWT_swense_regular/run_waveFoamSWENSE
bsub < bash.sh          # 16 ranks: decomposePar, then mpirun -np 16 waveFoamSWENSE -parallel
```

The job script template is `doc/swense_port/job_template_waveFoamSWENSE.sh`.
The isoAdvector run uses the same script with waveIsoFoamSWENSE.

## 5. Not run

Top-level `Allwmake` (it runs `git submodule update --init --remote`, which would
move the checked-out submodules). Its solver loop already matches the new
`waveFoamSWENSE` and `waveIsoFoamSWENSE` directories.
