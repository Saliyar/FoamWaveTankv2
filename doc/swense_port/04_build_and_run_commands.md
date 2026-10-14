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

FFTW: the system library `/usr/lib64/libfftw3.so.3` is used (decided with Sithik).
The cluster module `fftw3/3.3.10-openmpi-5.0.6-gcc-12.4.0` could not be used: it
requires mpi/5.0.6 (OpenFOAM loads 5.0.7), and its library contains AVX-512 code
that crashes on the XeonE5-2660v3 nodes. The system build exports the Fortran
interface (dfftw_*) and contains no AVX-512 instructions. libGrid2Grid.so carries
`/usr/lib64` as RPATH (`-Wl,--disable-new-dtags -Wl,-rpath,/usr/lib64`), so no
LD_LIBRARY_PATH change is needed. The solvers do not link FFTW directly; they get it
through libGrid2Grid.so (`-Wl,--unresolved-symbols=ignore-in-shared-libs`), so the
ThirdParty FFTW on LD_LIBRARY_PATH is never loaded.

Portable build: the OpenFOAM module sets `-march=native` (login-node CPU). The
environment used for builds overrides FOAM_EXTRA_CXXFLAGS with `-march=x86-64-v3`,
which runs on the XeonE5-2660v3 and EPYC nodes. Checked: libwaves2Foam,
libwaves2FoamSampling, libswenseWaveTank, libGrid2Grid and both solvers contain no
AVX-512 instructions (objdump).

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
