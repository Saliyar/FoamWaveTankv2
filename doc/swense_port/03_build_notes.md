# Phase C build notes: swenseWaveFoam

Status: built and checked on the incident fields. Not yet run as a flow case
(that is Phase D). The full top-level Allwmake has not been run (see below).

## Grid2Grid for OpenFOAM v2412

FoamWaveTankv2 had no libGrid2Grid.so for v2412. It was built from a scratch
copy of ThirdParty/Grid2Grid (the repository copy was not edited):

- Compiler: gfortran 11.5 (the same toolchain as the v2412 build in this shell).
- Fortran sources compiled in the order of src/libGrid2Grid/CMakeLists.txt, with
  -cpp -fPIC -O2, then linked with -shared against libbspline.a, libfymc.a and
  the system libfftw3.so.3.
- CMake's Fortran dependency step failed in the scratch copy (module copy), so the
  objects were built by hand. The symbols that FoamWaveTank's Grid2Grid wrapper
  uses are exported (checked with nm -D).
- Installed to: /zhome/19/7/191420/OpenFOAM/sitali-v2412/platforms/linux64GccDPInt32Opt/lib/libGrid2Grid.so

Runtime FFTW: the Grid2Grid library needs the Fortran FFTW symbols (dfftw_*).
The ThirdParty FFTW in the v2412 environment does not export them, and it comes
first on LD_LIBRARY_PATH. The system /usr/lib64/libfftw3.so.3 does export them.
Runs and tests therefore need:

    export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

after sourcing the OpenFOAM and FoamWaveTankv2 bashrc files.

## Library and solver

- src/swenseWaveTank: libswenseWaveTank (Air2Air copied from the author's own
  files, plus swenseIncident, swenseSurface, swenseRelaxation and the two wave
  boundary conditions swenseWaveVelocity and swenseWaveAlpha).
- applications/solvers/swenseWaveFoam: solver. Uses OpenFOAM v2412 VoF alphaEqn
  (alphaEqnSubCycle.H), pimpleControl, constrainHbyA, constrainPressure and
  adjustPhi; the Grid2Grid forcing is created with FoamWaveTankv2's
  createExternalWaveForcing.H.
- The solver links -lwaves2Foam (for Grid2Grid, convexPolyhedral) and -lfftw3.

Not run: the top-level Allwmake. It starts with `git submodule update --init
--remote`, which would move the checked-out submodules to their latest remote
commits. The two Allwmake files were edited to register the new library and
solver, and the library and solver were built with wmake.

## pc solver settings in the reference run (Generaltesting_order3)

Solver lines in log.foamStarSwense are "GAMG: Solving for pc", so pc was solved
by GAMG. The fvSolution entry that matches pc is the regular expression
"(p_rgh|pc)" (not "pcorr.*", which is PCG and does not match pc):
GAMG, smoother DIC, agglomerator faceAreaPair, mergeLevels 1,
nCellsInCoarsestLevel 10, cacheAgglomeration true, tolerance 1e-8, relTol 0.
The final corrector uses "(p_rghFinal|pcFinal)", which copies p_rgh with relTol 0.
The current fvSolution file was modified after the log was written, so these
settings are the ones in the file now; they match the log's solver names.

## Cut-cell check (wetted fraction)

The wetted cell fraction alphaInc comes from swenseSurface (convexPolyhedral with
the surface d = z - eta). Checks on the Generaltesting_order3 mesh, t = 3 s:

- Flat surface (test-only subclass, z = mid-height of the cell): max error
  2.6e-14. The cut-cell routine is exact for a plane.
- Curved wave surface, 9 cells near the interface and in the bulk: compared with a
  column integral over the cell bounding box. The column integral converges
  (identical at 80, 200 and 600 columns per direction). The largest difference is
  1.8e-4 in alphaInc, at an interface cell, with a consistent positive bias.
  This is larger than round-off and is reported, not corrected.

## Incident values (water, air, interface faces)

Compared with an independent evaluation of the spec formulas with the same
Grid2Grid object and an independently built Air2AirWrapper:

- Water cells: |dU| = 0, relative |dpInc| = 0.
- Air cells: |dU| = 0, relative |dpInc| = 0.
- Interface faces (owner and neighbour of different class): relative |dpfInc| = 0.

These agree exactly because both sides call the same functions with the same
inputs. They check the wrapper logic (hydrostatic term, frame, interface rule),
not the Grid2Grid or Air2Air numbers themselves.

## Conventions used in the port

- Grid2Grid frame x_g = R(theta)(x - referencePosition), theta = referenceTheta.
  foamStar's wave frame with refDirection d equals referenceTheta = -atan2(d_y, d_x).
- Grid2Grid time is t + referenceTime. The Air2Air time uses the same offset.
  foamStar's wave time is t - refTime, so referenceTime = -refTime.
- Grid2Grid Tsoft must be 0 (the solver stops otherwise).
- Water hydrostatic term uses z_w = (x - referencePosition).z.
