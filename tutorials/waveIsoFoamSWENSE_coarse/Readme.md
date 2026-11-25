# waveIsoFoamSWENSE_coarse

Coarse SWENSE regular-wave case with isoAdvector alpha transport, solver `waveIsoFoamSWENSE`.

- **Mesh:** 16262 cells (346 x 47 in x and z, one cell in y), made by `blockMesh`
  from `system/blockMeshDict` (free-surface band fsmin..fsmax = -0.29..0.29 m, 24 cells).
- **Time:** deltaT 0.0059 s, endTime 40 s (6780 steps), writeInterval at 40 s.
- **Grid2Grid:** S1 settings in `constant/Grid2Grid.dict` (zMax 1.0, nZMax 50).
- **Alpha transport:** isoAdvector for alpha1 and alphaI (both isoAdvector, same function).
  waveIsoFoamSWENSE uses the same settings for the incident fraction alphaI.
- **Gauges:** surface elevation from alpha at x = -8, -1, 0, 1, 8 m (`functions` in `system/controlDict`).

## Modes file (not in the repository)

The HOS modes file is about 141 MB, so it is not stored in the repository
(`.gitignore` skips `constant/modes_HOS_SWENSE*.dat`). Make it before running:

1. Put the raw HOS-NWT output `modes_HOS_SWENSE.dat` (512 modes, h = 5 m,
   147603456 bytes, md5 `cb788107a479f03eb755d524d6760e65`) in `constant/`.
   It is the HOS-NWT run of this wave; ask the maintainers for a copy.
2. From this case folder run (Python 3 with numpy, about 15 s):

       python3 tools/clean_modes_generator.py constant/modes_HOS_SWENSE.dat constant/modes_HOS_SWENSE_clean.dat

   This zeroes the modes with k_x above KCUT = 7.71 1/m (113 of 511 modes are kept).
   The result must have md5 `d7e765c71c387fc5bc411ed7a1e62aec`.

The file `constant/Grid2Grid.dict` points to `constant/modes_HOS_SWENSE_clean.dat`.

## Settings

The case files hold the validated pure-wave settings; the code defaults are different:

- `constant/waveProperties`: `pcAtmosphereConsistent yes`, `incidentHydrostatic discrete`,
  `alphaHExact true` (relaxation target = exact wetted fraction), `incidentResidualMax 0.02`.
- `system/controlDict`: `alphaIExactReset yes`. After the relaxation, alphaI is set to the
  same exact wetted fraction outside `plainBodyZone`, and alpha1 takes the same increment.
  It needs `alphaHExact true` (the solver stops otherwise) and replaces `levelCorrection`.
- `system/fvSolution`: isoAdvector tolerances 1e-8, `PIMPLE nOuterCorrectors 1`, `psiProj` solver.

## Running

    ./Allrun      # blockMesh, then waveIsoFoamSWENSE (serial; the 8-core coarse run took about 20 min, so expect several hours serial, not measured)
    ./Allclean    # removes time directories, mesh and logs

Parallel runs need a `decomposeParDict` and `decomposePar` first; that is not included here.

## Checks

The per-step line `SWENSE check` in the log must show max|alpha1 - alphaI| = 0 (round-off,
criterion 1e-12) and max|Uc| at round-off (criterion 1e-10) throughout the run. The line
`SWENSE alphaIExactReset` gives the reset volume of each step (net and abs, m3); it should
stay small and steady.
