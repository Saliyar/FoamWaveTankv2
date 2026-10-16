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

The HOS modes file is about 141 MB, so it is not stored in the repository.
Generate it before running:

1. Put the raw HOS-NWT output `modes_HOS_SWENSE.dat` (512 modes, h = 5 m) in `constant/`.
2. Run the cleaning step from the repository:

       python3 doc/swense_port/clean_modes_generator.py constant/modes_HOS_SWENSE.dat constant/modes_HOS_SWENSE_clean.dat

   This zeroes the modes with k_x above KCUT = 7.71 1/m (113 of 511 modes are kept).

The file `constant/Grid2Grid.dict` points to `constant/modes_HOS_SWENSE_clean.dat`.

## Running

    ./Allrun      # blockMesh, then waveIsoFoamSWENSE (serial; the 8-core coarse run took about 20 min, so expect several hours serial, not measured)
    ./Allclean    # removes time directories, mesh and logs

Parallel runs need a `decomposeParDict` and `decomposePar` first; that is not included here.

## Checks

The per-step line `SWENSE check` in the log must show max|alpha1 - alphaI| = 0 (round-off,
criterion 1e-12) and max|Uc| at round-off (criterion 1e-10) throughout the run.
