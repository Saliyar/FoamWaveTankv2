# Dynamic mesh for waveFoamSWENSE and waveIsoFoamSWENSE: report

Status: STOPPED after Stage 3. Static invariance PASSES. Zero-motion check FAILS against the plain reference (alpha and p above 10x). Stage 2 items 7 to 9 are implemented, but the sixDoF free-floating path (Decision 2) is NOT implemented. Stages 4, 5 and 6 are not run.

## Summary

1. Branch SWENSE_dynamic at a078d9b, worktree wt_SWENSE_dynamic. Local commits only, no push, no merge, production binaries and libraries untouched (build names _dyn and libswenseWaveTank_dyn).
2. Build: OpenFOAM v2506 (the production binaries are v2506 too), OFVERSION=2206 unchanged, 0 compile errors.
3. Stage 0 map (00_map.md) and Stage 1 spec (01_spec.md, with the continuity identity from Sithik) are done.
4. Stage 2 implemented: createDynamicFvMesh in both solvers, mesh.update at the first outer corrector (and each corrector with moveMeshOuterCorrectors), relative total flux phi = phiComp + phiInc - meshPhi, relative incident flux phiIRel for alphaI, gh/ghf update, Grid2Grid/Air2Air point rebuild only when points move, psiProj pinned to its own solver dictionary.
5. Stage 2 items 6 to 8 (added after Decision 1): moving body wall velocity Uc = Uwall - UInc with Uwall the normal face velocity from meshPhi; column skip for any column with a cell moved since t = 0 (log line per step); FatalError if the mesh moves inside a relaxation zone (weight above zero).
6. Static invariance (staticFvMesh, 20 steps): waveFoamSWENSE tutorial, waveIsoFoamSWENSE tutorial and box L1_sw are bit-identical to production a078d9b (max diff 0 in alpha, Uc, pc, U, p). PASS.
7. Zero motion (dynamicMotionSolverFvMesh, displacementLaplacian, zero displacement) vs plain solvers on the box case: SWENSE alpha 13.6x and p 31x the plain own static-vs-zero difference (FAIL, limit 10x); U 4.4x (PASS). Step-1 difference already visible (alpha1 - alphaI 2.6e-14 vs 2.2e-16). Term search (Decision 1) NOT done.
8. Decision 2 (sixDoFRigidBodyMotion, free floating, in both SWENSE solvers) NOT implemented. Stage 6 not run.

## Stage results

| item | status | numbers |
|---|---|---|
| 0 map | done | 00_map.md |
| 1 spec | done | 01_spec.md |
| 2 implement (prescribed motion path) | done | build 0 errors |
| 2 items 6-8 (body velocity, column skip, zone check) | done, not exercised by any run | - |
| 3 static invariance | PASS | bit-identical, 3 cases |
| 3 zero motion vs production | not a pass criterion now | box alpha 6.9e-13, U 1.76e-11, p 1.8e-8 |
| 3 zero motion vs plain (Decision 1) | FAIL | see table below |
| 4 moving mesh pure waves | not run | - |
| 5 prescribed body | not run | - |
| 6 free floating (sixDoF) | not implemented | - |

Zero-motion vs plain comparison, final time t = 0.118 (20 steps, deltaT 0.0059), absolute max differences, box L1_sw case:

| field | SWENSE static vs zero | plain waveFoam static vs zero | plain waveIsoFoam static vs zero | SWENSE / plain waveFoam | limit (10x plain) | result |
|---|---|---|---|---|---|---|
| alpha.water | 6.93e-13 | 5.11e-14 | 7.02e-14 | 13.6 | 5.1e-13 | FAIL |
| U | 1.76e-11 | 4.00e-12 | 4.42e-12 | 4.4 | 4.0e-11 | PASS |
| p | 1.80e-8 | 5.82e-10 | 4.44e-10 | 31 | 5.8e-9 | FAIL |

Growth over the 20 steps was not checked: only the final time was written. Open issue.

Coarse tutorials (no plain reference), zero vs static with the new build (equal to zero vs production because static is bit-identical):

| case | alpha.water | U | pc | p |
|---|---|---|---|---|
| waveFoamSWENSE coarse | 1.1e-12 | 9.4e-14 (abs) | 3.7e-11 | 3.6e-11 |
| waveIsoFoamSWENSE coarse | 2.2e-10 | 3.1e-13 | 0 | 3.6e-11 |

Plain reference setup (documented deviation): the plain solvers cannot read SWENSE BCs, so the box fields were copied to a scratch case with swenseWaveVelocity and swenseWaveAlpha replaced by zeroGradient, Uc, pc and alphaI removed, and an empty relaxationNames list added to waveProperties. Plain waveFoam/waveIsoFoam runs therefore have no wave forcing. They serve only as a mesh-motion reference.

## Deviations and environment

- Coarse tutorials are stale: their fvSolution lacks psiProj and cAlpha entries that the production binaries require. Runs used the box case fvSolution and fvSchemes (accepted deviation, see Decision 1 in the brief). The tutorials are NOT edited.
- HOS modes file for the coarse tutorials linked from SWENSE_testing diag new10 (not in the repository).
- The plain reference has no wave forcing (see above).

## Runs (wall times in seconds, serial)

| label | solver | case | mode | wall |
|---|---|---|---|---|
| s3_wfs_prod | waveFoamSWENSE prod | coarse | prod | 9 |
| s3_wfs_static | waveFoamSWENSE_dyn | coarse | static | 8 to 9 |
| s3_wfs_zero | waveFoamSWENSE_dyn | coarse | zero | 12 to 13 |
| s3_wis_prod | waveIsoFoamSWENSE prod | coarse | prod | 8 |
| s3_wis_static | waveIsoFoamSWENSE_dyn | coarse | static | 7 to 8 |
| s3_wis_zero | waveIsoFoamSWENSE_dyn | coarse | zero | 10 to 12 |
| s3_box_prod | waveFoamSWENSE prod | box | prod | 8 |
| s3_box_static | waveFoamSWENSE_dyn | box | static | 8 to 9 |
| s3_box_zero | waveFoamSWENSE_dyn | box | zero | 14 to 15 |
| p_wf_static | waveFoam prod | box (plain copy) | static | 3 |
| p_wf_zero | waveFoam prod | box (plain copy) | zero | 9 |
| p_wis_static | waveIsoFoam prod | box (plain copy) | static | 3 |
| p_wis_zero | waveIsoFoam prod | box (plain copy) | zero | 9 |

Earlier failed invocations (missing dictionary entries, wrong motion dictionary syntax, a run started before the build finished) are not results and are not counted.

## Open issues

1. Decision 1 term search (not done): switch one term at a time over 20 steps, in order: meshPhi (set to zero even when moving), gh/ghf recompute, the V0/V factor of the moved-mesh time derivative, the Ucf/phiComp rebuild, the phiInc projection. Step-1 divergence in alpha1 - alphaI (2.6e-14) points at the alpha step or the mesh-flux path.
2. Growth check over 20 steps needs fields at intermediate times (writeInterval 1).
3. Decision 2 (sixDoFRigidBodyMotion for waveFoamSWENSE and waveIsoFoamSWENSE, reading total p and U, rho by name) not implemented. OpenFOAM v2506 has the library (src/sixDoFRigidBodyMotion). Do not use lib6DoFDynamicMeshForcedMotion.so or FloatStepper libraries.
4. Grid2Grid setupPoints re-call on moved points is an API call, not verified against the "once per mesh" comment in Grid2Grid.H (no change inside Grid2Grid made).
5. Body velocity (item 6) and column skip (item 7) are not exercised: no moving body or moving columns in any run yet.
6. The alpha equation has no explicit V0/V factor for the moved mesh (spec item 3 requires it); the SWENSE VoF advection uses the standard advective form. This may be the source of the alpha discrepancy.
