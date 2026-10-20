# Dynamic mesh for waveFoamSWENSE and waveIsoFoamSWENSE: report

Status: STOPPED in Stage 3 (static invariance). Stage 3 zero-motion test FAILED. Stages 4 to 6 not run.

1. Branch SWENSE_dynamic at a078d9b, worktree wt_SWENSE_dynamic. No push, no merge, no change to production binaries.
2. Stage 0 map: diag/dynamic_mesh/00_map.md. Stage 1 spec: doc/swense_dynamic/01_spec.md.
3. R3 decision (from Sithik): R3 stays as written, also for deforming meshes. Continuity identity shown in the spec: relative and absolute continuity are the same statement once dV/dt = sum(meshPhi).
4. Stage 2 implemented (createDynamicFvMesh, outer-corrector mesh update, relative total flux phi = phiComp + phiInc - meshPhi, relative phiIRel for alphaI, gh/ghf update, Grid2Grid/Air2Air point rebuild on moving mesh only, psiProj solve pinned to psiProj). Build clean (0 errors).
5. Static invariance: waveFoamSWENSE (tutorial) and box L1_sw are bit-identical to production. waveIsoFoamSWENSE differs at 4e-16 (round-off).
6. Zero-motion (dynamicMotionSolverFvMesh, displacementLaplacian, zero displacement): NOT within 1e-12. Box pc differs by 1.8e-8 absolute (relative 4.7e-7). Differences appear at the first step (alpha1 - alphaI 2.6e-14 vs 2.2e-16).
7. Working hypothesis, not yet tested: round-off in the moved-mesh geometry (V0 versus V factor in the ddt terms), amplified to the pressure solver tolerance. The displacement written at t = 0.118 is exactly zero, so the mesh flux is zero; the difference is not a motion effect.
8. Attempt log for Stage 3 (R5 limit of 3): attempt 1 psiProj final-iteration fix (needed); attempt 2 Grid2Grid setupPoints only when points move (no effect on the numbers). Stop here.
9. Build environment: production binaries are OpenFOAM v2506 (checked with strings and ldd). OFVERSION=2206 kept unchanged in Make/options; build works.
10. Tutorials waveFoamSWENSE_coarse and waveIsoFoamSWENSE_coarse are stale: their fvSolution lacks psiProj and cAlpha entries that the production binaries require. Runs used the box case fvSolution and fvSchemes (superset) and the HOS modes file linked from SWENSE_testing diag new10 (not in repo). This is a deviation from the brief and should be fixed in the tutorials.
11. 6DoF libraries: the repo has no own 6DoF library. The motion libraries in the user platform lib directory (libfloaterMotion, libsixDoFRigidBodyStateWriter, libwaveMakerMotionSolver, libfloaterStateWriter, libmoorDynForFloatStepper) come from the FloatStepper submodule and are out of scope. lib6DoFDynamicMeshForcedMotion.so has no source in the repo; its compiled symbols show forcedRigidBodyMotionSolver.C. Per Sithik, forced motion comes from OpenFOAM, so prescribed motion uses dynamicMotionSolverFvMesh with displacementLaplacian. Free floating (Stage 6) is blocked until Sithik chooses a 6DoF library (OpenFOAM sixDoFRigidBodyMotion is available in v2506).
12. Not done: Stage 4 (moving-mesh pure waves, M1 and M2), Stage 5 (prescribed body motion; the body velocity extension of swenseBodyVelocity is not implemented), Stage 6 (free floating), column-level and relaxation-zone restrictions of spec items 7 and 8 (not implemented).

## Stage results

| stage | status | numbers |
|---|---|---|
| 0 map | done | 00_map.md |
| 1 spec | done | 01_spec.md (continuity identity shown) |
| 2 implement | done for prescribed-motion path; column and zone restrictions and body velocity open | build 0 errors |
| 3 static invariance, staticFvMesh | PASS for waveFoamSWENSE tutorial and box (bit identical); waveIsoFoamSWENSE 4.4e-16 alpha, 5.5e-16 U, 2.3e-16 p | see below |
| 3 static invariance, zero motion | FAIL | wfs alpha 1.1e-12, U 3.7e-11 rel, pc 3.7e-11 abs; wis alpha 1.3e-11, U 7.5e-12 rel; box alpha 6.9e-13, Uc 1.8e-11 abs, pc 1.8e-8 abs (4.7e-7 rel), U 5.1e-9 rel |
| 4 moving pure waves | not run | - |
| 5 prescribed body | not run | - |
| 6 free floating | blocked (no own 6DoF library) | - |

Stage 3 field comparison at t = 0.118 (20 steps of deltaT 0.0059), max|diff| and relative to max|field|:

| case | variant | alpha.water | Uc | pc | U | p |
|---|---|---|---|---|---|---|
| wfs coarse | static vs prod | 0 | 0 | 0 | 0 | 0 |
| wfs coarse | zero vs prod | 1.1e-12 | 9.4e-14 abs | 3.7e-11 abs | 9.4e-14 (3.7e-11 rel) | 3.6e-11 |
| wis coarse | static vs prod | 4.4e-16 | 0 | 0 | 1.4e-18 | 1.1e-11 abs |
| wis coarse | zero vs prod | 1.3e-11 | 0 | 0 | 1.9e-14 | 3.6e-11 |
| box L1_sw | static vs prod | 0 | 0 | 0 | 0 | 0 |
| box L1_sw | zero vs prod | 6.9e-13 | 1.8e-11 abs | 1.8e-8 abs | 1.8e-11 abs | 1.8e-8 abs |

Wall times (s, serial): wfs prod 9, static 9, zero 13; wis prod 8, static 8, zero 11 to 12; box prod 8, static 8, zero 15.

## Commits on SWENSE_dynamic

Commits are local on branch SWENSE_dynamic, not pushed. See git log for hashes.

## Runs

- Stage 3 pass 1 (9 jobs): wfs, wis and box prod, static and zero. Five failed in 0 to 5 s with missing dictionary entries (psiProj, cAlpha, psiProjFinal, displacementLaplacian syntax, missing mesh); fixed in the scripts and rerun.
- Stage 3 pass 2 (9 jobs) after the psiProj final-iteration fix: all rc 0.
- Stage 3 retry (3 zero-motion jobs) after the moved-points guard: rc 0, numbers unchanged.
- Early failed invocations are in the scratchpad logs and are not counted as results.

## Open issues

- Zero-motion Stage 3 failure: confirm the round-off hypothesis (V0/V factor in the moved-mesh ddt versus the static path) by a diagnostic on the existing run fields, or by a deliberate test, before any change. Decision for Sithik on R3 and on the pressure solver tolerance for the pass criterion.
- Stale coarse tutorials (fvSolution, fvSchemes) and the missing modes file in the repo.
- Column-level skip for moving columns (spec item 7), relaxation-zone FatalError (spec item 8), swenseBodyVelocity moving wall velocity (spec item 6): not implemented.
- sixDoF free floating (Stage 6): choice of library needed.
- Grid2Grid setupPoints re-call on moved points: API call only, not verified against the Grid2Grid basis semantics ("once per mesh" comment in Grid2Grid.H).
