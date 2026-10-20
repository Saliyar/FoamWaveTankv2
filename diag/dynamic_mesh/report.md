# Dynamic mesh for waveFoamSWENSE and waveIsoFoamSWENSE: report

Status: STOPPED after Stage 1 (spec). Stage 2 onward not started.

1. Branch SWENSE_dynamic at a078d9b, worktree wt_SWENSE_dynamic. No push, no merge.
2. Stage 0 map done: diag/dynamic_mesh/00_map.md (commit 0e0a267).
3. Stage 1 spec written: doc/swense_dynamic/01_spec.md (same branch, see commits below).
4. STOP reason: spec item conflicts with R3 (pure waves keep Uc exactly 0).
5. Under deforming motion (M2 morphing, body heave) cell volumes change, so
   the complementary flux must satisfy div(phiComp) = div(meshPhi) != 0 and Uc
   is not zero, even for pure waves. This is physics, not an incident error.
6. For rigid translation (M1) div(vMesh) = 0 per cell and R3 holds exactly.
7. Decision needed: (a) restate R3 as "Uc = 0 for rigid motion only", or
   (b) keep R3 and restrict moving tests to rigid and volume-preserving motion.
8. Environment: OpenFOAM loaded is v2506, the brief says v2412; Make/options say OFVERSION=2206.
9. Plain solvers: FWT tutorials waveFoam and waveIsoFoam use staticFvMesh. The only dynamic-mesh tutorial is FloaterWaveFoam (FloatStepper, out of scope).
10. Continuity: the brief says absolute flux; plain solvers use relative flux. The spec keeps relative.
11. Grid2Grid/Air2Air evaluation points are built once from t=0 geometry (swenseIncident.C:199). Must be rebuilt on a moving mesh.
12. Column level correction and relaxation weights are built once (swenseColumnLevel.H:49, swenseRelaxation.C:31). Spec restricts them to motion-free columns and zones.
13. gh and ghf are computed once (swenseSetup.H:219-220). Must be recomputed as plain waveFoam does (waveFoam.C:135-136).
14. The sixDoF force reads p and rho by name (sixDoFRigidBodyMotionSolver.C:238-248). SWENSE p and rho are the totals, so the names hold.
15. Nothing was built or run. No compute spent.

## Stage results

| stage | status | numbers |
|---|---|---|
| 0 map | done | 00_map.md |
| 1 spec | done, STOP on R3 | 01_spec.md |
| 2 implement | not started | - |
| 3 static invariance | not started | - |
| 4 moving mesh pure waves | not started | - |
| 5 prescribed body motion | not started | - |
| 6 free floating | not started | - |

## Commits on SWENSE_dynamic

| hash | files |
|---|---|
| 0e0a267 | diag/dynamic_mesh/00_map.md |
| (this commit) | doc/swense_dynamic/01_spec.md, diag/dynamic_mesh/report.md |

## Runs

None. Stage 0 and 1 are read-only and document-only.

## Open issues

- R3 conflict (see summary item 4 to 7): needs Sithik's decision before Stage 2.
- Build against OpenFOAM v2506 with -DOFVERSION=2206 in Make/options: check the build before Stage 2.
- Grid2Grid and Air2Air forward-only time buffer under repeated calls at the same time: not verified (spec item 9).
- sixDoF library key in dynamicMeshDict (motionSolverLibs) to be checked against the OpenFOAM v2506 tutorial.
- Not read in Stage 0: solver/swenseAlphaRelax.H, swenseIncident.C constructor body beyond the mapped lines, blendUP, and the Air2Air call sites.
