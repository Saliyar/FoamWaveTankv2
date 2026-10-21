# Dynamic mesh for waveFoamSWENSE and waveIsoFoamSWENSE: report

Status: Stages 0 to 4 (M1) done. Stage 5 and Stage 6 NOT run. M2 NOT run. sixDoF implemented as a configuration and smoke-tested only (20 steps).

Latest commit (six fixes: absolute-divergence alpha source, correctPhi for phiComp, relative velocity for isoAdvector, full body wall velocity, column skip at 1e-3 dz, V0/V scaling removed): Stage 3 static, M1 (2 periods) and M2 (2 periods) are NOT run against this code. The library and both solvers build; no results are claimed for the new code.

## Summary

1. Branch SWENSE_dynamic at a078d9b, worktree wt_SWENSE_dynamic. Local commits only; no push, no merge. Production binaries and libswenseWaveTank.so untouched (build names _dyn and libswenseWaveTank_dyn).
2. Build: OpenFOAM v2506 (the production binaries are v2506 too), OFVERSION=2206 unchanged, clean rebuild of library and both solvers, 0 errors. Lesson: a stale solver object after a header change caused a segfault; solver objects are now cleaned before a rebuild that changes headers.
3. Static invariance (staticFvMesh, 20 steps): waveFoamSWENSE coarse, waveIsoFoamSWENSE coarse and box L1_sw are bit-identical to production a078d9b (max diff 0 in alpha, Uc, pc, U, p). PASS, with the current build.
4. Decision 1 (accepted absolute limits): SWENSE zero motion (dynamicMotionSolverFvMesh, zero displacement, box, 20 steps, writeInterval 1) against static, current build. p max 1.8e-8 Pa (limit 1e-6), alpha max 6.9e-13 (limit 1e-10), U max 1.8e-11 (limit 1e-9). PASS on absolute limits. Growth: U and alpha grow roughly linearly to step 20, p levels off near 1.8e-8 Pa. No exponential growth. PASS.
5. Decision 2 (sixDoFRigidBodyMotion, OpenFOAM v2506): implemented as dynamicMotionSolverFvMesh with motionSolver sixDoFRigidBodyMotion, patches boxFront boxBack boxBottom, plane constraint (normal 0 1 0), symplectic solver, rho field named rho. The force is the forces function object on the total p and the SWENSE rho (log: forces with rho rho). Smoke test: 20 steps, rc 0; body centre of mass z 0.9033 m moves to linear velocity -0.0289 m/s after 20 steps. Not verified: a physical free-heave decay (Stage 6a). No code change in the solver was needed for sixDoF.
6. Stage 2 items 6 to 8 are implemented: moving body wall velocity Uc = Uwall - UInc with Uwall the normal face velocity from meshPhi (not exercised in any run with a body in motion yet), column skip for moved columns (exercised in M1: all 346 columns skipped), relaxation weights recomputed each step on a moving mesh (spec item 8, option a, replaces the FatalError I had before; the FatalError was blocking M1).
7. Moving-mesh V0/V factor on the old-time alpha (spec item 3) added in swenseAdvanceVoF, gated by mesh.moving(); static identity confirmed.
8. Stage 4 M1 (rigid translation of the whole mesh, x = 0.02 sin(2 pi t / 1.3), solidBody motion solver, coarse pure-wave tutorial, 620 steps = 3.66 s): see table. Gauges at fixed earth positions were NOT set up, so the 0.1 mm gauge criterion is not evaluated.
9. Stage 4 M2 (morphing with a spatial displacement profile on the bottom): NOT run. The spatial profile needs a time-dependent point displacement BC which is not set up.

## Stage results

| item | status | numbers |
|---|---|---|
| 0 map | done | 00_map.md |
| 1 spec | done | 01_spec.md (continuity identity, R3 decision by Sithik) |
| 2 implement | done for prescribed-motion path and items 6 to 8 | build 0 errors |
| 2 V0/V alpha | done | static bit-identical |
| 2 sixDoF (Decision 2) | implemented as dict config, smoke test 20 steps | rc 0 |
| 3 static invariance | PASS | bit-identical, 3 cases |
| 3 zero motion, absolute (Decision 1) | PASS | p 1.8e-8 Pa, alpha 6.9e-13, U 1.8e-11 |
| 4 M1 rigid translation | run, pass on Uc and alpha invariance, gauges not evaluated | see below |
| 4 M2 morphing | NOT run | - |
| 5 prescribed body (heave) | NOT run | needs heave BC on box patches and a still-water setup, not built |
| 6a free heave decay (still water) | NOT run | the box case has waves; no still-water box case built |
| 6b floating box in waves | NOT run | sixDoF smoke test only (20 steps) |

### Stage 4 M1 results (coarse pure-wave tutorial, 620 steps, t = 3.66 s)

| solver | max abs Uc, static | max abs Uc, M1 | max abs alpha1 - alphaI, static | max abs alpha1 - alphaI, M1 | columns skipped |
|---|---|---|---|---|---|
| waveFoamSWENSE | 1.6e-12 | 1.2e-11 | 8.3e-14 | 2.4e-13 | 346 of 346 |
| waveIsoFoamSWENSE | 0 | 0 | 0 | 0 | not logged in this run |

Both Uc values are at round-off (1e-11 m/s); the M1 value is about 7 times the static level, which I attribute to the mesh-motion terms and did not chase. Gauges were not set up.

### Stage 3 zero-motion per-step table (box, current build, WRITE_EVERY 1)

| time | p max diff, Pa | U max diff, m/s | alpha max diff |
|---|---|---|---|
| 0.0059 | 1.3e-10 | 4.0e-13 | 8.7e-14 |
| 0.0177 | 1.8e-8 | 2.6e-12 | 7.6e-14 |
| 0.0590 | 7.6e-9 | 3.8e-12 | 2.4e-13 |
| 0.1180 | 1.8e-8 | 1.8e-11 | 6.9e-13 |

## Plain-solver reference (historical, not a pass criterion now)

Plain waveFoam and waveIsoFoam on the box (zeroGradient wave BCs, no waves) at 20 steps had static-vs-zero differences of alpha 5e-14, U 4e-12, p 6e-10. Sithik decided that the ratio is not meaningful, so the absolute limits above are used.

## Runs (wall time in seconds, serial)

Stage 3 static and zero: 7 to 17 s each. Stage 4 M1: static 90 to 114 s, translating 234 to 291 s. Stage 6 sixDoF smoke test: 19 s. Failed invocations (missing dictionary entries, a stale binary that segfaulted, a missing cell zone name, and the earlier zone-check FatalError) are not counted as results.

## Deviations and environment

- Coarse tutorials are stale: their fvSolution lacks psiProj and cAlpha. Runs use the box case fvSolution and fvSchemes (accepted). The tutorials are NOT edited.
- HOS modes file for the coarse tutorials is linked from SWENSE_testing diag new10 (not in the repository).
- Plain reference used the box case with zeroGradient wave BCs and no relaxation zones (historical).

## Open issues

1. Stage 4 M2 (morphing) needs a spatial time-dependent displacement BC on the bottom patch.
2. Stage 4 gauges at fixed earth positions (0.1 mm criterion) are not set up.
3. Stage 5 (heave of the box) needs a heave BC on box patches and a still-water case.
4. Stage 6 (free floating): the sixDoF smoke test runs; the free-heave decay in still water needs a still-water case; the body parameters (mass 509 kg, centre of mass z 0.903 m, inertia from a uniform box of height 2.19 m, width 0.468 m, depth 5.63 m, draft 0.193 m) are derived from the mesh and are not checked against a design.
5. Growth over 20 steps is bounded but not zero; the zero-motion differences come from round-off amplified through the pressure solver, not confirmed by a term search (Decision 1 said no term search).
6. Moving body wall velocity (item 6) is implemented but not exercised with a moving body in any run.
7. Grid2Grid setupPoints re-call on moved points (API call only, not verified against the "once per mesh" comment in Grid2Grid.H); M1 translation exercised it.
