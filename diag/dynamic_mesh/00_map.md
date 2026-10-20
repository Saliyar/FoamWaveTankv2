# Stage 0 map: dynamic mesh in waveFoamSWENSE and waveIsoFoamSWENSE

Base: branch SWENSE_dynamic at a078d9b. No runs. All paths relative to the
worktree root. Line numbers from the files as they are at a078d9b.

Environment note: the loaded OpenFOAM is v2506 (WM_PROJECT_VERSION=v2506).
The Make/options of all four solvers say -DOFVERSION=2206 and do not name
OpenFOAM v2412. The brief's v2412 is not what is installed; this is an open
issue for the build stage.

## 0a. Plain solvers (waveFoam, waveIsoFoam)

Both plain solvers already call the OpenFOAM dynamic mesh API. Dynamic mesh is
selected only by constant/dynamicMeshDict. The tutorials in this repo set
`dynamicFvMesh staticFvMesh` (tutorials/waveFoam/constant/dynamicMeshDict,
tutorials/waveIsoFoam/constant/dynamicMeshDict). The only dynamic-mesh
tutorial in the repo is FloaterWaveFoam, which uses the FloatStepper
floaterMotion library and is out of scope.

| item | waveFoam | waveIsoFoam |
|---|---|---|
| createDynamicFvMesh | waveFoam.C:76 | waveIsoFoam.C (same include, before createFields) |
| moving-mesh controls (readDyMControls, createDyMControls) | OpenFOAM dynamicFvMesh headers, included by waveFoam.C:83, :100 | same headers |
| CourantNo / alphaCourantNo | waveFoam.C:91, :108-109 (inside `if (!LTS)` and per step, before ++runTime) | waveIsoFoam.C:103, :112-113 |
| mesh.update() placement | waveFoam.C:122-124: `if (pimple.firstIter() \|\| moveMeshOuterCorrectors) { mesh.update(); ... }` at the start of each outer corrector | waveIsoFoam.C:125 (same pattern) |
| topology change clear | waveFoam.C:130-133 (talphaPhi1Corr0.clear) | not used in waveIsoFoam |
| gh, ghf recompute | waveFoam.C:135-136 (`gh = g & (mesh.C() - referencePoint)`, `ghf` on Cf) | not checked in this pass |
| MRF update | waveFoam.C:138 | not checked |
| absolute flux from mesh velocity | waveFoam.C:144 `phi = mesh.Sf() & Uf();` | waveIsoFoam.C:156 same |
| correctPhi | waveFoam.C:140-146 (correctPhi.H, pEqn-side) | waveIsoFoam.C:152-158 |
| relative flux | waveFoam.C:149 `fvc::makeRelative(phi, U);` | waveIsoFoam.C:161 |
| mixture.correct after correctPhi | waveFoam.C:151 | not checked |
| checkMeshCourantNo | waveFoam.C:154-157 (meshCourantNo.H) | waveIsoFoam.C:166-168 |
| initial correctPhi | initCorrectPhi.H:3-19 (correctPhi.H once at start, with rAU) | waveIsoFoam.C: createUfIfPresent at :101 |
| pressure: makeRelative/makeAbsolute of phiHbyA | pEqn.H:19 (MRF), :23 `fvc::makeRelative(phiHbyA, U)`, :25 `fvc::makeAbsolute(phiHbyA, U)` | not checked |
| pEqn ddt correction | pEqn.H:17 `fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)` | not checked |
| momentum | UEqn.H:5 `fvm::ddt(rho, U) + fvm::div(rhoPhi, U)` | not checked |
| relative flux at the end of pEqn | pEqn.H:70 `fvc::makeRelative(phi, U)` | not checked |
| Uf | createUfIfPresent.H (OpenFOAM VoF helper) | same |

Continuity on the relative flux: the plain solvers keep phi relative to the
mesh (phi = absolute flux minus meshPhi, pEqn.H:70). This is the ALE
continuity. The revised brief's item 5 says "continuity is on the absolute
flux". The spec (Stage 1) resolves this by keeping the continuity on the
relative flux, as the plain solvers do.

## 0b. 6DoF and body-motion libraries used by the plain solvers

- The plain solvers link only -ldynamicFvMesh (Make/options, EXE_LIBS). They
  do not link a 6DoF library by name. The motion solver is chosen at run time
  by constant/dynamicMeshDict through the OpenFOAM dynamicFvMesh library.
- 6DoF solver: OpenFOAM v2506 src/sixDoFRigidBodyMotion (sixDoFRigidBodyMotion,
  sixDoFRigidBodyMotionSolver, sixDoFSolvers: Newmark, CrankNicolson,
  symplectic) and src/sixDoFRigidBodyState. Selected with `solver
  sixDoFRigidBodyMotion;` under a dynamicMotionSolverFvMesh (OpenFOAM
  overInterDyMFoam floatingBody tutorial is the reference dictionary). The key
  for loading the library (motionSolverLibs) must be checked against the
  tutorial when the test case is built (Stage 2).
- Force integration: sixDoFRigidBodyMotionSolver.C:189 (solve()). It is called
  inside mesh.update(). Force is `functionObjects::forces` with patches from the
  dict, `rhoInf` (= rhoInf_, default 1.0, read only when rho name is rhoInf,
  line 112-114) and `rho` = rhoName_ (default "rho", line 94). The forces object
  reads the pressure field named by its default pName "p" and the density field
  named by "rho". U is not an input of the sixDoF update.
  Lines: 238-248 build the forces dict and call `f.calcForcesMoments()`, then
  `motion_.update(firstIter, forceEff + mass*g, momentEff + ...)` (line 250+).
- Timing: forces are computed at the start of mesh.update(), that is at the
  top of each outer corrector (waveFoam.C:122-124 for plain waveFoam). They use
  the p and rho from the previous pressure solve. For plain waveFoam p is
  the solver's p field (created in createFields.H).
- Gravity for the body: g from the time registry (uniformDimensionedVectorField
  "g") or from the motion dict (line 171-185).
- Floater-specific libraries (floaterMotion, floaterMeshMotionSolver) are part
  of the FloatStepper submodule and are out of scope.

## 0c. SWENSE solvers: geometry reads and caches

Stage 0c is the table the Stage 2 changes are checked against. "Once" means
computed at construction or first step, and kept. "Every step" means
recomputed in each call.

| item | file:line | computed | moving-mesh action |
|---|---|---|---|
| gh, ghf (gravity potential at C and Cf) | solver/swenseSetup.H:219-220; used in solver/swensePcEqn.H:90 (p_rgh) | once | recompute on mesh.changing() (as waveFoam.C:135-136) |
| p_rgh = p - rho gh | solver/swenseSetup.H:229-234; swensePcEqn.H:90; swenseTimeLoop.H:105 | every step | follows gh, ok after recompute |
| phiInc initial = UfInc . Sf | solver/swenseSetup.H:104-109 | once | recompute with new Sf (only the start value) |
| phiInc per step = UfInc . Sf | solver/swenseTimeLoop.H:26 | every step | ok, uses the current Sf (mesh.update must come first) |
| phiInc projection psi, laplacian on current mesh | solver/swenseTimeLoop.H:29-57 | every step | ok, uses current mesh.magSf |
| phiComp (READ_IF_PRESENT, from fvc::flux(Uc)) | solver/swenseSetup.H:251-256 | once at start | restart: needs the relative/absolute convention decided in Stage 1 |
| Ucf (face Uc, normal part from phiComp) | solver/swenseSetup.H:271-279; rebuilt solver/swensePcEqn.H:79-83 | every pc solve | rebuild with the current Sf/magSf; phiComp convention must match |
| rhoPhi = interp(rho) phi | solver/swenseSetup.H:170-175, solver/swenseTimeLoop.H (via phi) | once, then updated where phi changes | must be relative (rhoPhi uses phi - meshPhi) |
| rhoPhiI = interp(rhoI) phiInc | solver/swenseSetup.H:260-265 | once | must be relative (phiInc - meshPhi) |
| swenseIncident pts_ (Grid2Grid/Air2Air evaluation points: all cell centres, then boundary face centres) | src/swenseWaveTank/swenseIncident.C:199-~260 (buildPointTables, uses mesh_.C() at :201, mesh_.Cf() at :202) | once (constructor) | MUST be recomputed after each mesh.update(), or the incident is evaluated at the t=0 positions. Index maps offCb_, offCfb_ rebuilt with it |
| swenseIncident per-point caches etaPt_, pExPt_, UPt_ (setSize pts_) | src/swenseWaveTank/swenseIncident.C:190-197 | once, filled per step | resize/refill follows pts_ |
| etaAt(ip,t) (fast path etaPt_ when etaTime_==t) | swenseIncident.C:370-378 | every call | depends on pts_ indices: must be refilled after the moving pts_ |
| weights(x, eta, alphaI, wWater, wAir) | swenseIncident.C:381-418 | every call | uses x from pts_: ok once pts_ is current |
| blendAlpha(t) | swenseIncident.C:485-519 (C at :495, Cf boundary at :503) | every step | uses current C and Cf, but etaAt(celli,t) indexes pts_: consistency needed |
| update(t, alphaI) | swenseIncident.C:521-~656 (C at :589-593, Cf at :631, alphaIf at :616) | every step | uses current C and Cf; must run after mesh.update() |
| psiWValue(ip,t) (psiWf used in phiSource0) | swenseIncident.C:658 | every step | index by pts_ |
| interfaceCellHeight(mesh) (delta thickness) | swenseIncident.C:26-53 | once (constructor) | Decision in Stage 2: recompute when the mesh moves, or document as static |
| blendThickness(db) | swenseIncident.C:54 | once | none |
| swenseRelaxation weight w (from patchWave distance to the wave patches) | src/swenseWaveTank/swenseRelaxation.C:31-58 (patchWave at :49, distance at :50) | once | choose restrict in Stage 1 item 8 (zero motion inside zone, FatalError otherwise) |
| relaxation of Uc: `Uc *= (1 - w)` | solver/swenseTimeLoop.H:100 | every step | uses w from above |
| swenseColumnLevel columns: cell bottoms/tops, dzCell, zBotCell, zTopCell, colOfCell, nbUp, nbDown, colX, colY, rowOfCol, cellsBelow | solver/swenseColumnLevel.H:49-75, 245-246, 345 (static bool colInit) | once (static) | restrict: eligible only if displacement since t=0 is below 1e-12 m (Stage 1 item 7). Log skipped columns each step |
| column skip box (levelCorrectionSkipMin/Max, by column x,y) | solver/swenseColumnLevel.H:98-124 | once | keep working; combine with the motion restriction |
| column level volume clip: sum(mag(a1Pre - a1Clip)*mesh.V()) | solver/swenseColumnLevel.H:371 | every step | uses current V: ok |
| swenseBodyVelocity: Uc = -UInc on body patch | src/swenseWaveTank/fvPatchFields/swenseBodyVelocity/swenseBodyVelocityFvPatchVectorField.C:76-80 | every evaluation | Stage 1 item 6: Uc = Uwall - UInc on moving body |
| swenseWaveAlpha BC: Cf (blend on boundary faces) | src/swenseWaveTank/fvPatchFields/swenseWaveAlpha/swenseWaveAlphaFvPatchScalarField.C:95 | every evaluation | ok after mesh.update() |
| swenseWaveVelocity BC: magSf (wet area check) | src/swenseWaveTank/fvPatchFields/swenseWaveVelocity/swenseWaveVelocityFvPatchVectorField.C:112-122 | every evaluation | ok |
| swenseVoFAdvection: phic = cAlpha mag(phi/magSf) | src/swenseWaveTank/swenseVoFAdvection.C:74 | every call | uses the flux passed in: must be the relative flux |
| swenseVoFAdvection: deltaN with average(V) (1e-8/pow(average(V),1/3)) | swenseVoFAdvection.C:104 | every call | ok; uses current V |
| swenseVoFAdvection: interface normal (gradAlphaf & Sf) | swenseVoFAdvection.C:141 | every call | ok |
| swenseVoFAdvection: alphaPhi10 = phi upwind interp (alpha) | swenseVoFAdvection.C:108-130 | every call | the flux is the relative flux |
| swenseIncident constructor: pInc_, UfInc_ fields | swenseIncident.C:112-124 | once (created), filled per update | refilled by update |
| Air2Air / Grid2Grid external wave: externalWave->step() | solver/swenseTimeLoop.H:15 | every step | Grid2Grid record update. Must check repeated calls at the same time (Stage 1 item 9) |
| swenseBudgetLog (logging only, budgetLog switch) | solver/swenseBudgetLog.H:29-38, 51-65, 90, 131-135 (mesh.C, points, V) | every step when on | logging only; uses current geometry; bins by x |

Not read in this pass (must be read before Stage 2): solver/swenseAlphaRelax.H
(the relaxation densities and its geometry use), swenseIncident.C:63-198 (the
constructor body beyond the lines above), swenseIncident.C:420-484
(blendUP), swenseIncident.C:521-656 (update body beyond the lines above), and
Air2Air calls.

## 0d. Time derivatives and convection terms

Notation: phi absolute total flux, phiI = phiInc absolute incident flux,
phiC = phiComp absolute complementary flux, meshPhi the swept mesh flux,
rhoPhi = interp(rho) phi, rhoPhiI = interp(rhoI) phiInc.

| term | file:line | flux used | moving-mesh form needed |
|---|---|---|---|
| alpha1 advection (advective form, div(phi) source) | src/swenseWaveTank/swenseVoFAdvection.C:37-44 (source divPhi from phi at :43); restart from oldTime at :46 | phi (absolute) | relative flux phi - meshPhi, plus the volume ratio V0/V on the old-time value (as OpenFOAM alphaEqn does with Vsc0/Vsc) |
| alphaI advection (same function) | solver/swenseTimeLoop.H:26, :55-56; alphaStep.H:6 (waveFoamSWENSE) | phiInc (absolute) | phiInc - meshPhi, same function, same V0/V factor as alpha1 |
| alpha1 and alphaI advection in isoAdvector | applications/waveIsoFoamSWENSE alphaStep.H:8-25 (Sp = fvc::div(phi), phiInc with SpI) | phi and phiInc | same as above, with the isoAdvector relative flux |
| rho convection and rho ddt in Uc: fvm::ddt(rho, Uc) | solver/swenseUcEqn.H:33 | Uc unknown, rho field | standard moving-mesh fvm::ddt (V0/V) |
| convection of Uc: fvm::div(rhoPhi, Uc) | solver/swenseUcEqn.H:34 | rhoPhi (absolute mass flux, from the rhoPhi definition at swenseSetup.H:170-175) | relative mass flux rhoPhi - interp(rho) meshPhi |
| mass source S = fvc::ddt(rho) + fvc::div(rhoPhi) (implicit -S Uc) | solver/swenseUcEqn.H:36 | rhoPhi absolute | relative flux and the same ddt as fvm::ddt(rho, Uc) |
| D term: D = rhoPhi - interp(rho/rhoI) rhoPhiI | solver/swenseUcEqn.H:17-22 | both absolute | relative: D_rel = D - meshPhi (rho_f - interp(rho/rhoI) rhoI_f). Not identically zero: the interpolation of the ratio differs from the product of interpolations. Keep the form, add the mesh term (Stage 1 item 4) |
| div(D, UInc) - (div D) UInc | solver/swenseUcEqn.H:25-29 | D (absolute) | D_rel |
| incident time derivatives | none. The incident fields (UInc, pInc, UfInc, alphaI target) are analytic/Grid2Grid/Air2Air functions of (x, t); grep found no fvc::ddt in src/swenseWaveTank/swenseIncident.C | n/a | the derivative of the incident field enters only through UInc in the D term and the time loop. Re-evaluate at the new positions (Stage 1 item 2-3) |
| phiSource0 = -psiWf snGrad(rho-rhoI) magSf | solver/swenseUcEqn.H:46-50 | geometric, no flux | recompute each step on the new mesh (snGrad uses new delta) |
| predictor pressure: fvc::snGrad(pc) magSf | solver/swenseUcEqn.H:60 | none | ok |
| pc equation: phiHbyA = interp(HbyA).Sf + interp(rAU rho) fvcDdtUfCorr(Uc, Ucf) | solver/swensePcEqn.H:7-17 | absolute | moving-mesh ddt correction (ddtUcFlow scheme) with absolute face velocity Ucf and the swept volume |
| pc continuity: fvm::laplacian(rAUf, pc) == fvc::div(phiHbyA) | solver/swensePcEqn.H:34-36 | phiHbyA (absolute) | relative continuity (plain waveFoam pEqn works on relative phi, pEqn.H:19-25, :70). See spec item 5 |
| phiComp = phiHbyA - pcEqn.flux() | solver/swensePcEqn.H:73 | absolute | consistent with relative continuity |
| phi == phiComp + phiInc | solver/swensePcEqn.H:85 | absolute | phiRel = phi - meshPhi, set in the same place |
| Uc reconstruct from face flux | solver/swensePcEqn.H:74 | absolute | ok |
| phiInc projection: laplacian(psi) = div(phiInc), phiInc -= snGrad(psi) magSf | solver/swenseTimeLoop.H:29-57 | absolute phiInc | fine on the new mesh (solved every step); its divergence must be the relative one if the continuity is relative (spec item 5) |
| restart: phiComp read (READ_IF_PRESENT) | solver/swenseSetup.H:251-256 | absolute | the restart reads the absolute flux; the convention must be documented in the file |
| restart: Ucf rebuilt from Uc and phiComp | solver/swenseSetup.H:271-279 | absolute face velocity | rebuild after mesh.update() with the current Sf |
| restart: alpha1, alphaI, rho.oldTime | solver/swenseSetup.H:124-145 | time directory | at a restart, the old-time value and the first mesh.update() must agree |

## Summary of the risks found

1. Grid2Grid/Air2Air evaluation points (pts_) are built once from the t=0
   cell centres and face centres (swenseIncident.C:199). This is the largest
   change: on a moving mesh the incident fields are evaluated at the wrong
   points unless pts_ is rebuilt after each mesh.update().
2. Column level correction (swenseColumnLevel.H:49-124) and the relaxation
   weight (swenseRelaxation.C:31-58) are built once from the t=0 geometry.
3. gh and ghf are computed once (swenseSetup.H:219-220). The plain solvers
   recompute them (waveFoam.C:135-136).
4. The D term and the mass source use the absolute mass flux; both must use
   the relative flux (swenseUcEqn.H:17-36).
5. Continuity: the brief says absolute flux; the plain solvers use relative
   flux. Stage 1 keeps the relative-flux continuity.
6. The sixDoF force uses the field named p and rho by default. In SWENSE p
   is rebuilt after each pc solve (swensePcEqn.H:89, swenseTimeLoop.H:104),
   and the forces dict reads rho by name (sixDoFRigidBodyMotionSolver.C:243).
   The SWENSE rho field is named rho (swenseSetup.H:139-144), so the default
   name works. p is rebuilt, so the force at the top of the next mesh.update()
   uses the last pc solve.
7. The environment is OpenFOAM v2506; the brief and Make/options name
   v2412 and OFVERSION=2206. Build check needed.
