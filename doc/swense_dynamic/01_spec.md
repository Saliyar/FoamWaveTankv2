# Spec: dynamic mesh for waveFoamSWENSE and waveIsoFoamSWENSE

Status: Stage 1. Based on diag/dynamic_mesh/00_map.md (commit 0e0a267).
Symbols: phi absolute total flux, phiI = phiInc absolute incident flux,
phiC = phiComp absolute complementary flux, meshPhi swept mesh flux,
phiRel = phi - meshPhi, V cell volume.

## 1. Fluxes

- phiInc and phiComp stay absolute (earth frame). phi = phiInc + phiComp.
- alpha1 advected with phiRel = phi - meshPhi.
- alphaI advected with phiIRel = phiInc - meshPhi, same function, scheme,
  limiter and settings as alpha1 (swenseAdvanceVoF, swenseVoFAdvection.C:25).
- In pure waves (phiComp = 0, alpha1 = alphaI) the two advections are identical
  under motion, so alpha1 = alphaI stays exact at round-off.

## 2. Incident fields at the moved geometry

After mesh.update() and before any equation, at the new time:
UInc, pInc, the alphaI target (blend weight), rhoI at the new cell centres,
UfInc and phiInc = UfInc . Sf at the new face centres and new Sf, and the
inlet and relaxation values at the new boundary face centres.
This needs the Grid2Grid and Air2Air evaluation points rebuilt from the current
C and Cf (map 0c, swenseIncident.C:199 buildPointTables). Without this the
incident is evaluated at the t = 0 points.

## 3. Time derivatives

The incident fields have no fvc::ddt in src/swenseWaveTank (map 0d). The
incident enters the momentum equation only through UInc in the D term and in
the time loop. So no incident time derivative needs a choice of (a) or (b). The
time derivatives present are on the unknowns:
- fvm::ddt(rho, Uc) and the mass source fvc::ddt(rho): standard moving-mesh
  form, with the same V0/V factor as the alpha advection.
- alpha advection: restart from the old-time value with the V0/V factor (as
  the standard OpenFOAM VoF alphaEqn does with Vsc0 and Vsc).
With meshPhi = 0 these reduce to the static forms, bit-identical.

## 4. Mixed terms and the pure-wave source

Uc equation (swenseUcEqn.H):
- The mass flux in fvm::div(rhoPhi, Uc) and in S = fvc::ddt(rho) + div(rhoPhi)
  is relative: rhoPhiRel = rhoPhi - interp(rho) meshPhi.
- D = rhoPhi - interp(rho/rhoI) rhoPhiI. With relative fluxes
  D_rel = D - meshPhi (rho_f - interp(rho/rhoI) rhoI_f).
  The mesh term is not identically zero, because interp(rho/rhoI) differs
  from interp(rho)/interp(rhoI). It is kept.
- div(D_rel, UInc) - div(D_rel) UInc uses D_rel.

Pure-wave source, discrete form: in pure waves alpha1 = alphaI, so rho = rhoI
cell by cell (same function, swenseRho). Then rho/rhoI = 1 exactly, and
rhoPhi = rhoPhiI when phi = phiI. So D = 0 identically, hence
D_rel = 0, div(D_rel, UInc) = 0, div(D_rel) UInc = 0, phiSource0 = 0 (it
contains snGrad(rho - rhoI) = 0). With Uc = 0 the Uc equation then has zero
source. This holds in the discrete form on a rigid translation (meshPhi with
zero divergence per cell). On a motion that changes cell volumes it holds only
if phiComp = 0 is consistent with the continuity below (see item 5).

## 5. Continuity and correctPhi (CONFLICT, see the STOP section)

The brief says continuity is on the absolute flux. The plain solvers use the
relative flux (waveFoam pEqn.H:19-25, :70). For an incompressible fluid in a
deforming cell, the relative continuity is sum over faces of phiRel = 0 per
cell, i.e. div(phi) = div(meshPhi) = dV/dt in absolute terms. The
spec keeps the relative form, which is what the plain solvers do.

After mesh motion (mesh.changing() and correctPhi on): rebuild phiInc at the
new faces, rebuild phiComp from Ucf and project with the pc correction so that
the discrete relative continuity holds. Projection of phiInc keeps the existing
psi treatment (swenseTimeLoop.H:29-57).

## 6. Body boundary

Uc = Uwall - UInc on the moving body patch, with Uwall the face velocity of the
mesh (movingWallVelocity convention, face motion from meshPhi). Relative total
flux through the body faces is zero: phiComp = meshPhi - phiInc there.
With a static mesh Uwall = 0 and this reduces to Uc = -UInc
(swenseBodyVelocityFvPatchVectorField.C:76-80), bit-identical.
alphaI treatment on body faces is unchanged.

## 7. Column level correction

The columns are built once from the t = 0 geometry (swenseColumnLevel.H:49-75).
On a moving mesh a column is eligible only if every cell in it has max point
displacement since t = 0 below 1e-12 m. Other columns are skipped as in the
existing levelCorrectionSkipMin/Max box. The number of skipped columns is
logged each step. The skip-box keyword keeps working.

## 8. Relaxation and farfield zones

Choice: restrict. The zone weights (swenseRelaxation.C:31-58) are computed
once. On a moving mesh a zone with max displacement above 1e-12 m is a
FatalError. Reason: the weights come from the wave-patch distance, which is
costly to recompute each step, and the wave zones should not move in the
planned runs.

## 9. Forces for body motion

- The sixDoF force (sixDoFRigidBodyMotionSolver.C:238-248) reads p and rho by
  name. SWENSE p = pc + pStar (swenseSetup.H:222-227, swensePcEqn.H:89),
  rebuilt after every pc solve. The rho field is the SWENSE rho.
  So the names hold the totals. The force is computed inside mesh.update(),
  so it sees the last pc solve, which is what the outer corrector needs.
- Incident fields are re-evaluated after every motion update, inside the outer
  correctors. Grid2Grid and Air2Air must be checked for a forward-only time
  buffer (externalWave->step() at swenseTimeLoop.H:15 is called once per time
  step). A repeated call at the same time must not advance the buffer. Not yet
  verified; Stage 2 item.

## 10. Scope

waveFoamSWENSE gets the same dynamic mesh features as waveFoam (dynamicFvMesh
selected by dynamicMeshDict, mesh.update in the outer corrector loop, correctPhi,
meshPhi, 6DoF through dynamicMotionSolverFvMesh). waveIsoFoamSWENSE gets the
same as waveIsoFoam. FloatStepper is out of scope.

## STOP: a spec item conflicts with R3

R3 says: pure waves must keep Uc exactly 0.

Under mesh motion this cannot hold in general. Take M2 (displacementLaplacian,
bottom patch sliding tangentially with d_x = A(x) sin(omega t), A = 0.5 dx
sin(pi (x - x0)/(x1 - x0))). The mesh velocity has divergence
div(vMesh) = dA/dx omega cos(omega t) = 0.5 dx (pi/L) omega cos(pi (x - x0)/L)
cos(omega t), with L = x1 - x0, which is nonzero except on the nodes of the
cosine. Interior cell volumes therefore change. Relative continuity needs
sum phiRel = 0 per cell, so with phiInc divergence-free (pure waves) the
complementary flux must satisfy div(phiComp) = div(meshPhi) != 0, which gives
Uc != 0 at the order of dV/dt / V times the relative velocity. This is physics,
not an incident-model error, and it is not round-off. The same holds for the
body heave in Stage 5 (the box moves and cells change volume).

For a rigid translation (M1) div(vMesh) = 0 per cell, and R3 can hold exactly.

Decision needed (STOP here per the brief):
(a) restate R3 as: Uc = 0 for rigid motion only; under deforming motion Uc is
    driven by div(meshPhi) and is reported, not required to vanish;
(b) keep R3 as written and restrict the moving-mesh tests to rigid translation
    and volume-preserving motion (div vMesh = 0 per cell).

Stages 2 to 6 are not started. The implementation depends on this choice (it
changes the pass criteria of Stage 4 and of the continuity projection).
