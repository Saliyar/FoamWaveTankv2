# Phase B: specification of the SWENSE solver for FoamWaveTankv2

Status: specification. Words and equations only. No code from the other code base is quoted. foamStar references are study file:line references (study copy, outside this repository). FoamWaveTankv2 references are real paths in this repository.

Notation: U total velocity, Uc complementary velocity, UInc incident velocity, p total pressure, pc complementary pressure, pInc incident pressure, pStar = (rho/rho1) pInc, alpha = water volume fraction, alphaInc incident water fraction, w relaxation weight, eta(x, y, t) free-surface elevation, z vertical coordinate, h = mean water level (seaLevel), g gravity magnitude.

## B1. What the study code computes

### B1.1 Fields
- Incident (from the incident-field update each step): UInc (cells and faces), pInc (cells and faces), alphaInc (cells), UfInc (faces). Study: createSwenseFields.H:24-51.
- Derived: pStar = (rho/rho1) pInc (study foamStarSwense.C:246).
- Read and written: Uc, pc (study createSwenseFields.H:58-103). Uc is a vector, pc a scalar. Both are written.
- Total fields, reconstructed each step: U = Uc + UInc, p = pc + pStar (study pSwenseEqn.H:96, 99). p_rgh = p - rho gh (study pSwenseEqn.H:100).
- Flux: phi = phiComp + phiInc, phiInc = UfInc . Sf (study updateSwenseField.H:5-6, pSwenseEqn.H:94).
- Mixture: alpha1 (water), alpha2 = 1 - alpha1, rho = alpha1 rho1 + (1 - alpha1) rho2 (study createSwenseFields.H:181-194). rhoPhi = interpolate(rho) phi (createSwenseFields.H:209-220).
- Gravity: gh = g . x - ghRef, ghf = g . xf - ghRef (study createSwenseFields.H:265-290). ghRef from a reference height hRef.

Read vs written (case files): U, p_rgh, alpha.water, pc, Uc (MUST_READ, study createSwenseFields.H:59-70, 75-86, 120-131). Written: U, Uc, pc, alpha.water, rho, p_rgh (AUTO_WRITE). p is NO_WRITE.

### B1.2 One time step, in order (study foamStarSwense.C:148-287)
1. Time control and time increment.
2. Clear the wave-generator cache.
3. Update the incident fields (relaxation zone update), then phiInc and phi (study updateSwenseField.H:1-6).
4. Mesh update (not used in the reference case: static mesh; mesh.update is still called but does nothing).
5. Alpha solve (study alphaEqn.H): explicit-in-time flux correction with compression; then relaxation alpha = (1 - w) alpha + w alphaInc (alphaEqn.H:64); then rho and mixture update (alphaEqn.H:66-74).
6. pStar = (rho/rho1) pInc (foamStarSwense.C:246).
7. Momentum predictor for Uc (USwenseEqn.H).
8. Pressure corrector loop (pSwenseEqn.H), with nCorrectors iterations; the last corrector reconstructs Uc and phiComp.
9. Turbulence: laminar in the reference case (constant/turbulenceProperties simulationType laminar), so no turbulence update.
10. Relaxation of the complementary velocity: Uc = (1 - w) Uc (foamStarSwense.C:267-269).
11. Write and advance.

Outer PIMPLE loop: nOuterCorrectors 8 in the reference case (fvSolution PIMPLE). The alpha, momentum and pressure steps are repeated per outer corrector (study foamStarSwense.C:181-266).

### B1.3 Equations, term by term
Momentum for Uc (study USwenseEqn.H:14-24):

  Euler ddt of (rho Uc) + div(rhoPhi, Uc) [Gauss linearUpwindV]
  + rho div(phiCompAbs, UInc) [Gauss linear]
  + div(devRhoReff(Uc)) [laminar: viscous term]
  = fvOptions - D Uc - rho refAcc

with D = rho alpha2^16 (1 + |U - <alpha1 U>|)^adCoeff, and adCoeff the airDamping keyword (default 0). The reference case does not set airDamping, so D = 0.

Predictor (study USwenseEqn.H:31-42): momentum equation solved with right-hand side -grad(pc) - grad(rho) pfInc / rho1, projected on faces (reconstruct).

Pressure (study pSwenseEqn.H):
- rAU = 1 / A(UEqn); HbyA = constrainHbyA(rAU H(UEqn), Uc, pc).
- phiHbyA = interpolate(HbyA) . Sf + interpolate(rAU rho) ddtUfCorr(Uc, Ucf).
- phiSource = -(snGrad(rho) / rho1) pfInc |Sf| rAUf (study pSwenseEqn.H:25-30); phiHbyA += phiSource.
- constrainPressure(pc, Uc, phiHbyA, rAUf).
- Solve laplacian(rAUf, pc) = div(phiHbyA), with a reference cell if pc needs one.
- Reconstruct: phiComp = phiHbyA - flux(pc); Uc = HbyA + rAU reconstruct((phiSource - flux(pc)) / rAUf).
- Update Ucf; phi = phiComp + phiInc; U = Uc + UInc; p = pc + pStar.

The optional term includeDivUInc (case value false) adds div(phiInc) to the pressure right-hand side. It is off in the reference case, and the port keeps that switch with the same default.

Alpha (study alphaEqnEul.H:1-77; alphaEqnCNEul.H for the Crank-Nicolson variant, which the reference case does not use: ddt default Euler):
- Euler implicit transport of alpha1 with upwind flux, then nAlphaCorr MULES corrections (reference: nAlphaCorr 2, cAlpha 1, MULESCorr yes) with the compression flux phir = phic nHatf.
- rhoPhi = phiAlpha (rho1 - rho2) + phi rho2.

### B1.4 Incident velocity (water and air)
The SWENSE incident update (study swenseRelaxationZones.C:192-327):

Classification: a point (cell centre or face centre) is water if z < eta(x, y, t), else air. Cells and faces are classified separately. Boundary patches are classified the same way with patch values.

Time: the incident fields use the simulation time t. The wave time is t - refTime. No ramp factor and no current are applied in the SWENSE incident path (the dispatch used by the relaxation zone takes only the alpha classification; study waveGenerator.C HEAD dispatch lines 629-699). This differs from the ramped, cutoff-limited velocity used by the waveVelocity boundary condition (section B1.8). The port keeps this difference.

Position: wave-frame coordinates are the rotated, shifted point, x_w = R (x - refPos - U0 t) with R from refDirection (study waveGeneratorI.H:100).

Water velocity: UInc = U_HOS(x_w, t_w), from the HOS / Grid2Grid record.

Air velocity: UInc = U_Air2Air(x_w, y_w, z - h), where h = seaLevel. The air solution is evaluated at z measured from the mean water level, not from the local free surface (study waveGenerator.C HEAD lines 674-684).

Boundary faces: the same classification with patch face centres.

Level set (not used in the reference case): levelSetInc = gDirection . x + eta(x, y, t).

### B1.5 Incident pressure (water and air)
- Water: pInc = rho_water (pd - g z), where pd is the dynamic pressure from the HOS / Grid2Grid record and g the gravity magnitude. The hydrostatic term is included, with z the vertical coordinate of the point (study hos.C:349). No ramp factor is applied to pd in this path.
- Air: pInc = pAir(x_w, y_w, z - h) from Air2Air (linear pressure model, configured in waveProperties Air2Air block).
- Both are used at cell centres and at face centres, with the same classification as the velocity.

### B1.6 Interface-pressure correction (relaxation zone, study swenseRelaxationZones.C:253-267, uncommitted in study)
On each internal face whose two neighbouring cells have different classification (one water, one air), the incident face pressure is replaced by the water-side pressure at the free surface above the face centre:

  pfInc(f) = pWater(x_f, y_f, eta(x_f, y_f, t)).

Purpose: the water pressure includes hydrostatics and the air pressure does not, so this removes the mismatch at the interface. The port implements this as a new behaviour in the incident-field class (section B2), not as a copy of the study hunk.

### B1.7 Transition band (z-eta blend, study waveGenerator.C uncommitted hunks)
Applies whenever Air2Air is configured. Band full width h_b = 0.0242 m. With dz = z - eta(x, y, t) and |dz| < h_b / 2:

  wW = clip(0.5 - dz / h_b, 0, 1)
  UInc = wW U_water + (1 - wW) U_air
  pInc = wW p_water + (1 - wW) p_air

Outside the band the plain classification of section B1.4 applies. Whether the reference run contained this band is open (see 01_mapping_report.md section 12); the port implements it because the study build used for comparison contains it.

### B1.8 Relaxation
- Weight w: the farfield scheme of the reference case (waveProperties relaxationName farfield, relaxationScheme farfieldList, farfieldPatchNames inlet and outlet, farfieldDistances 8.76 8.76 m, blendingDistances 8.76 8.76 m, spatialExponent default 3.5).
  - For each farfield patch p with distance d_p (distance to the patch, computed by a patch-wave distance) and blend width b_p:
    - main_p = 1 if d_p < d_farfield - b_p, else 0;
    - blend_p = 1 if d_p < d_farfield, times (1 - main_p);
    - scaled_p = blend_p (1 - (blend_p d_p - (d_farfield - b_p)) / b_p);
    - s = max over all patches of (main_p + scaled_p);
    - w = (exp(s^k) - 1) / (exp(1) - 1), with k = spatialExponent = 3.5.
  - Note: this is computed from the cell distance to the farfield patches. It is 1 in the farfield, 0 away from it, and blends in between.
- Relaxed quantities: Uc is multiplied by (1 - w) once per time step (after the PIMPLE loop). Alpha is relaxed toward alphaInc with weight w (before the rho update).
- Where w is applied: in the same time step as the incident update.

### B1.9 Boundary conditions in Generaltesting_order3 (study 0/ files)
- U: inlet and outlet waveVelocity; atmosphere pressureInletOutletVelocity (value 0); bottom slip; front and down empty.
- Uc: inlet, outlet and bottom fixedValue (zero); atmosphere pressureInletOutletVelocity; empty sides.
- alpha.water: inlet and outlet waveAlpha (nonuniform value list of 47 entries, overwritten at run time); bottom zeroGradient; atmosphere inletOutlet (value 0).
- pc: inlet and outlet fixedFluxPressure (inlet) and zeroGradient (outlet); bottom fixedFluxPressure; atmosphere totalPressure (p0 = 0).
- p_rgh: inlet, outlet and bottom fixedFluxPressure; atmosphere totalPressure.

waveVelocity (study waveVelocityFvPatchVectorField.C:132-201): for each boundary face, if the face is dry, U = wind velocity; if partly or fully wet, U = wave velocity at the centre of the wetted area of the face, and the mixed-type value fraction is 1. The wetted area is computed by the wetFace geometry (clipping the face at eta).

waveAlpha (study waveAlphaFvPatchScalarField.C:148-175): refValue = wetted area / face area, clipped to [0, 1].

### B1.10 Pressure variable and gravity
- The solved variable is pc. p is reconstructed as pc + pStar.
- constrainHbyA and constrainPressure: the same maths in OpenFOAM 5.0 and v2412 (only the API changed: updateCoeffs became updateSnGrad). Use the v2412 versions.
- Gravity from constant/g; hRef from constant (default 0); ghRef = g . (direction) hRef.

### B1.11 Case keys read (Generaltesting_order3)
- PIMPLE: momentumPredictor, nOuterCorrectors 8, nCorrectors 3, nNonOrthogonalCorrector 0, correctPhi no, moveMeshOuterCorrectors no, includeDivUInc false, fsiTol and fsiMaxIter (multi-body keys; inactive here), skipIter (default 1), airDamping (absent: 0).
- waveProperties: waveType HOS, Grid2GridInputFile, rampTime 0, refPosition, refDirection, refTime 0, zCutOff 0.5, extrapolateCubic on, layerThickness 2, layerdzValue 0.02, relaxationName farfield, relaxationScheme farfieldList, farfieldPatchNames, farfieldDistances, blendingDistances, Air2Air {Ha 2.0, order 3, pressureModel linear, verbose yes, cutoff auto, cutoffC 1, modesFile}.
- fvSchemes: ddt default Euler; div(rhoPhi,Uc) Gauss linearUpwindV; div(phi,alpha) Gauss vanLeer01; div(phirb,alpha) Gauss interfaceCompression; div(phiCompAbsolute,UInc) Gauss linear.
- fvSolution: alpha.water.* (nAlphaCorr 2, cAlpha 1, icAlpha 0, MULESCorr yes, nLimiterIter 5, alphaApplyPrevCorr no, smoothSolver symGaussSeidel tol 1e-10); pcorr.* present.
- OPEN (must be settled in Phase C): the solver dictionaries named pc and pcFinal (the names the study solver asks for) are not present in the reference fvSolution as far as Phase B could read; only pcorr.* is. Phase C will check how the name is resolved, and the first run will show it.

## B2. How each piece is built in FoamWaveTank

Legend: REUSE = use as is; EXTEND = change an existing class in a backward-compatible way; NEW = new class or file in the SWENSE library.

| Piece | FoamWaveTank class or file it uses or extends | New code needed |
|---|---|---|
| Solver skeleton (time loop, PIMPLE, alpha, pressure) | applications/solvers/waveFoam/waveFoam.C:98-181 and createFields.H: EXTEND as the starting point (copy the structure, not the file set) | NEW solver swenseWaveFoam: the Uc/pc split, the momentum predictor with the incident pressure gradient, the pressure equation with phiSource, and the reconstruction of U and p (B1.2-B1.3) |
| Alpha solve (MULES, compression, nAlphaCorr) | waveFoam alpha solve (alphaEqnSubCycle.H) and the MULES/cAlpha machinery of the two-phase mixture | REUSE if numerically identical to B1.3 alpha step; otherwise NEW alpha step in the solver |
| Relaxation zone base | src/waves2Foam/relaxationZone/relaxationZone.H:55-109 (base class, relaxationNames, correct(), relaxationWeightsMomentum) | EXTEND: a relaxation class that relaxes Uc toward zero and alpha toward alphaInc, with the same weight w (B1.8) |
| Farfield weight | src/waves2Foam/relaxationZone/relaxationScheme (spatial, spatialInterpolation, empty) and relaxationWeight (exponential, thirdOrderPolynomial, freePolynomial, CourantCorrection). The exponential form there is 1 - (exp(sigma^e) - 1)/(e - 1), which is not the farfield form of B1.8. | NEW scheme class swenseFarfield with the B1.8 weight; REUSE the relaxationZone plumbing |
| Incident field: water | src/waves2Foam/waveTheories/externalWaveForcing/Grid2Grid/Grid2Grid.{H,C}: eta, U, pExcess (Grid2Grid.H:153-167). pExcess multiplies by factor(time), a sine soft-start with Tsoft (Grid2Grid.C:120-128, 259). | EXTEND or wrapper: the SWENSE wrapper calls the Grid2Grid Fortran functions for velocity and dynamic pressure pd without the soft-start (require Tsoft = 0 in the case, else fatal), and adds the hydrostatic term to give pInc. No change to Grid2Grid itself. |
| Incident field: air | none (FoamWaveTank has no Air2Air) | COPY Air2Air (section B5) and wrap it |
| Incident dispatch (classification, blend, interface fix, cells/faces/patches) | none | NEW class swenseIncident (B1.4-B1.7) |
| Wetted cells/faces | src/waves2Foam/convexPolyhedral (localCell, localFace, convexPolyhedral) may clip cells at a surface | Candidate REUSE for wetted geometry; to be proven identical to the study wetted-cell behaviour in Phase C, otherwise NEW |
| Boundary condition waveVelocity (SWENSE) | src/waves2Foam/fvPatchFields/waveVelocityNew (TypeName waveVelocity, study-equivalent name clash) | NEW BC swenseWaveVelocity implementing B1.9 with wetted face area (new name, see B3) |
| Boundary condition waveAlpha (SWENSE) | src/waves2Foam/fvPatchFields/waveAlphaNew (TypeName waveAlpha) | NEW BC swenseWaveAlpha implementing B1.9 (new name) |
| Constrain pressure / HbyA | OpenFOAM v2412 cfdTools (constrainHbyA, constrainPressure) | REUSE, no port |
| Gravity and hRef | OpenFOAM (uniformDimensionedVectorField g, readGravitationalAcceleration) | REUSE |
| Function objects (probes, fieldMinMax, volume average, alpha sum) | OpenFOAM v2412 standard function objects | REUSE |
| Build | top-level Allwmake (lists solvers and libraries) | EXTEND: add the new library and solver to the build |

Points where the FoamWaveTank structure forces a difference from the study code (for the Phase B checkpoint):
1. The FoamWaveTank Grid2Grid pressure has a soft-start factor; the port must run with Tsoft = 0 or bypass it.
2. The FoamWaveTank exponential weight is a different function from the study farfield weight; the port needs its own class (not an edit of the existing class).
3. FoamWaveTank's waveVelocity and waveAlpha BCs use the FoamWaveTank wave theory with its own ramp and cutoff behaviour. The port's BCs take the study behaviour (B1.9) under new names.

## B3. Names
- Solver: swenseWaveFoam (default), in applications/solvers/swenseWaveFoam.
- Library: libswenseWaveTank (default), in src/swenseWaveTank.
- New classes (provisional, FoamWaveTank style): swenseIncident (incident water/air dispatch and interface correction), swenseRelaxation (relaxes Uc and alpha), swenseFarfieldScheme (farfield weight), swenseWaveVelocityFvPatchVectorField (BC type swenseWaveVelocity), swenseWaveAlphaFvPatchScalarField (BC type swenseWaveAlpha), swenseAir2Air (wrapper, see B5).
- Case translation: the case keys for the BCs change from waveVelocity / waveAlpha to swenseWaveVelocity / swenseWaveAlpha. Recorded in the Phase D translation table.
- Because the new BC types have their own names, a case can load waves2Foam and libswenseWaveTank in one process without a type clash.

## B4. Grid2Grid for v2412
- FoamWaveTankv2 has no libGrid2Grid.so in the v2412 user library directory (checked: the v2412 platform lib directory has no Grid2Grid library).
- The source is in ThirdParty/Grid2Grid (untracked in this repository; not edited). The FoamWaveTank build script for Grid2Grid expects libfftw3 static and the Grid2Grid CMake build.
- Plan (Phase C step 1): build Grid2Grid in a scratch copy with the v2412 toolchain (gfortran from the v2412 environment), with FFTW3 static, and install libGrid2Grid.so in the v2412 user library directory, so that the library links with -L$(FOAM_USER_LIBBIN) -lGrid2Grid. Check with ldd.
- Install location: $FOAM_USER_LIBBIN for the v2412 user, i.e. /zhome/19/7/191420/OpenFOAM/sitali-v2412/platforms/linux64GccDPInt32Opt/lib.
- The swense library Make/options will link -lGrid2Grid, and the solver will load it as a normal shared library.

## B5. Air2Air
Files to copy (the Air2Air set, section 0.3 of the task). Git authorship of each file is shown; all are recorded under the author name "Saliyar" in the foamStar history. Sithik to confirm that this author name is his before Phase C copies the files.

| File (study path, copied as is) | Git authors (study history) | Role |
|---|---|---|
| waves/waveTheory/Air2AirEnhanced.H and .C | Saliyar (6 commits) | Air2Air solver: order 1/2/3, modal coefficients, pressure model |
| waves/waveTheory/Air2AirWrapper.H and .C | Saliyar (2 commits) | Interface used by the wave generator (velocity, pressure, update) |
| waves/waveTheory/analyticAirCoeffs.H | Saliyar (1 commit) | Analytic order 1/2/3 coefficients (header only in the study working tree; .C deleted there, present in its history) |
| waves/waveTheory/hosSpectralCutoff.H and .C | Saliyar (1 commit) | Spectral cut-off of the HOS air coefficients |
| waves/waveTheory/cosSinTransform.H and .C | Saliyar (1 commit) | DCT-II / DST-II on midpoint grids |
| waves/waveTheory/HOSModesReader.H and .C | Saliyar (2 commits) | Reader for the modes file (constant/modes_HOS_SWENSE.dat) |

Not compiled in the study working tree (not in the study build): hosSpectralCutoff, cosSinTransform, HOSModesReader, analyticAirCoeffs.C. Phase C copies them only if the Air2Air code needs them; otherwise they are left out. Phase C decides after checking the include closure of Air2AirEnhanced.

Study-code dependencies to remove when copying: the study version header, the study waveTheory base class and the study macros (for example the version macros and the study build string). These are replaced by FoamWaveTank's own classes and by the v2412 API. The Air2Air code itself does not change behaviour.

Not copied (company files, behaviour re-implemented in B2): waveGenerator (dispatch, z-eta blend, interface), swenseRelaxationZones (incident update and interface correction), the HOS wave theory, the wet cells, the farfield scheme, the BCs, the solver.

## Open points carried to Checkpoint B and Phase C
1. Reference fvSolution pressure dictionary names (B1.11).
2. Whether the reference run contained the transition band (B1.7).
3. Confirm that git author "Saliyar" is Sithik, for the Air2Air copy (B5).
4. Grid2Grid soft-start must be off (B2, row Incident water).
5. Whether FoamWaveTank's wetted-geometry classes (convexPolyhedral) match the study wetted-cell behaviour (B2).
