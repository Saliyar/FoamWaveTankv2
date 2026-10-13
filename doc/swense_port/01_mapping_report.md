# Phase 1: study of the foamStar SWENSE solver (OpenFOAM 5.0) for the FoamWaveTankv2 port

Status: study report. Behaviour is described in words and equations. No source code from the other code base is included here. Notation: foamStar file:line references point to the study copy, which is outside this repository.

Scope note: the code that the reference case ran is a working tree of foamStar at commit 16cd4910 plus uncommitted changes. Its exact contents are kept outside the repository, in the notes folder, and are not reproduced here.

## 1. Versions

- The reference case `Generaltesting_order3` ran with a build of the study code whose uncommitted changes are: an interface-pressure correction on the relaxation zone, an experimental transition band between water and air velocities and pressures, and a reworked Air2Air solver.
- The working-tree Air2Air code (Sithik's own) is kept in the notes folder.
- Blend status: the transition band has no switch. It is active whenever Air2Air is configured, and `Generaltesting_order3` is configured with Air2Air (waveProperties, Air2Air block). The run started before the band existed in any installed library build, so its results may not contain it. This is recorded as an open point for Phase 3.

## 2. Fields and time step

Fields (study file: applications/solvers/foamStarSwense/createSwenseFields.H):
- Incident: UInc (cell and face), pInc (cell and face), alphaInc, from the relaxation zone update.
- pStar = pInc * alphaInc.
- Complementary: Uc (read, written), pc (read, written, reference cell fixed by PIMPLE dictionary), Ucf = interpolate(Uc).
- Total, reconstructed each step: U = Uc + UInc, p = pc + pStar, p_rgh = p - rho * gh.
- Flux: phi = phiComp + phiInc, phiInc = UfInc . Sf.
- Mixture: alpha1 (water), alpha2 = 1 - alpha1, rho = alpha1 rho1 + (1 - alpha1) rho2, rhoPhi = interpolate(rho) phi.
- Gravity: gh = g . x - ghRef, ghf on faces.

One time step (foamStarSwense.C, main loop):
1. Update relaxation zone (incident fields) and recompute phiInc, phi.
2. Optional mesh update (not used in the reference case).
3. Alpha transport: MULES-style explicit flux correction with compression, then alpha relaxed toward alphaInc with weight w, rho and mixture updated.
4. pStar = (rho / rho1) pInc.
5. Momentum predictor (USwenseEqn.H).
6. Pressure corrector loop (pSwenseEqn.H), with reconstruction of Uc and phiComp.
7. Relaxation of Uc: Uc is multiplied by (1 - w) after the PIMPLE loop (foamStarSwense.C:267-269).

## 3. Equations

Total: U = Uc + UInc, p = pc + pStar.

Momentum (USwenseEqn.H): for the complementary velocity Uc,

  d(rho Uc)/dt + div(rhoPhi Uc) + rho div(phiCompAbs, UInc) + div(devRhoReff(Uc)) = fvOptions - D Uc - rho refAcc

with the air damping D = rho alpha2^16 (1 + |U - <alpha1 U>|)^adCoeff, switched on only when adCoeff > 0. Predictor: the pressure gradient term on the right is -grad(pc) - grad(rho) pfInc / rho1 projected on faces.

Pressure (pSwenseEqn.H): with rAU = 1 / A(UEqn), HbyA = constrained(rAU H(UEqn)),
  phiHbyA = interpolate(HbyA) . Sf + interpolate(rAU rho) ddtCorr(Uc),
  phiSource = -snGrad(rho) / rho1 * pfInc * |Sf| * rAUf,
  phiHbyA += phiSource,
  laplacian(rAUf, pc) = div(phiHbyA)  (with reference cell),
  phiComp = phiHbyA - flux(pc),  Uc = HbyA + rAU reconstruct((phiSource - flux(pc)) / rAUf),
  phi = phiComp + phiInc, U = Uc + UInc, p = pc + pStar.

Alpha: Euler or Crank-Nicolson time discretisation with off-centring from the ddt scheme; upwind-type high-resolution flux with compression; MULES-style corrections (nAlphaCorr). rhoPhi = phiAlpha (rho1 - rho2) + phi rho2.

## 4. Incident field (water and air)

Relaxation zone update, cell and face values (foamStar swenseRelaxationZones.C):
- Classification: a point is water if z < eta(x, y) (the test is on the point's z versus the free surface elevation). Cells and faces are classified separately.
- Water: velocity and pressure from the HOS / Grid2Grid wave theory.
- Air: velocity and pressure from Air2Air, evaluated at z measured from the mean water level seaLevel (not from the local free surface).
- Boundary patches: classified the same way, with patch values.
- levelSetInc = gDirection . x + eta(x) (used only in level-set mode).
- alphaInc (cells, not level-set mode): wetted volume fraction of the cell from the wetCell class: 1 if fully wet, 0 if fully dry, otherwise wetVolume / cellVolume.

Water-air dispatch in the wave generator (waveGenerator): velocity and pressure are split by the alpha value passed in (alpha > 0.5 is water). In the transition band the study code blends the two sides (section 5).

## 5. Transition band (z-eta blend) and interface-pressure correction

These two items were uncommitted in the study code. They are described here as equations so the port can implement them new.

Transition band, velocity: with a band width h = 0.0242 m and dz = z - eta, for |dz| < h/2 the weight of water is
  wW = clip(0.5 - dz/h, 0, 1),
and the velocity is wW U_water + (1 - wW) U_air, where U_water is the HOS velocity at the point and U_air is the Air2Air velocity at (x, y, z - seaLevel).

Transition band, pressure: the same weight wW and the same band, with p = wW p_water + (1 - wW) p_air.

Interface-pressure correction (relaxation zone): on each face whose two adjoining cells have different classification (one water, one air), pfInc is replaced by the water-side pressure evaluated at (x_f, y_f, eta(x_f, y_f)) as a value on the free surface. Purpose: the HOS water pressure includes hydrostatics and the Air2Air pressure does not, so this removes the mismatch at the interface.

## 6. Water pressure definition

HOS water pressure (foamStar hos.C): p = rho_water (pd - g_mag z), where pd is the dynamic pressure from the HOS record. So the incident water pressure includes the hydrostatic term. This is the form to keep in the port (pInc = rho (pd - g z)).

## 7. Relaxation

- Weight w: from a relaxation scheme (farfieldList in the reference case). Clipped to [0, 1]. Applied to the zone cells; zone extents and the farfield distances come from the waveProperties relaxationZone block (farfieldDistances, blendingDistances, farfieldPatchNames inlet and outlet).
- Relaxed quantities: Uc *= (1 - w) after each time step; alpha = (1 - w) alpha + w alphaInc before rho and mixture are updated.
- Fully wet and fully dry cells: alphaInc is 1 or 0; the interface-pressure correction (section 5) acts on faces where classification changes.

## 8. Boundary conditions used by Generaltesting_order3

- U: inlet and outlet waveVelocity; atmosphere pressureInletOutletVelocity; bottom slip; front and down empty.
- Uc: inlet, outlet, bottom fixedValue (zero); atmosphere pressureInletOutletVelocity; empty sides.
- alpha.water: inlet and outlet waveAlpha (nonuniform value list); bottom zeroGradient; atmosphere inletOutlet.
- pc: inlet and outlet fixedFluxPressure / zeroGradient; bottom fixedFluxPressure; atmosphere totalPressure.
- p_rgh: inlet, outlet, bottom fixedFluxPressure; atmosphere totalPressure.

The SWENSE-specific BCs (waveVelocitySwense and similar) are not used by this case.

## 9. Pressure and gravity

- The solved variable is pc (complementary pressure). p is reconstructed as pc + pStar.
- constrainHbyA and constrainPressure behave the same in OpenFOAM 5.0 and v2412 (the API changed: updateCoeffs became updateSnGrad). They are used from OpenFOAM v2412 directly.
- Gravity: g from constant/g; ghRef from hRef; gh and ghf as above.

## 10. Case input keywords (Generaltesting_order3)

- PIMPLE (fvSolution): momentumPredictor yes, nOuterCorrectors 8, nCorrectors 3, nNonOrthogonalCorrector 0, correctPhi no, moveMeshOuterCorrectors no, includeDivUInc false, fsiTol, fsiMaxIter, nCorrPISO-like counts.
- waveProperties: waveType HOS, Grid2GridInputFile, rampTime 0, refPosition, refDirection, refTime, zCutOff 0.5, extrapolateCubic on, layerThickness 2, layerdzValue 0.02, relaxationName farfield, relaxationScheme farfieldList, farfieldPatchNames (inlet outlet), farfieldDistances (8.76 8.76), blendingDistances (8.76 8.76), Air2Air block (Ha 2.0, order 3, pressureModel linear, verbose, cutoff auto, cutoffC 1, modesFile constant/modes_HOS_SWENSE.dat).
- fvSchemes: ddt Euler; div(rhoPhi,Uc) Gauss linearUpwindV; div(phi,alpha) Gauss vanLeer01; div(phirb,alpha) Gauss interfaceCompression; div(phiCompAbsolute,UInc) Gauss linear.

## 11. Helper classes in the study code (classification for the port)

Core classes (study code): swenseRelaxationZones, the wave generator, the HOS wave theory, the Air2Air solver (wrapper, enhanced solver, analytic coefficients, spectral cut-off, cosSin transform, modes reader), the wetted-cell classes, the relaxation scheme, the relaxation shape, the wave velocity and wave alpha boundary conditions, the pimple control, and the adjustAlphaPhi helper.

Air2Air files are Sithik's own and may be copied. Everything else is studied and re-implemented.

## 12. Open points for later phases

- Whether the reference run included the transition band (section 5). Affects Phase 3 comparison.
- Which binary produced Generaltesting (foamStarSwense or foamStarSwenseCSV). Not needed for the port.
- Grid2Grid: the same source is available in FoamWaveTankv2 ThirdParty and must be built for v2412.
- Name overlap: waves2Foam registers the boundary condition names waveVelocity and waveAlpha. SWENSE BCs need new names if both libraries are loaded.
