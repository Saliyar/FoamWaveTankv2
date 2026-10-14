# SWENSE port verification checklist (filled)

Date: 2026-10-14
Reference: foamStarSWENSE HEAD (16cd4910 + z−η blend + interface fix), described in words and file:line only.
Port: FoamWaveTankv2, branch `SWENSE` (commits up to 21d65af).
Reviewer: AI-assisted review (Claude), checked against the code and the run outputs.

Verdict: **Needs revision.** The port is faithful for what it implements (binary dispatch, interface-pressure rule, SWENSE momentum/pressure equations, relaxation, boundary conditions, η against HOS). It is not yet a faithful port of the blended reference: the z−η blend is not implemented, the VoF alpha overshoot is 1.04, and the complementary velocity grows strongly.

---

## Physics Implementation

### Incident Field
- [x] Grid2Grid/Air2Air dispatch location identified: `src/swenseWaveTank/swenseIncident.C`, `update()` (cells ~l.216–229, boundary ~l.242–248, faces ~l.261–268, boundary faces ~l.280–286).
- [ ] z−η blending equation verified: **NOT IMPLEMENTED.** Dispatch is a hard switch on `z − eta < 0`. The blend (weight clamp(0.5 − (z−η)/0.024, 0, 1)) exists only in the uncommitted foamStar working tree (waveGenerator.C ~l.700–721 and ~l.809–821), and was not ported.
- [x] Pressure dispatch: binary (same hard switch as velocity; `swenseIncident.C` l.223–228, l.263–264). foamStar's pressure dispatch is binary as well (its blend also applies to pressure, which is not ported).
- [x] Hydrostatic pressure included in p_I: `waterP()` = ρ_w(pExcess − g·z_w), `swenseIncident.C` ~l.140–150. Air pressure from Air2Air, no hydrostatic term added (`airP()` ~l.169).

### Complementary Velocity Equation (u_C)
- [x] Continuity: div(u_C) is not tracked as a source. Binary dispatch has no blend source, so none is needed for the port as written. (Consistent with the port; differs from the blended reference.)
- [x] Momentum RHS terms (`src/swenseWaveTank/solver/swenseUcEqn.H`):
  - [x] Pressure gradient: −∇p_C/ρ, in the predictor as `−fvc::snGrad(pc)` (l.34), and in the corrector through `phiHbyA` (`swensePcEqn.H`).
  - [x] Interface pressure term: `phiSource = −∇ρ/ρ1 · pfInc · |Sf| · rAUf` (`swensePcEqn.H:25–30`), and the predictor term `−∇ρ·pfInc/ρ1` (`swenseUcEqn.H:35`).
  - [x] Advection: `fvm::div(rhoPhi, Uc)` (`swenseUcEqn.H:17`).
  - [x] Time derivative: `fvm::ddt(rho, Uc)` (`swenseUcEqn.H:16`).
  - [x] Viscous: `turbulence->divDevRhoReff(rho, Uc)` (`swenseUcEqn.H:20`), laminar in the reference case.
  - [x] Damping: not included (air damping is zero in the reference case; `adCoeff` absent).
- [ ] Missing blend terms: **not implemented**, because the blend is not implemented.
  - [ ] div(u_C) source from ∇a: not present.
  - [ ] Momentum source from (u_C·∇a)dU: not present.
  - [ ] Modified pressure term with ∇a: not present.

### VOF and Interface Handling
- [x] alpha_I transport: `alphaInc` is not advected. It is a target: `alpha1 = (1−w)·alpha1 + w·alphaInc` (`swenseAlphaRelax.H:5`). The advected field is `alpha1`, with the total SWENSE flux `phi = phiComp + phiInc`.
- [x] Max alpha: **VoF 1.040 (spike over the run; ~1.02 at 30 s). isoAdvection 1.000.** Target 1.00 — VoF does not meet it.
- [x] Sharp density step: ρ = limited α·ρ1 + (1−α)·ρ2 (`swenseAlphaRelax.H:6–12`). Same as foamStar's sharp step at α = 0.5 after MULES.
- [x] Relaxation zone alpha blending: `alpha1 = (1−w)·alpha1 + w·alphaInc` after the alpha step (`swenseAlphaRelax.H:5`).

### Boundary Conditions
- [x] Inlet/outlet: incident velocity BC applied: `swenseWaveVelocity` sets the wave velocity at the wetted face centre (wind on dry faces) (`src/swenseWaveTank/fvPatchFields/swenseWaveVelocity/…C`, `updateCoeffs`).
- [x] Free-surface pressure at the outlet: `pc` uses `fixedFluxPressure` / `totalPressure` as in the reference `0/` files; no change.
- [x] Walls: `slip` on the bottom, `empty` on the sides, as in the reference. No interface-specific wall condition.

---

## Numerical Validation

### Surface Elevation (η)
- [x] vs HOS, NRMSE (normalised by H/2 = 0.0725 m), 30 s, 5 gauges:
  - x = −8: VoF 0.007, iso 0.007
  - x = −1: VoF 0.057, iso 0.066
  - x = 0: VoF 0.058, iso 0.080
  - x = 1: VoF 0.048, iso 0.076
  - x = 8: VoF 0.006, iso 0.005
  Target < 0.1: **pass at all gauges for both solvers.**
- [x] Lag: ≤ 0.008 s, about one time step (0.0059 s). Target ≤ 1 step: **pass.**
- [x] 5-gauge comparison: **pass** (both solvers).

Note: an earlier phase-lag finding (~0.35 s) was a column-mapping error in the comparison script (gauges compared with the HOS point 1 m away); corrected, and no physical lag remains.

### Complementary Velocity (Uc)
- [x] max|Uc| at t = 10 s: VoF 0.0277 m/s; isoAdvection 0.0564 m/s.
- [x] max|Uc| at t = 30 s: VoF 0.231 m/s; isoAdvection 0.421 m/s.
  Baseline references in the checklist: unblended 0.218, blended 0.214. VoF is within 10% of the blended number but the port is binary, so this is not a like-for-like match; isoAdvection is not within 10%.
- [ ] Growth 0–30 s: VoF 734%, isoAdvection 647% (from 10 s to 30 s). Target < 23%: **does not meet.** Growth is consistent with the reference behaviour the checklist expects ("Uc growth still occurring"), but the size needs explaining before the new formulation.
- [x] Ripples present: yes (visible in the 30 s max|Uc| trace; the trace oscillates with the wave period).

### Conservation
- [x] Water volume oscillation (relative change of Σα): VoF −0.65% to +0.82%; isoAdvection −0.92% to +0.67%. Target < 0.7%: **VoF and isoAdvection do not meet it** (both are oscillations at the wave period, not a drift).
- [x] Continuity error: local ~4.7e-5 per step, global ~1.6e-8; cumulative −0.0034 at 30 s. **Stable, not growing.**

---

## Code Structure Equivalence

### Files/Classes
- [x] Incident field class(es): FWT `swenseIncident` (`src/swenseWaveTank/swenseIncident.{H,C}`). Air2Air: `Air2AirWrapper`, `Air2AirEnhanced`, `analyticAirCoeffs.H` (copied, Sithik's). foamStar equivalents: waveGenerator, Air2Air classes.
- [x] Complementary solver: FWT `swenseUcEqn.H` (momentum) and `swensePcEqn.H` (pressure), shared by `waveFoamSWENSE` and `waveIsoFoamSWENSE` (`src/swenseWaveTank/solver/`). foamStar equivalents: USwenseEqn.H, pSwenseEqn.H.
- [x] Relaxation zones: FWT `swenseRelaxation` (weight, `src/swenseWaveTank/swenseRelaxation.{H,C}`) and `swenseAlphaRelax.H`. foamStar: swenseRelaxationZones, relaxationSchemeFarfieldList.
- [x] VOF handling: **reused FWT/OpenFOAM v2412 MULES** (`alphaEqnSubCycle.H`) in `waveFoamSWENSE`; FWT's own isoAdvection alpha step (copied from `waveIsoFoam`) in `waveIsoFoamSWENSE`.

### Deviations from foamStar (documented)
1. **z−η velocity and pressure blend not implemented.** Dispatch is binary. Reason: the blend was added in an uncommitted foamStar working-tree change; it was never ported. Impact: the port is not directly comparable with the blended reference.
2. **Grid2Grid comes from FWT's `externalWaveForcing`**, with soft start off (`Tsoft = 0`). Frame and time offset follow the FWT convention (`referenceTheta = −angle(refDirection)`, `referenceTime = −refTime`).
3. **Relaxation weight** is the farfield formula (exponential blend, normalised), as in foamStar. It is not a smoothstep.
4. **Solver stack:** OpenFOAM v2412 `pimpleControl`, VoF MULES and isoAdvection instead of foamStar's MULES/alpha code; FWT `constrainHbyA` and `constrainPressure`.
5. **Multi-body, forward speed and air damping omitted** (zero in the reference case).
6. **System FFTW** (`/usr/lib64`) via RPATH, and **portable build flags** (`-march=x86-64-v3`).

---

## Open Questions / Issues
- [ ] Mismatches vs foamStar: (a) z−η blend missing; (b) VoF alpha overshoot 1.04 (1.02 at 30 s); (c) Uc growth 734% over 10–30 s (target < 23%); (d) water-volume oscillation slightly above 0.7% for both solvers.
- [ ] Solver instabilities: none fatal. Alpha overshoot grows over time in the VoF run; no NaNs, no failures in either 30 s run.
- [x] Uc growth still occurring: yes. Expected until the new formulation, per the checklist.
- [ ] Root cause of the Uc growth and of the VoF overshoot is not established. The relaxation step is applied correctly (weight and target check out), so the next checks are: (i) whether the relaxation undoes itself each outer iteration (the alpha solve restarts from `alpha.oldTime()`), and (ii) the MULES limiter behaviour in the relaxation band.

---

## Review Status

**Reviewer:** AI-assisted review (Claude), for Sithik
**Date:** 2026-10-14
**Verdict:** ☐ Pass (faithful port) | ☐ Pass with notes | ☑ Needs revision | ☐ Major issues

**Notes:**
The port matches the reference for the parts it implements: incident dispatch structure and hydrostatic pressure, the interface-pressure rule, the SWENSE momentum and pressure equations with phiSource, relaxation order, boundary conditions, and the surface elevation against HOS (NRMSE ≤ 0.08, lag ≤ 1 step). It does not yet reproduce the blended reference: the z−η blend is missing, the VoF alpha overshoot (1.04) and Uc growth (about 700% from 10 to 30 s) are above target, and the water-volume oscillation is slightly above 0.7%. Recommended before the new formulation: port the blend (with its source terms) so the comparison is like-for-like, and isolate the VoF overshoot.

---

## Next Steps
- [ ] Address the blend (port z−η dispatch and the three source terms), or record the decision not to.
- [ ] Isolate the VoF alpha overshoot (relaxation vs MULES limiter).
- [ ] Explain the Uc growth against the reference run before starting the new formulation.
- [ ] Run new-formulation prompt (when ready).
- [ ] Compare the 100 s run vs the foamStar baseline.
