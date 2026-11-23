# Body-zone fix (plainBodyZone, lam = 0): result 23 Nov 2026

Branch `swense-zone-fix` from 30df70a (not pushed).
Commits: 84d1cac (edge modes min/max/ramp + zone source-split log),
1eba198 (candidate c, `plainBodyZoneFace solenoidal`). All new switches are
off by default, so the defaults behave as 30df70a.

**Outcome: no candidate passes the gate. Stopped after step 2; steps 3 and 4 not run.**

## Set-up

- Coarse empty tank, blendThickness 0.0438666667, rampTsoft 12, 0-23 s.
- Reference F3b: the same with no zone. F3 used blend 0.0433, so it was rerun with 0.0438667.
- Runs are in `/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix/`
  (F3b, Z0, Za, Zb, Zc).
- Gate, on gauges x = ±1, ±1.5, ±2, last 2 periods before 23 s; Uc taken in water, 0.5 < |x| < 1.5, outside the zone, 12-23 s:
  - water Uc outside the zone ≤ 2× F3b;
  - wave error vs Stokes II ≤ F3b + 1 mm;
  - wiggles (harmonics n ≥ 4) ≤ 0.3 mm.

## Step 1: Z0, where the Uc outside the zone comes from

- **First crossing of 0.05 m/s:** at t = 0.077 s, in still water (the 12 s ramp has barely started).
  - The cell is 487, at x = −0.525 m and z = −0.013 m: the free-surface row, layer 1 outside the zone.
  - Its zone-side face is 970 (x = −0.5), with lamf 0 and faceLim 0.011.
- **Source split** (N/m³):

| Where | Cell closure (1−λ) rI | Face term (1−λf) G_f | Rest |
|---|--:|--:|--:|
| Around cell 487 (max) | 23.7 | **48.1** | 0.83 |
| Edge band (∫\|.\|dV) | 0.21 | **1.54** | 0.0015 |
| Rest of zone (∫\|.\|dV) | 1.73 | **13.6** | 0.0005 |

- **Dominant term: the face term G_f.**
  - In every zone column it is a vertical pair, −810 / +830 N/m³, in the rows just below and above the surface cell.
  - It is the hydrostatic mismatch between the analytic pInc and the smeared rhoI.
  - pc balances it inside the zone, but the columns outside have none. So pc jumps sideways at x = ±0.5 and drives Uc there.

## Step 2: candidates on Z0

| Run | Edge / face option | Water Uc outside, max / median (m/s) | Wave error (mm) | Wiggles (mm) | Gate |
|---|---|--:|--:|--:|---|
| F3b | no zone | 0 / 0 | 3.71 | 0.20 | reference |
| Z0 | min (as 30df70a) | 0.67 / 0.35 | 5.58 | 2.09 | fail |
| Za | a: lamf = max on edge faces | 0.59 / 0.34 | 6.02 | 1.84 | fail |
| Zb | b: λ ramp over 3 cells | 0.51 / 0.35 | 5.67 | 1.94 | fail |
| Zc | c: gradient part of G_f removed | 0.73 / 0.29 | **11.1** | 1.94 | fail |

- **Uc gate.** F3b's Uc is exactly 0 (an empty SWENSE tank forms no complement), so "≤ 2× F3" means Uc = 0. No zone with λ = 0 can meet that.
- **Candidate a** changes almost nothing: at the crossing Uc is 0.05674 against 0.05677. The edge faces carry only 0.011.
- **Candidate b** lowers the Uc peak by 25 %. Wave error and wiggles stay unchanged.
- **Candidate c** fixes the still-water source but makes the waves worse:
  - Uc outside at 1 s is 0.014 m/s, against Z0's 0.376.
  - By 6-8 s, as the wave builds, it is back to 0.25 m/s.
  - Gauges at ±1 and ±1.5 are 13-20 mm off.
  - Max air Uc is 27 m/s, against 4.8 in Z0, and dt is 14× smaller (run wall time 65 min, against 5).
  - The remaining source is the cell closure, which is wave-driven.

## Water Uc by x band (max, m/s)

| t (s) | Z0, zone | Z0, 0.5-1 | Z0, 1-1.5 | Zc, zone | Zc, 0.5-1 | Zc, 1-1.5 |
|--:|--:|--:|--:|--:|--:|--:|
| 1 | 0.26 | 0.38 | 0.02 | 0.02 | 0.01 | 0.01 |
| 4 | 0.13 | 0.28 | 0.16 | 0.14 | 0.16 | 0.06 |
| 8 | 0.18 | 0.33 | 0.20 | 0.24 | 0.24 | 0.26 |
| 16 | 0.20 | 0.35 | 0.28 | 0.27 | 0.19 | 0.19 |
| 23 | 0.22 | 0.39 | 0.36 | 0.28 | 0.16 | 0.31 |

## What this means

- Any zone with λ = 0 that cuts the free surface injects the incident residual's discrete error (hydrostatic and wave-driven) into Uc.
- Moving or smoothing the edge does not remove it.
- This agrees with phaseC.md C6 loop 4: a smaller zone gave no gain.
- A fix would have to make the incident residual zero for the discrete incident wave itself, not reshape the zone.

# Discrete incident hydrostatics (`incidentHydrostatic discrete`), 23 Nov 2026

Commit d74387f (switch, default `analytic` = unchanged). Static mesh only; serial or whole-column decomposition.

## Step 0: code check

- **(a) Part V is not the same thing.**
  - Part V solves `laplacian(rAUf, deltaHeld) = div(Phi[r_I])` over the whole r_I flux (cell part plus rAUf G_f) and applies it one step late ([swenseBalanceFlux.H:30-35](../src/swenseWaveTank/solver/swenseBalanceFlux.H#L30-L35), [swenseBalanceSolve.H:39-46](../src/swenseWaveTank/solver/swenseBalanceSolve.H#L39-L46)).
  - That is a global Poisson projection. It does not rebuild pInc column by column, and it leaves the solenoidal part. So I continued.
- **(b) The cell closure has no hydrostatic part.**
  - It is `(1-lam)*(ddt(rhoI UI) + div(rhoPhiI UI) - visc)` ([swenseUcEqn.H:192-196](../src/swenseWaveTank/solver/swenseUcEqn.H#L192-L196)).
  - Pressure and hydrostatics enter only through the face term faceLim. That term is reconstructed in the Uc predictor ([swenseUcEqn.H](../src/swenseWaveTank/solver/swenseUcEqn.H), `fvc::reconstruct(... - faceLim)`) and enters pc as a flux ([swensePcEqn.H:29-33](../src/swenseWaveTank/solver/swensePcEqn.H#L29-L33)).
  - It is not fvc::grad, so the "reconstruct when on" change was not needed.

## Step 1: the switch

- **pDyn** = p_side − rho_side g (eta − z).
  - p_side is the water model pressure where z < eta_HOS, and the Air2Air pressure above.
  - So pDyn = 0 in still water.
- **pHydI_rgh** is marched down each column: `P(lower) = P(upper) - ghf_f (rhoI(lower) - rhoI(upper))`.
  - It uses the same ghf as G_f.
  - The top cell is set to the analytic rho_side g (eta − z) − rhoI gh.
  - Boundary faces use the same relation from the owner cell.
- **pInc** = pHydI_rgh + rhoI gh + pDyn.
- This pInc is used in G_f, in the closure face term and in p. Files: [swenseHydroDiscrete.H](../src/swenseWaveTank/solver/swenseHydroDiscrete.H) and `swenseIncident::sidePressure`.

## Step 2: gate S (still water, coarse, zone on, 2 s)

The wave is scaled to about 3e-12 with rampTsoft 1e12, because Air2Air refuses H = 0.

| Run | max \|G_f\| (N/m³) | max cell closure | max \|Uc\| (m/s) | max \|pc\| (Pa) |
|---|--:|--:|--:|--:|
| S0, switch off | ±800 pair (as Z0) | ~0 | 0.35 | (as Z0) |
| S1, discrete | 3.3e-9 | 3.1e-9 | 1.5e-12 | 7e-10 |

**Gate S passes.** Everything is at round-off. The 3e-9 left is the 3e-12-scaled wave itself.

## Step 3: gate W (waves, zone on, 23 s; gauges ±1, ±1.5, ±2)

| Run | Water Uc outside, max / median (m/s) | Wave error (mm) | Wiggles (mm) | Gate |
|---|--:|--:|--:|---|
| F3b, no zone | 0 / 0 | 3.71 | 0.20 | reference |
| Z0, zone, analytic | 0.67 / 0.35 | 5.58 | 2.09 | fail |
| ZD, zone, discrete | 0.35 / 0.11 | 4.93 | 0.83 | fail |

- **Gate W fails.** Wave error is F3b + 1.2 mm, against the 1 mm limit. Wiggles are 0.83 mm, against 0.3. Uc is not zero.
- **ZD improves on Z0.** Uc outside the zone is halved (3× lower median), wiggles drop by 60 % and the error is 0.65 mm better.
- **First crossing of 0.05 m/s:** at 2.67 s, against 0.077 s in Z0. It is now wave-driven, not hydrostatic.
- **Source split, ZD, 12-23 s** (median N/m³):

| Where | Cell closure | Face term | Rest |
|---|--:|--:|--:|
| Around the max-Uc cell | 80 | **278** | 126 |
| Edge band (∫) | 5.9 | **9.4** | 0.26 |
| Zone (∫) | 54 | **60** | 0.97 |

- **What is left.**
  - The face term is still the largest, but the cell closure is now close to it.
  - In waves G_f reaches 900-1100 N/m³ on both vertical and horizontal faces.
  - That is the dynamic and lateral residual of the smeared wave interface, not hydrostatics.

## Step 4: F3b with the switch on

Not required, because W did not come within 1 mm. It ran in parallel anyway: F3bD is identical to F3b (Uc = 0; 3.706 mm and 0.195 mm). The switch does not touch a zone-free tank.

Runs (`zonefix/`): S0, S1, ZD, F3bD.

# Width, moving mesh, heave decay (switch on), 23 Nov 2026

Commit bdcffce: `incidentHydrostatic discrete` on a moving mesh.
- Columns are found once at set-up through the top and bottom faces.
- Each step the march uses the face's current gh_f and the moved incident points.
- Columns under the box wall start from the analytic value. The box mesh has 104 column starts: 100 plus 4 under the box.
- The static-mesh FatalError is removed.

**Outcome:** M0 is not round-off at 23 s, but the cause is the moving-mesh path, not the switch. Step 3 fails (no width within 0.5 %). Stopped before step 4.

## Step 1: zone width (empty tank, coarse, switch on, 23 s; gauges ±1, ±1.5, ±2)

| Run | Zone | Water Uc outside, max / median (m/s) | Wave error (mm) | Wiggles (mm) |
|---|---|--:|--:|--:|
| F3b | none | 0 / 0 | 3.71 | 0.20 |
| ZD | 1.0 m (\|x\| < 0.5) | 0.35 / 0.11 | 4.93 | 0.83 |
| W07 | 0.7 m (\|x\| < 0.35) = box motion range + 2 cells | 0.24 / 0.08 | 4.51 | 0.71 |
| W05 | 0.5 m (\|x\| < 0.25) | 0.21 / 0.07 | **4.22** | 0.72 |

- The motion range comes from test 2: the box spans x = −0.214 to 0.240 m, plus 2 cells of 0.05 m.
- **The error falls as the zone shrinks.** W05 is within 1 mm of F3b (F3b + 0.52 mm).
- Wiggles stay at about 0.7 mm, against the 0.3 mm limit, and Uc is not zero.

## Step 2: gate M0 (solidBody linearMotion with velocity 0, so mesh.moving() is true)

| Pair | Wave error (mm) | Largest gauge difference at 23 s (mm) | First Uc difference |
|---|--:|--:|---|
| M0 against ZD (switch on) | 4.86 / 4.93 | 0.20 | step 3, 2e-11 relative |
| M0off against Z0 (switch off) | 5.39 / 5.58 | 1.1 | step 3, 5e-10 relative |

- G_f (maxG) is bit-identical between M0 and ZD through step 50. The switch's own terms do not change on the moving path.
- The moving-mesh path alone (switch off) departs from static at about 1e-10 and grows to about 1 mm by 23 s.
- With the switch on it departs less (0.2 mm). Strictly, M0 is not round-off at 23 s, and this predates the switch.

## Step 3: still-water heave decay (moored box, test 2 settings, 0.15 m/s heave kick at t = 0, 6 s)

Fit of z = z0 + A e^(−ζωt) cos(ω_d t + φ) over 0.2-6 s:

| Run | Equilibrium z (mm) | Period (s) | ζ | A (mm) | Against plain fine: T / ζ / A |
|---|--:|--:|--:|--:|---|
| Plain fine (HPF) | −16.5 | 0.798 | 0.101 | 7.13 | reference |
| Plain medium (HPM) | −16.6 | 0.799 | 0.095 | 7.07 | +0.2 % / −6 % / −1 % |
| Plain coarse (HP) | −8.4 | 0.806 | 0.028 | 5.55 | +1.0 % / −72 % / −22 % |
| SWENSE zone 1.0 m (H10) | −16.7 | 0.808 | 0.050 | 5.86 | +1.3 % / −50 % / −18 % |
| SWENSE zone 0.7 m (H07) | −16.8 | 0.799 | 0.065 | 6.42 | +0.1 % / −36 % / −10 % |
| SWENSE zone 0.5 m (H05) | −16.8 | 0.799 | 0.093 | 7.54 | +0.1 % / **−7.5 %** / +6 % |
| SWENSE zone 1.0 m, switch off (H10off) | −17.2 | 1.15 | (poor fit) | | irregular, beating |

- **Step 3 fails.** No width is within 0.5 % in damping, whether against plain coarse or plain fine.
- **Plain coarse is a poor reference.** Its setFields start puts the water about 9 mm high, so the box floats 8 mm high, and its damping is a quarter of plain fine's.
- **Smaller zones come closer to plain fine.** H05 has the period within 0.1 % and the damping within 7.5 %. Plain medium is itself 6 % off fine.
- **The switch matters for the moving box.** Without it (H10off) the decay is irregular, with a stretched second cycle.

Runs (`zonefix/`): W05, W07, M0, M0off, heave/{HP, HPM, HPF, H10, H07, H05, H10off}.

# Step 4: moored box, zone 0.5 m, switch on (23 Nov 2026)

**Set-up.** Run B05 (`zonefix/B05`):
- coarse mesh, test 2 settings (blend 0.0438667, rampTsoft 12), 23 s;
- bodyZone |x| < 0.25, z < 0.1118;
- `incidentHydrostatic discrete` and `logZoneEdge` on (commit 47b76d2).

Plain fine and plain medium already existed (`meshes/`).

**Zone edge.** The smallest distance from the box wall to the zone edge is exactly 1 cell all run.
- It never drops below 0.99 in 4850 steps. The nearest edge is the zone top face, one row (0.049 m) above the box top.
- The zone cells inside innerDistance move rigidly with the box, so the gap does not change.

## Motion amplitude, last 5 periods (13.0-23 s): error against plain fine | against experiment

| Run | Surge | Heave | Pitch | Tension |
|---|--:|--:|--:|--:|
| Plain fine | 0 \| +11.0 % | 0 \| +6.7 % | 0 \| +8.1 % | 0 \| −10.8 % |
| Plain medium | −0.2 \| +10.8 % | −2.1 \| +4.5 % | +1.3 \| +9.5 % | −0.6 \| −11.4 % |
| Plain coarse | −3.0 \| +7.7 % | −7.2 \| −1.0 % | −39.7 \| −34.8 % | +5.0 \| −6.4 % |
| SWENSE test 2 (zone 1.0 m, switch off) | +1.6 \| +12.8 % | +4.8 \| +11.8 % | −19.0 \| −12.4 % | +14.0 \| +1.7 % |
| SWENSE B05 (zone 0.5 m, switch on) | +7.2 \| +19.0 % | +4.2 \| +11.2 % | **−62.3 \| −59.3 %** | **+32.4 \| +18.1 %** |

## Wave elevation against plain fine, rms (mm), same window

Elevation is the column integral of alpha from the written fields.

| Run | x = −1 | x = −0.25 | x = +0.25 | x = +1 | Mean at ±1 | Mean at the box (±0.25) |
|---|--:|--:|--:|--:|--:|--:|
| Plain medium | 1.5 | 6.4 | 5.4 | 2.3 | 1.9 | 5.9 |
| Plain coarse | 9.4 | 8.0 | 11.5 | 9.8 | 9.6 | 9.7 |
| SWENSE test 2 | 8.8 | 12.7 | 10.8 | 10.0 | 9.4 | 11.7 |
| SWENSE B05 | 7.6 | 23.6 | 19.6 | 8.7 | 8.2 | **21.6** |

**Result.** B05 is worse than test 2 near the box:
- pitch −62 % against plain fine;
- tension +32 %;
- wave at the box 21.6 mm against 11.7 mm.

It is only slightly better away from the box (8.2 mm against 9.4 mm at ±1). The still-water heave decay (step 3) does not carry over to the wave case: the small zone hurts pitch badly.

# Terms of rI in the band, from written fields (no solver run), 23 Nov 2026

**Method.** Script [diag/rI_terms.py](rI_terms.py) (Python on the written fields).
- Cases: ZD (empty tank, zone 1.0 m) and B05 (box, zone 0.5 m), both with the switch on.
- Times: 15, 19 and 22.9 s. 22.9 s instead of 23 s because the ddt central difference needs t ± 0.1 s.
- Band: 0.01 < alphaI < 0.99. Region: ZD outside the zone at 0.5 < |x| < 1.5; B05 within 0.3 m of the box faces.

**Terms** (per unit volume; B05 includes the moving mesh):
- **ddt**: (rhoI UI V) central difference over ±0.1 s, divided by V.
- **conv**: Gauss linear with flux interpolate(rhoI)·(UI_f·Sf − meshPhi). This approximates the solver's limitedLinearV and projected incident flux.
- **Face part**: G_f = snGrad(pHydI_rgh + pDyn) + ghf snGrad(rhoI), split three ways and each reconstructed:
  - pHyd: snGrad(pHydI);
  - gravity: the rest of the hydrostatic part;
  - pDyn: snGrad(pDyn).
- **visc**: laplacian(mu, UI).

## 1. Terms of rI in band cells (N/m³)

| Term | ZD max | ZD rms | B05 max | B05 rms |
|---|--:|--:|--:|--:|
| ddt(rhoI uI) | 1558 | 726 | 1556 | 875 |
| div(mI uI) | 1272 | 543 | 989 | 477 |
| pHyd (face) | 7245 | 5348 | 7631 | 5201 |
| gravity (face) | 7237 | 5340 | 7667 | 5152 |
| pHyd + gravity | 794 | 292 | 1179 | 576 |
| pDyn (face) | 408 | 145 | 291 | 171 |
| viscous | 0.1 | 0.0 | 0.0 | 0.0 |
| ddt + conv | 709 | 339 | 1005 | 578 |
| ddt + conv + pHyd + gravity | 661 | 290 | 850 | 409 |
| **total rI** | 476 | **231** | 888 | **397** |
| total, x component | 452 | 154 | 887 | 379 |
| total, z component | 276 | 173 | 290 | 119 |

## 2. Column level error etaHOS − etaI (etaI = bottom + Σ alphaI dz on the moved mesh)

| Case | Columns | max \|.\| (mm) | rms (mm) | mean (mm) |
|---|--:|--:|--:|--:|
| ZD | 120 | 7.8 | 4.4 | −1.4 |
| B05 | 39 | 44.6 | 19.0 | −15.7 |

In B05, alphaI is 20-35 mm above the HOS surface in the columns next to the box (19 s); away from the box it is within 2-8 mm.

## 3. Band G_f against rho_w g (etaHOS − etaI)/dz, column by column

The band G_f per column is the largest-|.| vertical face value in the band, signed upward.

| Case | Slope | R² | n |
|---|--:|--:|--:|
| ZD | −0.05 | 0.02 | 120 |
| B05 | 0.06 | 0.23 | 39 |

## Verdict

- **pDyn is not the cause.** It is the smallest pressure term (rms 145-171), the slope is about 0.05 rather than 1, and R² is 0.02-0.23.
- **Dominant term: unbalanced incident inertia.** ddt + conv (rms 339 / 578) is balanced only partly by the lateral hydrostatic pressure of the discrete alphaI surface (pHyd + gravity, 292 / 576), and pDyn takes off little more.
- In B05 the residual is mainly horizontal (x 379 against z 119) and sits where alphaI is 16 mm (mean) off the HOS surface. The analytic acceleration follows the HOS surface, while the hydrostatic slope follows the advected alphaI.

# alphaI drift at the box wall (no solver run), 23 Nov 2026

## 1. Code

- **alphaI on the box is zeroGradient** (`zonefix/B05/0.orig/alphaI`, patch box), so its face value is the wall cell's own alphaI.
- **The incident flux through the box wall is not zero.**
  - The raw flux is `phiRaw = UfInc & Sf` ([swenseIncidentStep.H:19](../src/swenseWaveTank/solver/swenseIncidentStep.H#L19)). UfInc on boundary faces is the blended incident velocity at the face centre ([swenseIncident.C:935-943](../src/swenseWaveTank/swenseIncident.C#L935-L943)).
  - psiProj is zeroGradient on every patch except the atmosphere ([swenseSetup.H:394-405](../src/swenseWaveTank/solver/swenseSetup.H#L394-L405)), so the projection flux on the box is 0 and `phiInc_box = phiRaw_box`.
  - `phiIRel = phiInc - meshPhi` ([swenseIncidentStep.H:57](../src/swenseWaveTank/solver/swenseIncidentStep.H#L57)) therefore carries the incident wave through the wall, as if there were no body.
- **exactWallFlux does not change this.** It only sets Uc on the box so that phiComp + phiInc = meshPhi ([swenseBodyVelocityFvPatchVectorField.C:121-137](../src/swenseWaveTank/fvPatchFields/swenseBodyVelocity/swenseBodyVelocityFvPatchVectorField.C#L121-L137)). phiInc and phiIRel are unchanged.
- **Where phiIRel flows from the wall into the fluid**, alphaI is advected by `advectorI(alphaI, phiIRel, UInc)` (createSolverFields.H:7; alphaStep.H:14-19).
  - isoAdvector starts every face from `upwind(phi).flux(alpha)*dt` (OpenFOAM isoAdvectionTemplates.C:456). On a boundary face that uses the patch value, which is the cell's own alphaI because the patch is zeroGradient.
  - The isoface correction is applied only to outflow boundary faces (`phiP >= 0`, isoAdvection.C:383).
  - So the wall imports alphaI equal to the cell's own value on inflow and exports through the isoface on outflow: **the box wall acts as a source and sink of alphaI.**

## 2. Fixed box against moving box: etaHOS − etaI next to the box

- **Fixed:** S6a_rev1 (Rev1, box locked by point and orientation constraints, meshPhi uniform 0, HOS wave). Box |x| < 0.234; samples 0.02-0.30 m from each box face. Script diag/eta_wall_x.py.
- **Moving:** B05 (columns within 0.3 m of the box). Script diag/eta_wall.py.

| Run | t (s) | Next to box, rms (mm) | Mean (mm) | Mean, x < box / x > box (mm) | Control x = ±1 (mm) |
|---|--:|--:|--:|--:|--:|
| S6a, fixed | 8.0 | 0.7 | −0.6 | −0.6 / −0.5 | +0.7 / +0.9 |
| S6a, fixed | 15.0 | 2.1 | −0.2 | −0.7 / +0.2 | −1.2 / +2.6 |
| S6a, fixed | 19.0 | 8.5 | +7.3 | +4.6 / +10.0 | −0.8 / +6.3 |
| S6a, fixed | 23.0 | 13.1 | +11.6 | +6.0 / +17.2 | −2.5 / −0.2 |
| B05, moving | 15 | 12.3 | −11.1 | −8 / −16 | rms 4.8 |
| B05, moving | 19 | 18.5 | −16.2 | −10 / −27 | rms 4.4 |
| B05, moving | 22.9 | 24.3 | −19.8 | −9 / −37 | rms 5.1 |

**Verdict: (a) the wall boundary.**
- The drift also grows with the box fixed: S6a rms 0.7 → 13 mm next to the box, while the controls stay at or below 6 mm. It is largest on the side x > box in both runs.
- The moving mesh roughly doubles it (B05 24 mm against S6a 13 mm at 23 s), but it is not the origin.
- The sign differs between the runs (S6a alphaI below HOS, B05 above). S6a is a different wave (HOS) and a different box.

## 3. B05 per-column alphaI budget over one period (diag/wall_budget.py)

- **Wall flux** = Σ over box faces of −alphaI_cell·(UInc_b·Sf − meshPhi_b), trapezoid over the 0.1 s writes.
- **Total** = Δ Σ alphaI V. **Internal** = total − wall.
- Values are given as level over one column (0.05 × 0.2 m), in mm per period.

| Column (set-up x) | 15-17 s: total / wall / internal | 19-21 s: total / wall / internal |
|---|---|---|
| −0.125, side, x < box | +2.2 / −120.5 / +122.7 | +1.1 / −110.2 / +111.2 |
| −0.075, over/under box | +4.5 / +8.4 / −3.9 | +3.3 / +10.7 / −7.4 |
| −0.025, over/under box | +4.1 / +3.2 / +0.9 | +2.5 / +5.1 / −2.6 |
| +0.025, over/under box | +3.7 / −2.0 / +5.7 | +1.8 / −0.4 / +2.2 |
| +0.075, over/under box | +3.3 / −7.1 / +10.4 | +1.0 / −6.0 / +6.9 |
| +0.125, side, x > box | +5.8 / +164.2 / −158.4 | +4.6 / +153.5 / −148.9 |

- The wall faces carry the incident wave's net mass transport through the box: 110-160 mm of level per period, out on the side x < box and in on the side x > box. The internal faces pass almost all of it on.
- What remains is a gain of 1-6 mm per period in every column next to the box. That is the size of the drift growth (B05 side x > box: −16 → −37 mm in 4 periods).
- The wall flux uses the cell alphaI, not isoAdvector's isoface value on outflow faces, so the split is an estimate. The totals are exact from the fields.
