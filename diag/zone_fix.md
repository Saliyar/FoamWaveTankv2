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

# alphaI boundary condition swenseIncidentAlpha on the box (23 Nov 2026)

## 1. Code

- **New BC** `swenseIncidentAlpha` ([swenseIncidentAlphaFvPatchScalarField.C](../src/swenseWaveTank/fvPatchFields/swenseIncidentAlpha/swenseIncidentAlphaFvPatchScalarField.C)), a mixed type.
  - Faces with `phiIRel < 0` (from the wall into the fluid) get valueFraction 1 and refValue = wetted area / face area below the ramped incident eta ([C:172-176](../src/swenseWaveTank/fvPatchFields/swenseIncidentAlpha/swenseIncidentAlphaFvPatchScalarField.C#L172-L176)).
  - The cut is `swenseSurface::wetArea`, the exact edge-cut of the face on the current `mesh.points()`, so it follows the moving mesh. It is evaluated at t^{n+1/2} = t − dt/2 ([C:146](../src/swenseWaveTank/fvPatchFields/swenseIncidentAlpha/swenseIncidentAlphaFvPatchScalarField.C#L146)).
  - Outflow faces get valueFraction 0 and refGrad 0, which is zeroGradient ([C:187](../src/swenseWaveTank/fvPatchFields/swenseIncidentAlpha/swenseIncidentAlphaFvPatchScalarField.C#L187)), so isoAdvector's outflow treatment is unchanged.
- **Switch** `incidentAlphaWall` in waveProperties, default off ([swenseSetup.H:805](../src/swenseWaveTank/solver/swenseSetup.H#L805)).
  - On: every body patch of alphaI must be swenseIncidentAlpha, otherwise FatalError.
  - Off: the BC is exactly zeroGradient.
- **Evaluation.** alphaStep.H evaluates the body patches after `alphaI = alphaI.oldTime()` and before `advectorI.advect` ([alphaStep.H:22-58](../applications/solvers/waveIsoFoamSWENSE/alphaStep.H#L22-L58)). It logs the inflow faces and the alphaI flux in and out.
- **Where isoAdvector takes the inflow value** (OpenFOAM v2506):
  - `dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_)*dt` (isoAdvectionTemplates.C:456). On a non-coupled boundary face the interpolate takes the patch value, `sfbf[pi] = vf.boundaryField()[pi]` (surfaceInterpolationScheme.C:194).
  - `timeIntegratedFlux` overwrites only boundary faces with `phiP >= 0` (isoAdvection.C:383).
  - `limitFluxes`/`boundFlux` change only a cell's downwind faces (isoAdvectionTemplates.C:329).
  - So on an inflow face dVf = phiIRel·alphaI_b·dt, with alphaI_b the BC value.
- **Bit check.** With the switch off, B05 rerun to 0.3 s with the new build is identical to B05 (alphaI, alpha.water, Uc, pc and points cmp-equal).
- **Build** `private_ia`, the frozen copy used by L0, L1 and M1.

## 2. Gate L: B05 with the box locked (staticFvMesh), 23 s

Runs `zonefix/L0` (switch off) and `zonefix/L1` (switch on). Values are etaHOS − etaI in columns within 0.3 m of the box (diag/eta_wall.py).

| t (s) | L0, off: rms / mean (mm) | L1, on: rms / mean (mm) | L1 controls x = ±1: rms (mm) |
|--:|--:|--:|--:|
| 8 | 3.7 / −0.0 | **2.1** / −1.2 | 1.9 |
| 15 | 7.2 / +2.9 | **2.2** / −1.7 | 2.6 |
| 19 | 12.3 / +7.0 | **2.6** / −2.0 | 1.2 |
| 23 | 17.8 / +11.5 | **2.3** / −1.7 | 1.1 |

- **L passes.**
  - L1 stays at 2-3 mm with no growth; the largest value at 23 s is 4 mm. That is below the 6 mm control level.
  - It is still above L1's own controls (1-2.6 mm).
- **L0 confirms (a) on the B05 geometry.** The drift grows to 17.8 mm with the box fixed, and the side x > box reaches +15 to +34 mm at 23 s.

## 3. Gate M: B05 (moving, zone 0.5 m, discrete hydrostatics) with the BC, 23 s

Run `zonefix/M1`.

| t (s) | B05: rms / mean (mm) | M1: rms / mean (mm) | M1 controls: rms (mm) |
|--:|--:|--:|--:|
| 8 | — | 2.8 / −2.5 | 2.9 |
| 15 | 12.3 / −11.1 | 6.0 / −5.0 | 4.4 |
| 19 | 18.5 / −16.2 | **7.1** / −5.9 | 4.2 |
| 23 | 24.5 / −20.5 | **8.0** / −6.4 | 4.6 |

Per column at 23 s:

| Side | B05 (mm) | M1 (mm) |
|---|--:|--:|
| x > box | −27 to −45 | −2 to +1 |
| x < box | −9 to −11 | −8 to −12 |

**M fails the drift gate** (8.0 mm against 6 mm at 23 s).
- The BC removes the error on the side x > box.
- The side x < box keeps the B05 error. With the box locked that side is about −1 mm (L1), so this residual comes with the mesh motion.
- On that side the wall faces are mostly outflow (phiIRel into the wall), which the BC leaves to isoAdvector.

**Remaining wall budget**, M1 per column, level in mm per period. The wall part uses the written box patch value.

| Column | 15-17 s: total / wall / internal | 19-21 s: total / wall / internal |
|---|---|---|
| −0.125, side, x < box | +2.3 / −157.1 / +159.4 | +2.1 / −157.0 / +159.1 |
| −0.075, over/under | +1.3 / +10.2 / −8.8 | +0.5 / +10.2 / −9.6 |
| −0.025, over/under | +1.1 / +4.6 / −3.4 | +0.6 / +4.9 / −4.3 |
| +0.025, over/under | +0.9 / −0.9 / +1.8 | +0.7 / −0.1 / +0.8 |
| +0.075, over/under | +0.7 / −6.0 / +6.7 | +0.7 / −4.8 / +5.5 |
| +0.125, side, x > box | +1.5 / +137.1 / −135.6 | +1.2 / +129.6 / −128.4 |

- The net gain falls from +1 to +6 mm per period (B05) to +0.5 to +2.3.
- The largest remaining gain is the column on the side x < box.

**Motion amplitude, last 5 periods (13.04-23 s), error against plain fine | against experiment**

| Run | Surge | Heave | Pitch | Tension |
|---|--:|--:|--:|--:|
| SWENSE test 2 (zone 1.0 m) | +1.6 \| +12.8 % | +4.8 \| +11.8 % | −19.0 \| −12.4 % | +14.0 \| +1.7 % |
| SWENSE B05 | +7.2 \| +19.0 % | +4.2 \| +11.2 % | −62.3 \| −59.3 % | +32.4 \| +18.1 % |
| **SWENSE M1 (B05 + BC)** | +12.0 \| +24.3 % | +7.2 \| +14.4 % | **−39.5 \| −34.6 %** | +26.0 \| +12.4 % |

**Wave against plain fine, rms (mm)**

| Run | x = −1 | x = −0.25 | x = +0.25 | x = +1 | Mean at ±1 | Mean at the box |
|---|--:|--:|--:|--:|--:|--:|
| SWENSE test 2 | 8.8 | 12.7 | 10.8 | 10.0 | 9.4 | 11.7 |
| SWENSE B05 | 7.6 | 23.6 | 19.6 | 8.7 | 8.2 | 21.6 |
| **SWENSE M1** | 7.9 | 11.4 | 10.3 | 6.4 | 7.1 | **10.9** |

- The wave at the box halves (21.6 to 10.9 mm), now just below test 2.
- Pitch improves from −62 % to −40 % but is still worse than test 2 (−19 %).
- Surge gets worse (+12 % against +7 %).

## 4. Zone 1.0 m

Not run, because M fails the drift gate. The case `zonefix/M10` (bodyZone |x| < 0.5, BC on) is prepared but not submitted.

Scripts: diag/eta_wall.py (times as an argument), diag/wall_budget.py (case as an argument, uses the patch value), diag/s4.py (motions and waves; extra runs as arguments).

# Isoface velocity, incidentResetInZone, M10 (23 Nov 2026)

## 1. Isoface velocity on a moving mesh (code, no run)

isoAdvector moves the isoface with the **absolute** velocity it is given; the mesh velocity is never subtracted.
- `interpolationCellPoint<vector> UInterp(U_)` (OpenFOAM v2506, isoAdvection.C:248) and `Un0[i] = UInterp.interpolate(x0, celli) & n0` (isoAdvection.C:278).
- The face-crossing times use the current points, `((mesh_.points()[fi] - x0) & n0)/Un0` (cutFaceAdvect.C:237). The volume is scaled by the flux passed in, `dVf = phi/magSf*timeIntegratedArea(..., Un0)` (cutFaceAdvect.C:264).
- So the flux is relative to the mesh while the isoface velocity is absolute, in all three cases:

| Solver | Field | Flux | Velocity |
|---|---|---|---|
| waveIsoFoamSWENSE | alphaI | `phiIRel` = phiInc − meshPhi | `incident.UInc()`, absolute (createSolverFields.H:7) |
| waveIsoFoamSWENSE | alpha1 | `phi` = phiComp + phiInc − meshPhi | `U == UInc + Uc`, absolute (createSolverFields.H:3, alphaStep.H:11) |
| plain waveIsoFoam | alpha1 | `phi` after `fvc::makeRelative(phi, U)` (waveIsoFoam.C:161) | `U`, absolute (createFields.H:136) |

## 2. Switch incidentResetInZone (waveProperties, default off)

[swenseResetInZone.H](../src/swenseWaveTank/solver/swenseResetInZone.H) is included after the alphaI step, the relaxation and the incident step, in both the static and the moving path ([swenseTimeLoop.H:48](../src/swenseWaveTank/solver/swenseTimeLoop.H#L48), [:103](../src/swenseWaveTank/solver/swenseTimeLoop.H#L103)). It needs plainBodyZone ([swenseSetup.H:833](../src/swenseWaveTank/solver/swenseSetup.H#L833)).

In the body-zone cells it does the following:
- **alphaI** = V_wet/V below the ramped incident eta at t^{n+1}, on the current cell points ([:40](../src/swenseWaveTank/solver/swenseResetInZone.H#L40)). rhoI and Drho are rebuilt.
- **uI = mI/rhoI and pInc** come from the cell blend with the new fraction (`swenseIncident::blendCells`, [:51](../src/swenseWaveTank/solver/swenseResetInZone.H#L51)). The projection's cell correction of the unreset state is kept: UInc_B = uI_B − R/rhoI_B, with R = rhoI_A (uI_A − UInc_A) ([:91](../src/swenseWaveTank/solver/swenseResetInZone.H#L91)).
- **Uc** += UInc_A − UInc_B ([:93](../src/swenseWaveTank/solver/swenseResetInZone.H#L93)).
- **pc** += (pIncUse + Drho gh)_A − (pIncUse + Drho gh)_B in every cell, because p = pIncUse + pc + Drho gh (swensePcEqn.H:111). With the discrete hydrostatics, a change of rhoI also changes pIncD below it in the column. pc boundary values move with their face cells on patches that do not fix the value; fixedFluxPressure cannot be evaluated outside the pc equation.
- **Not touched:** phiInc, UfInc, phiIRel, phi and rho (from alpha1). The face incident stays that of the unreset step, and the projection is not re-run, so UInc outside the zone is unchanged.
- In the zone (lam = 0) the face force −ghf snGrad(Drho) − snGrad(pc) − G_f adds up to −ghf snGrad(rho) − snGrad(p_rgh). It is therefore invariant under the shift.

## 3. Check

**Reset instant.** The solver logs the change of the totals at every reset. Over all 5671 resets of M10:
- max |ΔU_tot| = 1.8e-15 m/s;
- max |Δp_tot| = 4.6e-13 Pa (p_tot up to 3.4e3 Pa).

This is exact to round-off.

**After the first step** (dt 0.00096 s; M10 set-up, reset on against off, cases C1on/C1off, diag/cmp_fields.py):

| Variant | U, max \|on − off\| (m/s) | p, max \|on − off\| (Pa) | alpha1 |
|---|--:|--:|--:|
| as run (pc tol 1e-7, limitedLinearV) | 4.6e-10 | 1.1e-5 | 0 |
| pc, Uc tolerances 1e-14 (converged) | 2.2e-10 | 1.2e-5 | 0 |
| tolerances 1e-14 + `Gauss linear` for div(rhoPhi,U/Uc) | 1.2e-11 | 3.8e-7 | 0 |

- So 1e-12 holds at the reset, not after the first solve.
- The gap is not the linear solvers. Most of it is the `limitedLinearV` limiter, which is computed separately from Uc (implicit) and UInc (explicit), so their sum is not the scheme applied to U.
- What remains with linear convection (1e-10 relative) is a smaller operator that is not linear in the split; not identified.
- At 0.3 s (C3on/C3off) the totals differ by U 2.4e-4 m/s, p 0.19 Pa, alpha1 2.4e-4.

**Bit check, switches off.** B05 to 0.3 s with this build (`private_rz`) is identical to B05 (alphaI, alpha.water, U, Uc, p, pc cmp-equal).

## 4. M10: zone 1.0 m, discrete hydrostatics + incidentAlphaWall + reset, 23 s

Run `zonefix/M10`, coarse, test 2 settings, wall time 326 s.
- Reset size: |ΔalphaI| up to 0.011 (median 0.002 per step).
- |ΔUInc| up to 0.78 m/s (median 0.009, 99th percentile 0.27): one interface cell whose blend moves between the water and the Air2Air velocity. Uc takes the opposite shift.

**Drift**, etaHOS − etaI in columns within 0.3 m of the box, rms / mean (mm):

| Run | t (s) | All | x < box | x > box | Controls x = ±1, rms |
|---|--:|--:|--:|--:|--:|
| test 2 | 8 | 3.0 / −2.7 | 3.5 / −3.4 | 2.6 / −2.3 | 2.8 |
| test 2 | 15 | 9.9 / −8.7 | 5.8 / −5.8 | 14.2 / −13.5 | 3.3 |
| test 2 | 19 | 14.7 / −12.3 | 6.5 / −6.4 | 22.3 / −21.7 | 3.4 |
| test 2 | 23 | 18.3 / −14.9 | 6.8 / −6.8 | 28.2 / −27.9 | 4.1 |
| M1 | 8 | 2.8 / −2.5 | 4.1 / −3.8 | 1.7 / −1.6 | 2.9 |
| M1 | 15 | 6.0 / −5.0 | 7.4 / −7.3 | 1.0 / −0.9 | 4.4 |
| M1 | 19 | 7.1 / −5.9 | 8.9 / −8.8 | 1.0 / −0.6 | 4.2 |
| M1 | 23 | 8.0 / −6.4 | 9.9 / −9.8 | 0.9 / −0.4 | 4.6 |
| **M10** | 8 | 1.0 / −0.9 | 0.8 / −0.8 | 1.1 / −1.0 | 2.8 |
| **M10** | 15 | 1.3 / −1.2 | 1.2 / −1.1 | 1.3 / −1.3 | 3.7 |
| **M10** | 19 | 1.3 / −1.2 | 1.2 / −1.1 | 1.3 / −1.3 | 4.2 |
| **M10** | 23 | 1.2 / −1.2 | 1.2 / −1.1 | 1.4 / −1.3 | 5.1 |

- **The M10 drift is forced, not a test.** These columns lie inside the 1.0 m zone, where the reset sets alphaI to the incident surface every step. What it shows is that the reset is applied.
- The independent measures are the controls (outside the zone) and the motions and waves below.

**Motion amplitude, last 5 periods (13.04-23 s), error against plain fine | against experiment**

| Run | Surge | Heave | Pitch | Tension |
|---|--:|--:|--:|--:|
| Plain coarse | −3.0 \| +7.7 % | −7.2 \| −1.0 % | −39.7 \| −34.8 % | +5.0 \| −6.4 % |
| SWENSE test 2 (zone 1.0 m) | +1.6 \| +12.8 % | +4.8 \| +11.8 % | −19.0 \| −12.4 % | +14.0 \| +1.7 % |
| SWENSE M1 (zone 0.5 m, BC) | +12.0 \| +24.3 % | +7.2 \| +14.4 % | −39.5 \| −34.6 % | +26.0 \| +12.4 % |
| **SWENSE M10 (zone 1.0 m, BC + reset)** | +6.8 \| +18.6 % | **+2.9 \| +9.8 %** | **−17.4 \| −10.7 %** | +20.7 \| +7.6 % |

**Wave against plain fine, rms (mm)**

| Run | x = −1 | x = −0.25 | x = +0.25 | x = +1 | Mean at ±1 | Mean at the box |
|---|--:|--:|--:|--:|--:|--:|
| Plain coarse | 9.4 | 8.0 | 11.5 | 9.8 | 9.6 | 9.7 |
| SWENSE test 2 | 8.8 | 12.7 | 10.8 | 10.0 | 9.4 | 11.7 |
| SWENSE M1 | 7.9 | 11.4 | 10.3 | 6.4 | 7.1 | 10.9 |
| **SWENSE M10** | 6.6 | 10.0 | 8.4 | 4.6 | **5.6** | **9.2** |

**Result.**
- M10 has the smallest wave error of all coarse runs, SWENSE or plain: 5.6 mm at ±1 and 9.2 mm at the box.
- Heave and pitch are closest to plain fine of the coarse SWENSE runs (pitch −17 % against −19 % for test 2 and −40 % for plain coarse).
- Surge (+6.8 %) and tension (+20.7 %) remain worse than test 2 (+1.6 %, +14 %).
- M10 differs from test 2 in three switches at once (discrete hydrostatics, incidentAlphaWall, reset), so the gain cannot be assigned to the reset alone.

Scripts: diag/cmp_fields.py (on/off field difference), diag/eta_wall.py (side split, case path), diag/s4.py. Build `private_rz`.

## alphaHExact and the hybrid mesh (H_P, H_S)

### 1 Switch alphaHExact (waveProperties, default off)
- With the switch on, alphaH is the exact wetted fraction below the ramped etaHOS: the current cell polygon is cut by `swenseSurface::wetVolume`, the same routine incidentResetInZone uses. Code: `swenseIncident::exactAlpha` / `targetAlpha` / `alphaHExact` (src/swenseWaveTank/swenseIncident.C).
- The same target is used in all three places alphaH is targeted:
  - the relaxation target (swenseTimeLoop.H, both paths, plus the initial alphaI in swenseSetup.H);
  - the wall BC swenseWaveAlpha (exact wetted face fraction, `wetArea`);
  - the zone reset, which already used the cut.
- Log line: `SWENSE alphaHExact: max|exact - blend|` (solver/swenseAlphaHLog.H).
- With the switch off, OFF3 is bit-identical to B05 at 0.3 s (alphaI, alpha.water, U, Uc, p, pc).

| Check | max \|exact − blend\| |
|---|--:|
| Still water, coarse (SX0), initial | 2.2e-10 (blendThickness 0.0438666667 is dz = 0.1316/3 to 10 digits only) |
| Still water, hybrid (blendThickness 0.01316 exact), initial | 9.2e-15 |
| Waves, coarse M10X, 23 s | max 0.056, median per step 0.048 |
| Waves, hybrid H_S, 23 s | max 0.067, median 0.052, p95 0.063 |

- On coarse, the error jumps from 0.009 to 0.051 at t ≈ 6 s. That is when the ramped crest enters the 0.049 m rows above z = 0.053, and the blend assumes dz = delta.
- M10X (M10 plus alphaHExact) matches M10 to within 0.4 % on every motion and 0.1 mm on the waves. On its own the switch changes little.

### 2 Hybrid mesh (diag/mkHybrid.py, conformal blockMesh, no hanging nodes)
- **x:** the fine dx 0.0125 for |x| < 0.6 m (16 cells across the box), then 14 cells graded at q = 1.097 out to |x| = 0.975, then coarse dx 0.049 out to ±2.5. That gives 186 columns; the largest ratio is 1.097.
- **z:** the fine rows (5 × 0.0121 below the box, 10 × 0.0132 box rows, 5 × 0.0147 above) run from z = −0.079 (below the box bottom) to 0.127 (above the max crest).
  - Below them, 14 graded rows run down to the floor; above them, 15 graded rows reach z = 0.455; then 3 top rows of 0.015.
  - The same rows run across the tank, 52 in total. The largest ratio is 1.117, and only at the fine rows' own step from 0.0132 to 0.0147.
- **checkMesh:** OK, non-orthogonality 0, max aspect ratio 4.05.
- **Mesh motion:** moving points lie at |x| ≤ 0.49 (inside the fine part). Of the 1496 points in the top 3 rows (z ≥ 0.455), none moves.
- **Cells:** hybrid 9512; coarse 1088 (8.7×); fine 15040 (0.63×).

### 3 Gauges
- x = ±0.25 (box): linear between the two columns that straddle x at t = 0, moving with the mesh.
- x = ±1 (graded band) and **x = ±1.25 (coarse part, |x| ≥ 0.975, outside the SWENSE relaxation zones |x| > 1.5)**: linear in x between static columns.
- The box sampling changed from s4.py, where nearest-column tie-breaking could double the value (plain coarse at −0.25: 8.0 against 16.0 mm). The box and ±1 values below therefore differ from the M10 table above. diag/s5.py recomputes all rows.

### 4 Runs (23 s, serial)
- **Ramp:** 12 s in both runs, but the shape differs. H_P uses the waveModels inlet ramp (linear t/12, inlet only, no relaxation zones). H_S uses rampTsoft 12 (sin(π/2·t/12), everywhere).
- **H_S:** SWENSE with discrete hydrostatics, incidentAlphaWall, incidentResetInZone and alphaHExact on; zone 1.0 m; blendThickness 0.01316; dt 2e-4; build `private_hx`.
- **H_P:** plain waveIsoFoam (StokesII inlet), still-water start from setAlphaField (plane z = 0, normal (0 0 −1)).
- **Initial eta error** (diag/eta0.py, column integral, plus a sidedness check):
  - H_P: 7.8e-16 m.
  - H_S after one step (t = 2.4e-4 s, includes the ramped wave): 1.1e-6 m.
- **First H_P attempt failed.** With normal (0 0 1) the water sat above z = 0. The column integral cannot see this (−0.5 + 0.5 = 0 either way): the box fell and MoorDyn hit SIGFPE at 0.24 s. eta0.py now checks the sidedness.

### 5 Result against plain fine

| Run | Cells | Surge | Heave | Pitch | Tension | Box ±0.25 (mm) | ±1 (mm) | Coarse ±1.25 (mm) | Core-h |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Plain fine | 15040 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5.25 |
| Plain coarse | 1088 | −3.0 % | −7.2 % | −39.7 % | +5.0 % | 13.1 | 9.4 | 9.6 | 0.07 |
| **H_P** (plain, hybrid) | 9512 | **−0.6 %** | **−1.5 %** | **+0.7 %** | **−0.9 %** | **4.6** | **1.3** | **1.5** | 3.01 |
| H_S (SWENSE, hybrid) | 9512 | +6.0 % | +4.2 % | +0.6 % | +19.6 % | 11.9 | 5.3 | 6.9 | 1.14 |
| (SW M10, coarse) | 1088 | +6.8 % | +2.9 % | −17.4 % | +20.7 % | 11.2 | 5.3 | 6.4 | 0.09 |

Against experiment: plain fine +11.0 / +6.7 / +8.1 / −10.8 %; H_P +10.3 / +5.0 / +8.9 / −11.6 %; H_S +17.7 / +11.2 / +8.8 / +6.7 %.

**Result.**
- **Plain waveIsoFoam on the hybrid mesh reproduces plain fine:** every motion within 1.5 %, waves within 1.3 mm at ±1. It costs 0.57× the core-hours (3.0 against 5.25 h).
- **SWENSE on the same mesh gets the pitch right (+0.6 %), but surge and tension carry the same bias as on coarse** (+6 %, +20 %). Its waves are no better than M10 on the coarse mesh: 5.3 mm at ±1, 6.9 at ±1.25, 11.9 at the box.
- The SWENSE error is therefore not set by the mesh resolution at the box. A refined mesh fixes plain but not SWENSE.
- H_S is 2.6× cheaper than H_P (dt 2e-4 fixed against H_P's adaptive dt from 1.25e-3).

Scripts: diag/mkHybrid.py (blockMeshDict generator), diag/eta0.py (still-water level and sidedness), diag/s5.py (table; extra runs as label=path). Build `private_hx`.

## Harmonics, mean loads and RAOs from existing runs (no solver run; `diag/s6.py`)

Runs: plain fine, H_P and H_S. The fit window is 13–23 s, with per-period windows from 8 s. Each signal is fit by least squares to m + A1 at ω + A2 at 2ω (ω = 2π/1.993). The target is Stokes II, A1 = 59.75 mm and A2 = 12.94 mm, with phase kx + 3π/2 at absolute x. ESI StokesII and StokesIIIncident use the same origin and phase, so both solvers share this target.

**The gauge method changed.** eta is now −h plus ∑ α × (exact chord length of the vertical line x = x_g through each cell quad on the current mesh). Every column checks to 1.000000 m at 8 s and at 23 s. The old s5 method used zmax − zmin of each cell. That overstates the height of sheared cells in the moving zone and gave the ±0.25 gauges a +35 mm mean bias. So the s5 box rms numbers carry a pitch-dependent bias; the static gauges (±1, ±1.25) are unaffected.

**Waves, 13–23 s.** Columns: A1 in mm (ratio to target), phase relative to target, A2 in mm, mean in mm. "box" is the average of the ±0.25 gauges.

| Run | x = −1.25 | −1 | box | +1 | +1.25 |
|---|---|---|---|---|---|
| plain fine | 59.6 (0.997), −1°, 12.1, +1.1 | 61.2 (1.023), −2°, 12.1, +1.3 | 53.1 (0.961), −8°, 10.1, +0.4 | 60.6 (1.015), −6°, 15.8, +1.5 | 61.8 (1.034), −9°, 15.4, +1.5 |
| H_P | 59.8 (1.002), −2°, 11.6, +0.9 | 60.9 (1.020), −4°, 11.9, +1.1 | 52.3 (0.945), −8°, 10.1, +0.3 | 60.3 (1.008), −8°, 15.7, +1.3 | 60.8 (1.017), −11°, 15.0, +1.3 |
| H_S | 65.6 (1.097), +1°, 10.9, −2.9 | 64.5 (1.080), −3°, 13.9, −3.3 | 54.2 (0.981), −3°, 14.2, −2.8 | 63.4 (1.062), −5°, 15.9, +1.9 | 63.4 (1.061), −6°, 15.8, +2.3 |

At the box the A2 target is 9.2 mm.

**Per-period |A1|/target** (windows start at 8.00, 9.99, 11.99, 13.98, 15.97, 17.96 and 19.96 s):

| Run | Gauge | Ratio |
|---|---|---|
| plain fine | −1.25 | 0.996 0.998 0.999 0.997 0.997 0.996 0.999 |
| plain fine | box | 0.951 0.943 0.946 0.947 0.948 0.948 0.947 |
| plain fine | +1.25 | 1.027 1.028 1.030 1.030 1.031 1.034 1.034 |
| H_P | −1.25 | 0.695 0.859 0.987 0.998 1.000 1.002 1.002 |
| H_P | box | 0.617 0.777 0.916 0.945 0.935 0.930 0.931 |
| H_P | +1.25 | 0.582 0.754 0.919 1.014 1.020 1.021 1.020 |
| H_S | −1.25 | 0.999 1.089 1.095 1.096 1.098 1.090 1.101 |
| H_S | box | 0.890 0.974 0.984 0.986 0.985 0.913 0.982 |
| H_S | +1.25 | 0.908 1.018 1.049 1.051 1.050 1.062 1.072 |

H_P's ramp finishes at the 13.98 s window (its linear inlet ramp plus travel time). H_S has reached its plateau by the 9.99 s window, but that plateau sits 9–10 % above target upstream and 5–7 % above target downstream. The fit window (13 s onward) is therefore past both ramps.

**Means, 13–23 s.**

| Run | Surge (mm) | Tension − T0 (N) | Pitch (deg) | Heave (mm) |
|---|--:|--:|--:|--:|
| plain fine | +8.10 | +0.232 | +0.67 | +1.2 |
| H_P | +8.08 | +0.230 | +0.64 | +0.8 |
| H_S | +9.62 | +0.276 | +1.13 | −0.4 |

Gap from H_S to plain fine: +1.52 mm surge and +0.043 N tension. Their ratio, ΔT/Δx = 28.6 N/m, equals the mean stiffness of plain fine (0.232/8.10 mm = 28.7 N/m), so the gap is a steady extra force of about 0.043 N.

**Stokes drift with no return current.** In the closed tank the return current is U = E/(ρch) = g H²/(8ch) = 17.3 mm/s, and the SWENSE incident lacks it. The drag of that current on the box is ½ρ C_d (draft × 0.2 m) U² = 0.002–0.005 N for C_d = 1–2. The Doppler change in wave drift is about 4U/c × 0.23 N = 0.008 N. Both are 5–20 times smaller than 0.043 N. The tension ratio, 0.276/0.232 = 1.19, is close to the ratio of squared amplitudes upstream, (65.6/59.6)² = 1.21. The mean-load gap therefore follows H_S's 10 % larger upstream wave, not a missing return current.

**RAOs.** Each motion's A1 is divided by the box-gauge A1; phase is the motion minus the box gauge. Values in brackets are relative to plain fine.

| Run | Surge (m/m) | Heave (m/m) | Pitch (deg/m) | Tension (N/m) |
|---|---|---|---|---|
| plain fine | 2.271, −88° | 1.120, −2° | 45.7, +19° | 11.13, −54° |
| H_P | 2.295 (+1.1 %, +1.9°) | 1.121 (+0.2 %, +0.1°) | 46.6 (+1.8 %, +4.5°) | 11.30 (+1.5 %, +1.4°) |
| H_S | 2.383 (+4.9 %, +3.8°) | 1.134 (+1.3 %, −3.6°) | 52.5 (+14.9 %, −5.4°) | 12.96 (+16.4 %, +3.1°) |

Pitch is mostly 2nd harmonic. Plain fine has h1 2.43° and h2 3.21°; H_S has h1 2.85° (+17 %) and h2 2.64° (−18 %). Their peak-to-peak values agree (4.70° against 4.72°) only because the two errors offset. Tension in H_S has h1 0.703 N against 0.591 N, and h2 0.201 N against 0.152 N.

**Verdict.** Only heave meets the 3 % test. The mean-load gap is 5–20 times larger than any return-current estimate. So the difference is not explained by the incident definition alone. The RAOs that differ are pitch (+15 % at 1st harmonic, −5°, with 2nd harmonic −18 %) and tension (+16 %); surge is borderline at +4.9 %. H_P matches plain fine within 2 % on every RAO, so the mesh is not the cause.

## Incident amplitude: etaI against eta, the fluxConsistentU switch, H_S2

### A. Splitting H_S into incident and scattered parts (`diag/s7.py`, s6 line cut, 13–23 s)

First-harmonic amplitude relative to the Stokes II target, with phase relative to target in brackets:

| Run | Field | x = −1.25 | −1 | +1 | +1.25 |
|---|---|---|---|---|---|
| H_S | eta | 1.097 (+1°) | 1.080 (−3°) | 1.062 (−5°) | 1.061 (−6°) |
| H_S | etaI | 1.034 (+1°) | 1.033 (0°) | 1.037 (−1°) | 1.038 (0°) |
| H_S | etaS = eta − etaI | 0.063 | 0.069 | 0.087 | 0.108 |
| F3b (coarse, empty) | eta = etaI | 1.027 (+6°) | 1.026 (+5°) | 1.040 (+5°) | 1.039 (+6°) |
| H_P | \|eta − target\| | 0.042 | 0.067 | 0.136 | 0.191 |

etaI is 3–4 % high even in the empty tank, so neither verdict branch applies cleanly (≈1.0 or 1.06–1.10). The high incident is the part that can be fixed and gated, so the work went on to B. The excess is not built up by advection: in H_S, `diag/s7x.py` shows A1(etaI) jumping from 1.008 at x = −1.8 to 1.036 at x = −1.6 (the inner edge of the relaxation zone), then staying flat across the interior.

### B. Source of the excess

- **Amplitude input:** `src/waves2Foam/waveTheories/regular/StokesIIIncident/StokesIIIncident.C:78` reads H = waveHeight (0.1195) from the inlet dictionary. The log confirms H 0.1195, L 4.0396, K 1.5554 and depth 0.5.
- **eta:** the same file at :197–198 computes H/2 cos θ + kH²/4 (3 − σ²)/(4σ³) cos 2θ. The coefficients match ESI StokesIIWaveModel.C:71-73.
- **Velocity:** :250-256 (u, w) use H/2 ω cosh(kz)/sinh(kh) + (3/16) H² ω k cosh(2kz)/sinh⁴(kh), the same as ESI.
- **Ramp:** `swenseIncident.C:535` applies `swenseRamp::factor`, which is 1 after rampTsoft = 12 s.
- **alphaH:** the relaxation target, the wall-BC target and the zone-reset target all come from `targetAlpha` (the analytic eta, exact cut when alphaHExact is on; `solver/swenseTimeLoop.H:36, 98`).
- **alphaI:** in the interior (w = 0), alphaI is not set from the analytic eta. It is transported by the projected incident flux (`alphaStep.H`, then `solver/swenseAlphaRelax.H:9`).

**Cause.** Stokes II is consistent only to second order. Integrating its velocity up to the free surface gives a column flux whose first harmonic is too large: Q₁ = ∫₋ₕ^η u dz has Q₁/(c·H/2) = 1.0342 for this wave (analytic; 1.0 when integrated to z = 0). Transported alphaI then settles at η₁ = Q₁/c, which is 3.4 % above H/2, matching the measured 1.033–1.038.

**Fix (switch, default off).** `fluxConsistentU true;` in `waveTheoryForcingCoeffs` multiplies the first-order u and w terms by fU, chosen so that Q₁ = c·H/2. Both column integrals are analytic in z and averaged over 3600 phase points (`StokesIIIncident.C:135-164`). For BOX3, fU = 0.96639. eta and the pressure are not changed.

libwaves2Foam is now built privately from this source. Switch off, F3c0 against F3b: alpha and U at 15 s are bit-identical.

**Gate: coarse empty tank, 15 s, switch on (F3c), window 13–15 s.**

| Run | etaI = eta at x = −1 | x = +1 | Uc |
|---|--:|--:|--:|
| F3c (switch on) | 0.994 | 1.006 | exactly 0 |
| F3c0 (switch off) | 1.030 | 1.041 | exactly 0 |

**The gate passes.** The phase stays at +5° on the coarse mesh; it is not gated and is unchanged by the switch.

### C. H_S2 = H_S plus fluxConsistentU (hybrid mesh, 23 s requested)

**H_S2 crashed at 22.69 s.** Pitch departed at 22.1 s and reached −18.8° by 22.45 s, against a normal range of −3° to +7°. The time step collapsed to 1e-87, and MoorDyn raised a floating-point exception in Line::GetLineOutput. H_S, which has the same set-up with the switch off, ran to 23 s. Up to 22.1 s the H_S2 motions were regular and smaller than H_S. The run used 6319 s to 22.69 s, against 4124 s for H_S over 23 s.

All rows below use the window 13–20.97 s (4 whole periods).

**Waves:** first-harmonic amplitude relative to target, at x = −1.25 / −1 / box / +1 / +1.25:

| Run | −1.25 | −1 | box | +1 | +1.25 |
|---|--:|--:|--:|--:|--:|
| plain fine | 0.997 | 1.023 | 0.961 | 1.014 | 1.033 |
| H_P | 1.002 | 1.019 | 0.948 | 1.007 | 1.016 |
| H_S | 1.098 | 1.080 | 0.981 | 1.063 | 1.059 |
| H_S2 | 1.055 | 1.039 | 0.967 | 1.008 | 1.008 |

In H_S2, etaI is 0.998 / 0.998 / 1.003 / 1.002 at ±1.25 and ±1, and etaS is 0.06–0.08 (H_S had 0.07–0.11). H_S2 mean levels are −1.3 mm upstream and +0.6 to +0.9 mm downstream (H_S: −3.2 and +2.3 mm).

**Motions and means.** Amplitude columns are peak-to-peak errors against plain fine.

| Run | Surge | Heave | Pitch | Tension | Mean surge (mm) | Mean tension (N) | Mean pitch (deg) |
|---|--:|--:|--:|--:|--:|--:|--:|
| plain fine | 0 | 0 | 0 | 0 | 8.07 | 0.2325 | 0.68 |
| plain coarse | −3.0 % | −7.3 % | −40.3 % | +4.8 % | | | |
| H_P | −0.8 % | −1.5 % | 0.0 % | −1.1 % | 8.00 | 0.2288 | 0.64 |
| H_S | +6.5 % | +4.4 % | −0.7 % | +21.6 % | 10.37 | 0.2821 | 1.11 |
| H_S2 | −0.4 % | −0.5 % | −12.4 % | +6.0 % | 7.82 | 0.2308 | 0.85 |

**RAOs** (first harmonic over the box-gauge first harmonic), against plain fine:

| Run | Surge | Heave | Pitch | Tension |
|---|---|---|---|---|
| H_P | +0.6 %, +2° | 0.0 %, 0° | +0.7 %, +5° | +1.0 %, +1° |
| H_S | +5.3 %, +3° | +1.6 %, −4° | +12.3 %, −7° | +17.9 %, +2° |
| H_S2 | −0.5 %, +7° | −0.2 %, −1° | −8.0 %, +10° | +4.5 %, +6° |

**Conclusion.** The flux-consistent incident removes the mean-load gap (surge and tension within 3 % and 1 % of plain fine) and brings the surge and heave RAOs within 0.5 %. Two things remain. Pitch now lies on the other side of plain fine (RAO −8 %, +10°; peak-to-peak −12 %), and tension is +4.5 % in RAO and +6 % peak-to-peak. A scattered wave of 6–8 % also remains upstream. And H_S2 is unstable at 22.1 s.

## H_S2 blow-up at 22.1 s: what departs first (no solver run; `diag/s8.py`, `diag/s8f.py`)

**Method.** The data are the solver log at every step from 20.0 to 22.7 s (H_S: 2862 steps; H_S2: 14992 steps), the 4 MoorDyn line files, and the field writes every 0.1 s, with max-U cell locations taken from those writes. A quantity "departs" at the first time |H_S2 − H_S| exceeds 5× its rms over 20.0–21.5 s and stays above that for 0.02 s.

**Order of departure:**

| Quantity | Departs at (s) |
|---|--:|
| max\|Uc\| (all cells) | 21.751 |
| max\|Uc\| air (\|z\| < 0.29 band) | 21.784 |
| pitch rate | 22.151 |
| pitch | 22.234 |
| max\|Uc\| water | 22.241 |
| min segment tension, lines 3 and 4 | 22.261 |
| fairlead tension, lines 3 and 4 | 22.264 |

These never leave their band: Courant number (held at 0.5 by adjustTimeStep), interface Co, lines 1 and 2, alpha1/alphaI bounds (always [0, 1]), max\|alpha1 − alphaI\|, and dt. dt in H_S2 does fall from about 1e-3 to 9e-5 after 21.7 s, because the Co limit then acts on the air jet.

**The first departure is air Uc in the top cell row** (z = +0.492, directly below the atmosphere at z = 0.5) at x = +1.0 to +1.3: the coarse part outside the body zone, just inside the outlet relaxation zone.

| t (s) | Cell | x (m) | \|Uc\| (m/s) | Uc (x, z) (m/s) |
|--:|--:|--:|--:|---|
| 21.7 | 9483 | +1.098 | 13.6 | (−2, −14) |
| 21.8 | 9481 | +0.999 | 24.9 | (−14, −21) |
| 22.0 | 9485 | +1.196 | 31.7 | |
| 22.1 | 9486 | +1.245 | 69.7 | (−36, −60) |

The jet points down into the tank and occupies the top two rows (cells 9452 and 9455 are the row below), drifting with the wave crest. H_S's top row stays at 3–4 m/s.

**History of max air Uc per second.** H_S2: 5–6 m/s up to 19 s, 9.8 m/s in [20, 21), 23.9 m/s in [21, 22), then blow-up. H_S had a comparable excursion of 14.6 m/s in [19, 20) and recovered.

**After the jet.** The water side follows: max\|U\| in water reaches 0.88–1.47 m/s at x ≈ +1.49, z ≈ +0.075 (interface, α ≈ 0.55, at the zone edge) from 22.1 s. Lines 3 and 4 go slack (minimum segment tension ≤ 0) at 21.76 s in H_S2 and at 21.67 s in H_S, so slackness happens in both runs and is not the trigger. Lines 1 and 2 first go slack at 22.26 s, only in H_S2, after pitch is already past −6°.

**Pitch harmonics, 13–20.97 s.** Amplitude in degrees, with the change relative to plain fine in brackets:

| Run | h1 | h2 | h3 |
|---|---|---|---|
| plain fine | 2.430 | 3.206 | 0.422 |
| H_S | 2.784 (+14.6 %, −11.6°) | 2.665 (−16.9 %, −21.0°) | 0.755 (+79 %, −18.7°) |
| H_S2 | 2.249 (−7.4 %, +1.4°) | 2.396 (−25.3 %, −17.1°) | 0.599 (+42 %, −14.0°) |

The switch fixes the h1 phase but leaves h2 25 % low and 17° early.

## H_S2 restart 21.4–21.8 s: terms of the Uc equation in the air jet

Restart case `zonefix/hybrid/H_S2r` runs from the 21.4 s write of H_S2 to 21.8 s with the same switches. It uses solver build `private_uc`, which is private_fu plus the log-only switch `ucCellLog` ([swenseUcCellLog.H](../src/swenseWaveTank/solver/swenseUcCellLog.H), hooked in at the end of the last outer corrector in [swenseTimeLoop.H](../src/swenseWaveTank/solver/swenseTimeLoop.H)). The analysis script is [s9.py](s9.py).

The logged cells are the top two rows (z = 0.477 and 0.492) at x 0.9–1.4: 22 cells in total. Each term is taken on the final fields of the step, in N/m³, with x and z components. The terms are ddt, convection, viscous, the incident mismatch brackets of incTerm, the cell closure (1 − lam) r_I, the face term −R[(1 − lamf) G_f], −grad pc, the gravity term −R[gh_f snGrad(Drho)], and the cell residual R. The log also records lam, w, rho and rhoI, plus pc_b, Uc_b, phiComp and snGrad pc on the 11 atmosphere faces above the cells.

**Reproduction is close but not bit for bit.** On restart the solver rebuilds phi from the written fields, so the first max Co is 0.531 against 0.513 in H_S2. The time steps therefore differ from the first step on (21.401053 against 21.401087). MoorDyn also restarts from equilibrium (`MoorDyn_Init`).

| Quantity | Max \|H_S2r − H_S2\|, 21.4–21.6 s | Max \|H_S2r − H_S2\|, 21.4–21.8 s |
|---|--:|--:|
| pitch | 0.22° | 0.22° |
| pitch rate | 3.0°/s | 3.0°/s |
| max\|Uc\| (24.9 m/s at 21.8 s in both runs) | 0.19 m/s | 0.19 m/s |
| air max\|Uc\| | 1.3 m/s | 1.3 m/s |
| fairlead tension | up to 0.56 N | up to 0.56 N |

The jet grows exactly as in H_S2.

**What grows first.** Each quantity is the maximum over the 22 cells. The growth time is the first step at which it reaches 1.5× (and 2×) its value over 21.40–21.41 s and stays there for 0.01 s.

| Quantity | ×1.5 at (s) | ×2 at (s) | Value at 21.40 s | Max to 21.8 s | Cell |
|---|--:|--:|--:|--:|--:|
| **convection div(rhoPhi, Uc)** | **21.468** | **21.526** | 60 | 1690 | 9453 (x 1.147, z 0.477) |
| incident mismatch (inc) | 21.476 | 21.541 | 38 | 629 | 9482 |
| ddt(rho, Uc) | 21.497 | 21.532 | 25 | 949 | 9486 |
| cell residual R | 21.502 | 21.580 | 13 | 397 | 9484 |
| atmosphere inflow −phiComp/\|Sf\| (m/s) | 21.518 | 21.598 | 4.2 | 21.6 | top faces |
| \|Uc\| (m/s) | 21.518 | 21.612 | 4.5 | 24.9 | 9483 |
| \|Drho\|, alphaI | 21.523 | 21.599 | 0.47 kg/m³, 4.7e-4 | 1.71 kg/m³, 1.7e-3 | 9482 |
| gravity term | 21.537 | 21.621 | 101 | 325 | 9482 |
| −grad pc | 21.641 | 21.684 | 151 | 852 | 9452 |

These terms stay negligible: the closure (1 − lam) r_I is zero to 21.6 s (max 38 later), the face term peaks at 0.43, and viscous stays below 0.15.

**Readings:**
- **The first term to grow is convection, div(rhoPhi, Uc), at 21.47 s in cell 9453, the second row from the top.** The incident mismatch bracket follows 0.01 s later in cell 9482 (top row), then ddt.
- **The air is not clean in the incident field.** These are pure air cells (alpha1 ≤ 4e-6), but alphaI = 3.5e-4 to 1.7e-3 and rises through the window. That makes rhoI = 1.35 → 2.69 against rho = 1.00, so the gravity term −gh_f snGrad(Drho) is 70–330 N/m³ on air of density 1. The incident mismatch brackets, (rho − rhoI) times the incident ddt and convection, are non-zero for the same reason.
- **The jet enters through the top boundary.** The atmosphere boundary has pc = 0 (fixedValue) and zeroGradient Uc. Over 0.05 s windows, the largest inflow through the top faces rises 5.1 → 8.5 m/s (21.40–21.60) → 21.6 m/s (21.80), and snGrad pc at those faces rises 446 → 1860 Pa/m. In cell 9481 Uc_z changes from −0.3 to −20 m/s between 21.67 and 21.80 s, with convection as the largest LHS term (up to 1050 N/m³).
- w = 0 in all 22 cells for the whole run: these cells are outside the body zone, and at their x position the relaxation weight is still zero. lam = 1 to 21.6 s; from 21.614 s the ratio test (rhoI/rho > 2) sets lam = 0 in some top-row cells, which is when the closure term appears.

No fix.

## isoAdvector tolerances 1e-8 (A) and pc_b = -Drho_b gh_b on the atmosphere (B), 2026-11-24

**Changes.**
- A: fvSolution `"(alpha.water|alphaI).*"` isoFaceTol 1e-8, surfCellTol 1e-8 (was 1e-2 in H_S2 and H_P).
- B: `pcAtmosphereConsistent` (waveProperties Switch, default off, read in swenseSetup.H like incidentHydrostatic). In swensePcEqn.H, before constrainPressure, the atmosphere value of pc is set to -Drho_b gh_b, so p = pInc + pc + Drho gh is the same on that patch whatever alphaI is there. pc stays fixedValue.
- Build private_pa. Runs: H_S3 = H_S2 + A + B; H_P2 = H_P + A (plain waveIsoFoam).

**Gates.**
- 1, switch off, 0.3 s (G1_off): all 66 files in 0.1, 0.2 and 0.3 are byte-identical to H_S2.
- 2, pure waves without a body, switch on, 200 steps (G2_on): max|alpha1 - alphaI| = 0, max|Uc| = 0, max|pc| = 0 at every step; pc on the atmosphere is uniform 0. The gate is weak because the 12 s ramp is tiny within 0.2 s.

**Results, 13-23 s** (s10.py, s6.py; H_S2 crashed at 22.69 s, so its values are to 22.6 s).

| | H_S2 | H_S3 | H_P | H_P2 | target |
|---|--:|--:|--:|--:|--:|
| max alphaI, cells wholly above etaI + 2 dz, writes 12-23 s | 5.5e-3 | 2.5e-9 (all 111 writes ≤ 1e-6) | - | - | ≤ 1e-6 |
| max alpha1, same cells | 0.99 (crash) | 1.6e-4 at x ≈ 0, z ≈ 0.06 (next to the body); 1.2e-9 for z > 0.3 | 7.0e-3 | 1.1e-9 | ≤ 1e-6 |
| max air \|Uc\| (m/s) | 748e81 (crash), 6.5 to 20 s | 10.7 at 19.39 s; 6-10 per second after 9 s | - | - | < 5 |
| wall time | 6317 s to 22.69 s | 4733 s | 10849 s | 14341 s | |
| RAO surge / heave / pitch / tension vs plain fine | H_S: +4.9 / +1.3 / +14.9 / +16.4 % | +0.6 / -0.8 / **+11.5** / **+6.5** % | +1.1 / +0.2 / +1.8 / +1.5 % | +1.2 / -0.3 / **+4.2** / +1.5 % | 3 % |
| RAO phase pitch / tension vs plain fine | H_S: -5.4 / +3.1° | +10.3 / +6.8° | +4.5 / +1.4° | +4.4 / +1.7° | |
| mean surge mm / tension N / pitch deg / heave mm (plain fine +8.10 / 0.2324 / 0.673 / +1.20) | H_S: +9.62 / 0.2758 / 1.128 / -0.37 | +7.58 / 0.2338 / **1.038** / **-0.72** | +8.08 / 0.2300 / 0.636 / +0.76 | +7.49 / 0.2276 / 0.748 / +0.77 | 3 % |
| pitch h1 / h2 (plain fine 2.430 / 3.211 deg) | H_S: +17.3 % -10.6° / -17.8 % -23.1° | +12.0 % +1.7° / -14.9 % -23.5° | +0.2 % / -2.3 % | +2.8 % +4.2° / -1.4 % -3.5° | |
| wave A1 at the box gauge / x = -1.25 (ratio to target; plain fine 0.961 / 0.997) | H_S: 0.981 / 1.097 | 0.966 / 1.060 | 0.945 / 1.002 | 0.948 / 1.004 | |

**Readings.**
- A removes the stray incident alpha in the air: max alphaI above etaI + 2 dz drops from 5.5e-3 to 2.5e-9. In H_P2 the same tolerance takes stray alpha1 in air from 7e-3 to 1e-9.
- H_S3 runs to 23 s without the 21.5 s jet. The air Uc is still above target, though: per-second maxima are 6-10.7 m/s against 5-6.5 m/s in H_S2 before its crash. B (together with A) raises the air velocity in the first 20 s.
- H_S3 is closer to plain fine than H_S on surge, heave, tension and the upstream wave height. The pitch RAO is still +11.5 % and the mean pitch +0.37 deg too high. The pitch h2 is 15 % low and lags 23.5°, unchanged from H_S. So the pitch error does not come from the stray alphaI.
- The 1e-8 tolerance alone moves plain waveIsoFoam (H_P2 vs H_P) by up to +2.4 points on the pitch RAO and +0.11 deg on the mean pitch, and it makes the run 32 % slower.
