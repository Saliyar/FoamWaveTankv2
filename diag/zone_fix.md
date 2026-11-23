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
