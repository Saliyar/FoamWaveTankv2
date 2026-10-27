# Step 0: ramp map for the SWENSE incident fields (Tsoft > 0)

Status: map only, no code change, no run. Branch SWENSE_dynamic (worktree wt_SWENSE_dynamic).

## Why Tsoft is rejected today

- The FatalError is `src/swenseWaveTank/swenseIncident.C:191-196`, "swenseWaveFoam requires Grid2GridCoeffs Tsoft = 0".
- Introduced in commit `5887a211` (Saliyar, 2026-10-14, "Add SWENSE wave library"), comment "The SWENSE path needs the dynamic pressure without the soft start". No other rationale is recorded in the commit.
- `doc/swense_port/02_spec.md:67` says no ramp factor is applied in the SWENSE incident path (a deliberate study choice, B1.4-B1.5).
- `doc/swense_port/02_spec.md:145, 157` and `03_build_notes.md:88`: Grid2Grid `pExcess` carries a soft start, so the port requires `Tsoft = 0` or bypass.
- The analytic path (`waveTheoryForcing.C:67-78`) already requires the theory's own Tsoft to be 0, so an incident-level ramp would not double-ramp it.

## Quantities

| # | Quantity | Source (file:line) | Consumers | Ramp f(t)? | Time derivative needed? |
|---|---|---|---|---|---|
| 1 | water eta, Grid2Grid | `swenseIncident.C:472` (etaAt), `:679` (etaVal_), `swenseSurface.C:47`, `swenseWaveAlpha...C:102` | classification, blend, alphaH target, wetted geometry, alpha BC | yes: f*eta | eta_t: f' eta + f eta_t (not consumed) |
| 2 | water eta, waveTheoryForcing | same sites, theory `eta` | same | yes | same |
| 3 | water U, Grid2Grid | `swenseIncident.C:420` (waterU), `:670` (UPt_ via evaluatePoints) | blend, velocity BC | yes: f*U | none |
| 4 | water U, theory | same, theory `U` (`waveTheoryForcing.C:125`) | blend, `swenseWaveVelocity...C:128` | yes: f*U | none |
| 5 | water pExcess, Grid2Grid | `swenseIncident.C:434` (waterP), `:670` (pExPt_) | pWs_, interface pressure | **NO, see stop** | none |
| 6 | water pExcess, theory | `waveTheoryForcing.C:115` (returns p/rho + g z) | same | **NO, see stop** | none |
| 7 | air U, Pd (Air2Air) | `swenseIncident.C:457` (flowTables), `:715` (flow) | blend, pAs_ | yes: f*U_air, f*Pd; keep -rho_a g z unramped | Pd_t, U_air_t: f' q + f q_t (not consumed) |
| 8 | Air2Air coefficients A, A_t | `Air2AirEnhanced.C` (update at `swenseIncident.C:664`), `analyticAirCoeffs.H` | air outputs only | NOT on coefficients: A is nonlinear in eta (order 2, 3), so f*eta does not give f*A. Ramp applied to outputs (7) | A_t: record difference (`Air2AirEnhanced.H:20`), used inside Pd only |
| 9 | eta_t inputs to Air2Air | `modeFSt` (modes file) and `analyticAirCoeffs.H` (theory path) | coefficients | not ramped (coefficients unramped, outputs ramped) | used only inside coefficients |
| 10 | blend pressure references pWs_, pAs_ | `swenseIncident.C:704-716` | interface and blend pressure | consistent if 5-7 are ramped (computed from the same functions) | none |
| 11 | analytic ddt of UInc or pInc | none found in solver (only `alphaI.oldTime()` in `swenseBudgetLog.H:55`) | - | - | none consumed |
| 12 | relaxation-zone targets | alphaH = blendAlpha(etaAt) (`swenseIncident.C:583-616`); Uc -> 0 | relaxation | follows 1-2 | none |

## Stop

Items 5 and 6 cannot take f(t) consistently as a plain multiplier:

- Theory path (item 6): `waveTheoryForcing::pExcess` returns `p/rho + g z` (`waveTheoryForcing.C:115`, comment "swenseIncident wants p/rho + g z"). The hydrostatic term is inside pExcess, and `waterP` removes it with `- magG_*zw`. Multiplying pExcess by f would ramp the hydrostatic term, giving a spurious -(1-f) g z pressure at start-up.
- Grid2Grid path (item 5): the split of pExcess into dynamic and hydrostatic parts is the open inconsistency already recorded in `doc/swense_analytic/00_facts.md` section 3 (Grid2Grid returns Pa, swenseIncident treats it as per unit density). The ramp rule cannot be applied consistently until that is measured.

A correct rule needs the dynamic part of pExcess per forcing: for the theory path it is `pExcess - g z` (exact), for Grid2Grid it is open. Stopping here per the step-0 rule. No implementation (step 1) and no runs (step 2) done.

## Ramp law (for reference)

waves2Foam `Grid2Grid.C:126-128`: f(t) = sin(2 pi / (4 Tsoft) min(Tsoft, t)) = sin(pi/2 min(t, Tsoft)/Tsoft); f' = (pi/(2 Tsoft)) cos(...) for t < Tsoft, 0 after.
