# Momentum blend: change points (branch swense-momentum-blend, from SWENSE 04a7748)

Status: notes only. No code changed. Not committed.

Definition (checked symbolically, 1D, Euler + alpha transported by uI):
- W(alphaI) = alphaI rhoW / rhoI,  rhoI = alphaI rhoW + (1 - alphaI) rhoA
- uI = W uW + (1 - W) uA
- D_I uI = -grad(pI)/rhoI + g + F_I,  F_I = (pW - pA) grad(alphaI)/rhoI - W(1-W) dU.grad dU
- div uI = q_I = (rhoW rhoA / rhoI^2) dU.grad(alphaI)
- d_t rhoI + div(rhoI uI) = rhoI q_I (incident mass flux = rhoI uI exactly)

Change points in the code:
1. swenseIncident.C:486  weight W: currently clamp(alphaI). Change to W(alphaI) above.
2. swenseIncident.C blendUP (~520-580): velocity blend uses weights from weights(). Velocity weights become W and 1-W. Pressure part unchanged (one-sided, z - eta).
3. swenseIncident.C weights() (~480-510): drift counting uses W, so it follows (1).
4. swenseSetup.H:264: incident mass flux = interp(rhoI) * phiInc. Stays; becomes exact incident momentum flux.
5. swenseTimeLoop.H:26-57: projection of phiInc. Divergence is now q_I (small). Projection still runs.
6. swenseUcEqn.H:46-49: phiSource0 = -psiWf snGrad(rho - rhoI) |Sf|. Unchanged in form. Its incident residual is now only the W(1-W) dU.grad dU term.
7. swenseAlphaRelax.H:15, swenseSetup.H:147: rhoI from alphaI. Unchanged.

Checks before any run:
- No-body invariance: uC = 0 to round-off (paper Thm 1 / Prop 1).
- Pure wave: uC stays 0.
- Band residual: compare integrated F_I with current formulation (should match in the limit).

Open: whether this changes the box Fx gap at delta = 0.024 m. Needs one run, on your go.

## Update (implemented, switch off by default)
- waveProperties `momentumBlend` (default false): W = a rhoW/rhoI in weights(); with true the pressure is the linear absolute blend pI = a pW + (1-a) pA (blendUP), replacing the one-sided scratch pressure.
- Default false = previous scheme (W = alphaI). Commit a0b1ba0 had changed the default; this commit restores it behind the switch.
- Needs case keys (not in code): Air2Air `pressureModel bernoulli` for the full pA.
- Not changed: phiSource0 psi form (finding 1 and 3 of the theory review).
- Not compiled, not run.
