// Incident pressure balance, V1 (momentum_swense Part V): face flux of the
// incident residual r_I after the last pressure corrector of the step.
// Included from swensePcEqn.H when incidentPressureBalance is on; uses rAU and
// rAUf of that corrector, and the operators of swenseUcEqn.H (UI, muOf,
// viscUI, convScheme are in scope there).
//
//   r_I cell part  = ddt(rhoI UI) + div(rhoPhiI UI) - V[mu(alphaI); UI]
//                    (the incident inertia, convection and viscous terms of
//                    the closure, swenseUcEqn.H)
//   G_f            = ghf snGrad(rhoI) + snGrad(pInc - rhoI gh)
//   Phi_f[r_I]     = (interpolate(rAU r_I cell) . Sf) + rAUf G_f |Sf|
//
// The cell part goes through the interpolation that builds phiHbyA from
// rAU H (swensePcEqn.H: interpolate(HbyA) . Sf).

const volVectorField rIBalanceCell
(
    "rIBalanceCell",
    fvc::ddt(rhoI, UI)
  + fvc::div(rhoPhiI, UI, convScheme)
  - viscUI(muOf(alphaI)())
);

const surfaceScalarField GIBalance
(
    "GIBalance",
    ghf*fvc::snGrad(rhoI) + fvc::snGrad(incident.pInc() - rhoI*gh)
);

const surfaceScalarField PhiRI
(
    "PhiRI",
    (fvc::interpolate(rAU*rIBalanceCell) & mesh.Sf())
  + rAUf*GIBalance*mesh.magSf()
);
