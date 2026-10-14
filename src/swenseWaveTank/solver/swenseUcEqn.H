// Complementary velocity (document Rev1, eq. 16.4), conservative form:
//
//   d(rho uC)/dt + div(m uC) + div(D uI) - (div D) uI - S uC
//       = -grad(pC_hat) - pI grad(rho/rhoI) + div(mu(grad uC + grad uC^T))
//
// with S = d(rho)/dt + div(m) from the discrete VOF transport (document
// Rev2, 16.7 and its mass-source paragraph). The term -S uC is implicit.
//
// with m = rhoPhi (total mass flux, from the alpha transport) and
// D = rhoPhi - (rho/rhoI)_f rhoPhiI (eq. 16.5, from the alphaI transport).
// The interface force is a face flux (eq. 16.6), used in the predictor here
// and in the pressure equation (swensePcEqn.H). No cell-centred force.
//
// Air damping and forward-speed terms are absent (both zero in the reference
// case), so they are not included.

    // Mass flux D of the incident convection (eq. 16.5)
    const surfaceScalarField D
    (
        "D",
        rhoPhi - fvc::interpolate(rhoByRhoI)*rhoPhiI
    );

    // D grad(uI) written as div(D uI) - (div D) uI (eq. 16.4)
    const volVectorField DdotGradUI
    (
        fvc::div(D, incident.UInc(), "div(D,UInc)")
      - fvc::div(D)*incident.UInc()
    );

    fvVectorMatrix UcEqn
    (
        fvm::ddt(rho, Uc)
      + fvm::div(rhoPhi, Uc)
      + DdotGradUI
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), Uc)
      + turbulence->divDevRhoReff(rho, Uc)
    );

    UcEqn.relax();

    // Interface force per unit area divided by rA (eq. 16.6, without the
    // (rA)_f factor, which the pressure corrector applies):
    //     phiSource0 = -pI,f snGrad(rho/rhoI) |Sf|
    const surfaceScalarField phiSource0
    (
        "phiSource0",
        -incident.pfInc()*fvc::snGrad(rhoByRhoI)*mesh.magSf()
    );

    if (pimple.momentumPredictor())
    {
        solve
        (
            UcEqn
         ==
            fvc::reconstruct
            (
                phiSource0 - fvc::snGrad(pc)*mesh.magSf()
            )
        );
    }
