// Complementary velocity: momentum with the incident convection term.
// Air damping and forward-speed terms are absent (both zero in the reference
// case), so they are not included.
    surfaceScalarField phiCompAbs
    (
        "phiCompAbsolute",
        Ucf & mesh.Sf()
    );

    volVectorField fvcDivPhicUInc
    (
        fvc::div(phiCompAbs, incident.UInc(), "div(phiCompAbsolute,UInc)")
    );

    fvVectorMatrix UcEqn
    (
        fvm::ddt(rho, Uc)
      + fvm::div(rhoPhi, Uc)
      + rho*fvcDivPhicUInc
      + turbulence->divDevRhoReff(rho, Uc)
    );

    UcEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UcEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(pc)
                  - fvc::snGrad(rho)*incident.pfInc()/rho1
                )*mesh.magSf()
            )
        );
    }
