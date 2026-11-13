// Complementary momentum equation (momentum_swense_Rev1, eq. 32 and 33):
//
//   ddt(rho Uc) + C[F; Uc] - V[mu(alpha); Uc]
//       = -I + R[(-gh_f snGrad(Drho) - snGrad(pc)) |Sf|]
//
// with F = rhoPhi (mass flux of the alpha1 transport) and the explicit
// incident term
//
//   I = [ddt(rho UI) - ddt(rhoI UI)]
//     + [C[F; UI] - C[FI; UI]]
//     - [V[mu(alpha); UI] - V[mu(alphaI); UI]]
//
// UI is the corrected incident velocity (step 1d), FI = rhoPhiI. Every
// bracket is one operator applied twice and subtracted (M5): it is zero bit
// for bit when alpha1 = alphaI and rhoPhi = rhoPhiI. No incident term stands
// on its own. The incident pressure does not appear.

    const volVectorField& UI = incident.UInc();

    // Molecular viscosity of a volume fraction, one function for both
    auto muOf = [&](const volScalarField& s) -> tmp<volScalarField>
    {
        const volScalarField lim(min(max(s, scalar(0)), scalar(1)));

        return
            lim*rho1*mixture.nuModel1().nu()
          + (scalar(1) - lim)*rho2*mixture.nuModel2().nu();
    };

    // Explicit viscous operator on UI for a viscosity field (the form of the
    // solver's laminar stress: laplacian plus the transpose part)
    const volTensorField gradUI("gradUI", fvc::grad(UI));

    auto viscUI = [&](const volScalarField& mu) -> tmp<volVectorField>
    {
        return
            fvc::laplacian(mu, UI, "laplacian(mu,UInc)")
          + fvc::div(mu*dev2(T(gradUI)), "div(mu*dev2(T(grad(UInc))))");
    };

    // Convection of UI with the scheme of Uc, for the total and the incident
    // mass flux
    const word convScheme("div(rhoPhi,Uc)");

    const volVectorField incTerm
    (
        "incTerm",
        (fvc::ddt(rho, UI) - fvc::ddt(rhoI, UI))
      + (fvc::div(rhoPhi, UI, convScheme) - fvc::div(rhoPhiI, UI, convScheme))
      - (viscUI(muOf(alpha1)()) - viscUI(muOf(alphaI)()))
    );

    // Logging only (not used by any solve): residual of the incident momentum
    //   rI = [ddt(rhoI UI) + div(rhoPhiI UI)]
    //      + R[(snGrad(pInc - rhoI gh) + ghf snGrad(rhoI)) |Sf|]
    // written with the time directory (AUTO_WRITE)
    const volVectorField rI
    (
        IOobject
        (
            "rI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::ddt(rhoI, UI) + fvc::div(rhoPhiI, UI, convScheme))
      + fvc::reconstruct
        (
            (fvc::snGrad(incident.pInc() - rhoI*gh) + ghf*fvc::snGrad(rhoI))
           *mesh.magSf()
        )
    );

    // A local object is not written by the time loop: write it here
    rI.write();

    fvVectorMatrix UcEqn
    (
        fvm::ddt(rho, Uc)
      + fvm::div(rhoPhi, Uc)
      + turbulence->divDevRhoReff(rho, Uc)
     ==
      - incTerm
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
                (-ghf*fvc::snGrad(Drho) - fvc::snGrad(pc))*mesh.magSf()
            )
        );
    }
