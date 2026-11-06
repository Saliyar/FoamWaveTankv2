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

    // Interface force per unit area divided by rA, without the (rA)_f factor
    // that the pressure corrector applies. Interface force (rho - rhoI) grad(psi),
    // psi = (pW - pW at eta)/rho_w, in the interFoam balanced form (see theory doc Sec 16):
    //     phiSource0 = -psi_f snGrad(rho - rhoI) |Sf|
    const surfaceScalarField phiSource0
    (
        "phiSource0",
        -incident.psiWf()*fvc::snGrad(rho - rhoI)*mesh.magSf()
    );

    // Blend residual on the right-hand side of the Uc equation, theory eq. (compcons):
    //   right side = -grad(pC_hat) - pI grad(rho/rhoI) - rho F   (sigma = 0, advective form)
    // with F = F_I (eq. incmom) on the momentum blend and F = F^v (eq. Fv) on the
    // default velocity blend. In both cases:
    //   face part  : -(rho/rhoI)_f dP_f snGrad(alphaI) |Sf|,  dP = pW - pA (absolute)
    //   cell part  : momentum : +rho alphaI(1-alphaI) rho1 rho2/rhoI^2 dU.grad(dU)
    //                velocity : +rho alphaI(1-alphaI) [ (rho1-rho2)/rhoI (aW-aA) + dU.grad(dU) ]
    // aW - aA = -grad(pW)/rho1 + grad(pA)/rho2. The psi-form phiSource0 is kept (note Sec. 11).
    const volScalarField& pWc = incident.pWc();
    const volScalarField& pAc = incident.pAc();
    const volVectorField& dUc = incident.dUc();

    surfaceScalarField dPface("dPface", fvc::interpolate(pWc - pAc));
    if (incident.momentumBlend())
    {
        dPface = incident.dPf();
    }

    const surfaceScalarField blendPhi
    (
        "blendPhi",
        -fvc::interpolate(rho/rhoI)*dPface*fvc::snGrad(alphaI)*mesh.magSf()
    );

    const volVectorField dUgraddU("dUgraddU", dUc & fvc::grad(dUc));

    volVectorField blendCell("blendCell", rho*alphaI*(1.0 - alphaI)*dUgraddU);
    if (incident.momentumBlend())
    {
        blendCell =
            rho*alphaI*(1.0 - alphaI)*rho1*rho2/sqr(rhoI)*dUgraddU;
    }
    else
    {
        const volVectorField aDiff
        (
            "aDiff",
            -fvc::grad(pWc)/rho1 + fvc::grad(pAc)/rho2
        );
        blendCell +=
            rho*alphaI*(1.0 - alphaI)*((rho1 - rho2)/rhoI)*aDiff;
    }

    if (pimple.momentumPredictor())
    {
        solve
        (
            UcEqn
         ==
            fvc::reconstruct
            (
                phiSource0 + blendPhi - fvc::snGrad(pc)*mesh.magSf()
            )
          + blendCell
        );
    }
