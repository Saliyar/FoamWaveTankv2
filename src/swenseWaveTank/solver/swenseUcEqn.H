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

    // Switch b: incident time derivative as a plain Euler difference of the
    // cell content, without the moving-mesh volume factor V0/V
    auto ddtTerm = [&]() -> tmp<volVectorField>
    {
        if (!swenseNoMeshDdtUc)
        {
            return fvc::ddt(rho, UI) - fvc::ddt(rhoI, UI);
        }

        return tmp<volVectorField>
        (
            new volVectorField
            (
                "ddtPlain",
                ((rho - rhoI)*UI - (rho.oldTime() - rhoI.oldTime())*UI.oldTime())
               /mesh.time().deltaT()
            )
        );
    };

    volVectorField incTerm
    (
        "incTerm",
        ddtTerm()
      + (fvc::div(rhoPhi, UI, convScheme) - fvc::div(rhoPhiI, UI, convScheme))
      - (viscUI(muOf(alpha1)()) - viscUI(muOf(alphaI)()))
    );

    // Limited closure (waveProperties limitedClosure, default off).
    // lam = 0 in cells where max(rho,rhoI)/min(rho,rhoI) > 2, else 1;
    // lamf = min of owner and neighbour (boundary: lam of the owner cell).
    // With lam = 0 the incident residual is removed from the incident term
    // and from the face flux, so U = UI + Uc obeys the plain momentum
    // equation (p_rgh = pInc + pc - rhoI gh) with the operator of UI.
    volScalarField lam
    (
        IOobject("lam", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("lam", dimless, 1.0)
    );

    surfaceScalarField lamf
    (
        IOobject("lamf", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("lamf", dimless, 1.0)
    );

    if (limitedClosure)
    {
        label nZero = 0;

        forAll(lam, c)
        {
            const scalar rA = rho[c];
            const scalar rB = rhoI[c];

            if (limitedClosureContinuous)
            {
                // Continuous closure, linear in the mismatched mass:
                // 1 - lam = |rho - rhoI|/(rho1 - rho2) = |alpha1 - alphaI|.
                // nZero counts the cells with lam < 0.5
                lam[c] = 1.0 - min
                (
                    mag
                    (
                        min(max(alpha1[c], scalar(0)), scalar(1))
                      - min(max(alphaI[c], scalar(0)), scalar(1))
                    ),
                    scalar(1)
                );

                if (lam[c] < 0.5)
                {
                    ++nZero;
                }
            }
            else if (max(rA, rB)/min(rA, rB) > limitedClosureRatio)
            {
                lam[c] = 0.0;
                ++nZero;
            }
        }

        // Plain closure in the body zone (plainBodyZone): lam = 0 in every
        // zone cell; lam_f = 0 follows from the min below
        if (plainBodyZoneID >= 0)
        {
            for (const label c : mesh.cellZones()[plainBodyZoneID])
            {
                // nZero counts cells with lam < 0.5 (both closure modes)
                if (lam[c] >= 0.5)
                {
                    ++nZero;
                }

                lam[c] = 0.0;
            }
        }

        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            lamf[f] = min(lam[mesh.owner()[f]], lam[mesh.neighbour()[f]]);
        }

        // Boundary lam_f = lam of the owner cell (the patch value of lam is
        // never set and stayed 1)
        forAll(lamf.boundaryField(), patchi)
        {
            lamf.boundaryFieldRef()[patchi] =
                lam.boundaryField()[patchi].patchInternalField();
        }

        Info<< "SWENSE limitedClosure: lam = 0 cells " << nZero << endl;

        incTerm += (1.0 - lam)*
        (
            fvc::ddt(rhoI, UI) + fvc::div(rhoPhiI, UI, convScheme)
          - viscUI(muOf(alphaI)())
        );
    }

    // Logging only (not used by any solve): residual of the incident momentum
    //   rI = [ddt(rhoI UI) + div(rhoPhiI UI)]
    //      + R[(snGrad(pInc - rhoI gh) + ghf snGrad(rhoI)) |Sf|]
    // Computed and written only when waveProperties logIncidentResidual is on
    // (default off), and only at normal write times.
    if (logIncidentResidual && runTime.writeTime())
    {
        const volVectorField rI
        (
            IOobject
            (
                "rI",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (fvc::ddt(rhoI, UI) + fvc::div(rhoPhiI, UI, convScheme))
          + fvc::reconstruct
            (
                (fvc::snGrad(incident.pInc() - rhoI*gh) + ghf*fvc::snGrad(rhoI))
               *mesh.magSf()
            )
        );

        rI.write();
    }

    fvVectorMatrix UcEqn
    (
        fvm::ddt(rho, Uc)
      + fvm::div(rhoPhi, Uc)
      + turbulence->divDevRhoReff(rho, Uc)
     ==
      - incTerm
    );

    UcEqn.relax();

    // Face term of the limited closure; exactly zero when limitedClosure is
    // off, because lamf = 1 then. V4 of the incident pressure balance
    // (incidentPressureBalance on): r_I is built with the held incident
    // pressure pIncHeld = pInc - deltaHeld, deltaHeld of the previous step
    // (zero before the first balance). Off: pInc itself, unchanged.
    tmp<volScalarField> tpIncClosure(incident.pInc());

    if (incidentPressureBalance)
    {
        tpIncClosure = incident.pInc() - deltaHeld;
    }

    surfaceScalarField incFace
    (
        "incFace",
        fvc::snGrad(tpIncClosure() - rhoI*gh) + ghf*fvc::snGrad(rhoI)
    );

    // limitedClosureJump (default off): remove the jump of the two wave
    // models at the alphaI interface, snGrad(gamma) [mean(P_W - P_A)
    // + (rho1 - rho2) ghf] (swenseSetup.H)
    if (limitedClosureJump)
    {
        const volScalarField gammaJ
        (
            "gammaJ",
            min(max(alphaI, scalar(0)), scalar(1))
        );

        const volScalarField PJ
        (
            "PJ",
            incident.phaseJumpP() - (rho1 - rho2)*gh
        );

        incFace -=
            fvc::snGrad(gammaJ)
           *(linearInterpolate(PJ) + (rho1 - rho2)*ghf);
    }

    const surfaceScalarField faceLim
    (
        (1.0 - lamf)*incFace*mesh.magSf()
    );

    if (pimple.momentumPredictor())
    {
        solve
        (
            UcEqn
         ==
            fvc::reconstruct
            (
                (-ghf*fvc::snGrad(Drho) - fvc::snGrad(pc))*mesh.magSf()
              - faceLim
            )
        );
    }
