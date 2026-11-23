// Pressure equation for pc = p - pInc - Drho gh (momentum_swense_Rev1,
// eq. 34). It constrains the complementary flux alone:
//   sum_f phiComp_f = 0,
//   phiComp_f = phiH_f - rAU_f (gh_f snGrad(Drho) + snGrad(pc)) |Sf|.
// The incident flux is not in it (M6).
{
    volScalarField rAU("rAU", 1.0/UcEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UcEqn.H(), Uc, pc));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rAU*rho)
       *fv::ddtScheme<vector>::New
        (
            mesh,
            mesh.ddtScheme("ddtUcFlow")
        ).ref().fvcDdtUfCorr(Uc, Ucf)
    );

    if (pc.needReference())
    {
        adjustPhi(phiHbyA, Uc, pc);
    }

    // Restoring force of the displaced interface, -gh_f snGrad(Drho)
    surfaceScalarField phig
    (
        "phig",
        -ghf*fvc::snGrad(Drho)*rAUf*mesh.magSf() - faceLim*rAUf
    );

    phiHbyA += phig;

    constrainPressure(pc, Uc, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcEqn
        (
            fvm::laplacian(rAUf, pc) == fvc::div(phiHbyA)
        );

        if (pc.needReference())
        {
            pcEqn.setReference(pcRefCell, getRefCellValue(pc, pcRefCell));
        }

        if
        (
            (pimple.corrPISO() == 1 || pimple.corrPISO() == pimple.nCorrPISO())
         && pimple.finalNonOrthogonalIter()
        )
        {
            pcEqn.solve(mesh.solverDict("pcFinal"));
        }
        else
        {
            pcEqn.solve();
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phiComp = phiHbyA - pcEqn.flux();
            Uc = HbyA + rAU*fvc::reconstruct((phig - pcEqn.flux())/rAUf);
            Uc.correctBoundaryConditions();
        }
    }

    Ucf = fvc::interpolate(Uc);
    {
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Ucf += n*(phiComp/mesh.magSf() - (n & Ucf));
    }

    // Face coefficient for the flux correction after the next mesh motion
    // (swenseCorrectPhi.H); read nowhere else, so the static path is unchanged
    rAUfPrev = rAUf;

    // Incident pressure balance (incidentPressureBalance, default off;
    // momentum_swense Part V): after the last corrector of the step, the flux
    // of r_I (V1), the correction deltaHeld (V2) and the log (V5)
    const bool balanceNow =
        incidentPressureBalance
     && pimple.finalIter()
     && pimple.corrPISO() == pimple.nCorrPISO();

    if (balanceNow)
    {
        #include "swenseBalanceFlux.H"
        #include "swenseBalanceSolve.H"
        #include "swenseBalanceLog.H"
    }

    // Total fields (Section 12.4)
    phi == phiComp + phiInc - meshPhi;
    #include "continuityErrs.H"

    U == Uc + incident.UInc();

    if (balanceNow)
    {
        // V3: p with the held incident pressure pIncHeld = pInc - deltaHeld;
        // forces, sixDoF and output read this p
        p == (hydroDiscrete ? pIncD : incident.pInc()) - deltaHeld + pc + Drho*gh;
    }
    else
    {
        p == (hydroDiscrete ? pIncD : incident.pInc()) + pc + Drho*gh;
    }

    p_rgh == p - rho*gh;
}
