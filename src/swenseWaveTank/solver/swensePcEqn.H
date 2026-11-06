// Pressure corrector for the complementary pressure pc (SWENSE form)
{
    volScalarField rAU("rAU", 1.0/UcEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    // Kept for the flux correction after a mesh motion (swenseCorrectPhi.H)
    rAUfPrev = rAUf;
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

    // Blend residual as face fluxes times rAU: the same blendPhi and blendCell as
    // the Uc predictor (swenseUcEqn.H), the cell part through its face interpolate
    surfaceScalarField phiSource(phiSource0*rAUf);
    phiSource += (blendPhi + (fvc::interpolate(blendCell) & mesh.Sf()))*rAUf;

    phiHbyA += phiSource;

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
            if (includeDivUInc)
            {
                solve(pcEqn == fvc::div(phiInc), mesh.solverDict("pcFinal"));
            }
            else
            {
                pcEqn.solve(mesh.solverDict("pcFinal"));
            }
        }
        else
        {
            if (includeDivUInc)
            {
                solve(pcEqn == fvc::div(phiInc));
            }
            else
            {
                pcEqn.solve();
            }
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phiComp = phiHbyA - pcEqn.flux();
            Uc = HbyA + rAU*fvc::reconstruct((phiSource - pcEqn.flux())/rAUf);
            Uc.correctBoundaryConditions();
        }
    }

    Ucf = fvc::interpolate(Uc);
    {
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Ucf += n*(phiComp/mesh.magSf() - (n & Ucf));
    }

    // Total flux relative to the mesh (absolute phiComp + phiInc minus meshPhi)
    phi == phiComp + phiInc - meshPhi;
    #include "continuityErrs.H"

    U == Uc + incident.UInc();
    p == pc + pStar;
    p_rgh == p - rho*gh;
}
