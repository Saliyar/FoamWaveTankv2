// Flux correction after a mesh motion (waveFoam correctPhi, spec item 3).
// Included after swenseIncidentStep.H on the new mesh, only when the mesh
// changes and correctPhi is on in the PIMPLE dictionary. The complementary
// flux is taken from the face velocity, then corrected by the pressure
// projection pcCorr (OpenFOAM CorrectPhi form, with rAUf of the last UcEqn,
// 1 on the first step) so that div(phiComp + phiInc) = 0. The source is
// div(phiComp), plus div(phiInc) when includeDivUInc, as in swensePcEqn.H.
phiComp = Ucf & mesh.Sf();

// As OpenFOAM correctUphiBCs: where Uc fixes its value (the body wall, the
// inlet and outlet), the boundary flux is taken from the boundary velocity.
// On the body this gives phiComp = meshPhi - phiInc, so the total flux
// relative to the mesh is zero there.
Uc.correctBoundaryConditions();

forAll(Uc.boundaryField(), patchi)
{
    if (Uc.boundaryField()[patchi].fixesValue())
    {
        phiComp.boundaryFieldRef()[patchi] =
            Uc.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
}

{
    wordList pcCorrTypes
    (
        pc.boundaryField().size(),
        fvPatchFieldBase::zeroGradientType()
    );

    // Fixed value where pc fixes its value, as CorrectPhi does
    forAll(pc.boundaryField(), patchi)
    {
        if (pc.boundaryField()[patchi].fixesValue())
        {
            pcCorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcCorr
    (
        IOobject
        (
            "pcCorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(pc.dimensions(), Zero),
        pcCorrTypes
    );

    mesh.setFluxRequired(pcCorr.name());

    if (pcCorr.needReference())
    {
        adjustPhi(phiComp, Uc, pcCorr);
    }

    while (pimple.correctNonOrthogonal())
    {
        volScalarField divRhs(fvc::div(phiComp));

        if (includeDivUInc)
        {
            divRhs += fvc::div(phiInc);
        }

        fvScalarMatrix pcCorrEqn
        (
            fvm::laplacian(rAUfPrev, pcCorr) == divRhs
        );

        if (pcCorr.needReference())
        {
            pcCorrEqn.setReference(pcRefCell, scalar(0));
        }

        pcCorrEqn.solve
        (
            mesh.solverDict
            (
                pimple.finalNonOrthogonalIter() ? "pcFinal" : "pc"
            )
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phiComp -= pcCorrEqn.flux();
        }
    }
}

// Total flux relative to the mesh, from the corrected absolute fluxes
phi = phiComp + phiInc - meshPhi;

Info<< "SWENSE correctPhi: max|div(phiComp + phiInc)| = "
    << max(mag(fvc::div(phiComp + phiInc))).value() << endl;
