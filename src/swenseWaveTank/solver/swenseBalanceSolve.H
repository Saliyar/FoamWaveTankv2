// Incident pressure balance, V2 (momentum_swense Part V): the pressure
// correction deltaHeld that takes the divergence out of Phi[r_I]
// (swenseBalanceFlux.H):
//
//   laplacian(rAUf, deltaHeld) = div(Phi[r_I])
//
// with the matrix structure of the pc equation (same rAUf, same laplacian
// scheme, the pcFinal solver controls). Boundaries (types set in
// swenseSetup.H): fixedValue 0 where pc fixes its value; elsewhere the
// gradient is set so that Phi_b - rAUf snGrad(deltaHeld) |Sf| = 0.
// After the solve PhiRIStar = Phi[r_I] - rAUf snGrad(deltaHeld) |Sf| is the
// flux of r_I built with pIncHeld = pInc - deltaHeld (r_I*).

forAll(deltaHeld.boundaryField(), patchi)
{
    if (isA<fixedGradientFvPatchScalarField>(deltaHeld.boundaryField()[patchi]))
    {
        refCast<fixedGradientFvPatchScalarField>
        (
            deltaHeld.boundaryFieldRef()[patchi]
        ).gradient() =
            PhiRI.boundaryField()[patchi]
           /(
                rAUf.boundaryField()[patchi]
               *mesh.magSf().boundaryField()[patchi]
            );
    }
}

surfaceScalarField PhiRIStar("PhiRIStar", PhiRI);

for (label nonOrth = 0; nonOrth <= pimple.nNonOrthCorr(); ++nonOrth)
{
    fvScalarMatrix deltaEqn
    (
        fvm::laplacian(rAUf, deltaHeld) == fvc::div(PhiRI)
    );

    if (deltaHeld.needReference())
    {
        deltaEqn.setReference(pcRefCell, scalar(0));
    }

    deltaEqn.solve(mesh.solverDict("pcFinal"));

    if (nonOrth == pimple.nNonOrthCorr())
    {
        PhiRIStar = PhiRI - deltaEqn.flux();
    }
}
