// Incident fields at the current mesh and time, run after mesh.update() in
// every outer corrector that moves the mesh (and once per time step on a
// static mesh, as before). Incident fields are blended with the advected
// alphaI, which is the value at the start of this step (document Rev1, 16.3.5).
// The HOS target alphaH is the relaxation target of both volume fractions
// (16.3.3).
incident.update(runTime.value(), alphaI);

alphaH.primitiveFieldRef() =
    incident.blendAlpha(runTime.value()).primitiveField();

// Incident flux of alphaI (absolute) and the fluxes relative to the mesh
phiInc = incident.UfInc() & mesh.Sf();
phi = phiComp + phiInc - meshPhi;
phiIRel = phiInc - meshPhi;

// Projection of phiInc: laplacian(psi) = div(phiInc), psi fixed to zero on
// the atmosphere patch and zero-gradient elsewhere; phiInc -= snGrad(psi) |Sf|.
// Switch: projectPhiI in controlDict (default on).
if (runTime.controlDict().lookupOrDefault<bool>("projectPhiI", true))
{
    wordList psiTypes
    (
        mesh.boundaryMesh().size(),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(mesh.boundaryMesh(), pI)
    {
        if (mesh.boundaryMesh()[pI].name() == "atmosphere")
        {
            psiTypes[pI] = fixedValueFvPatchScalarField::typeName;
        }
    }
    // Target of div(phiInc). Momentum blend (waveProperties momentumBlend true):
    // div(rhoI uI) = rhoI qI, i.e. div(uI) = qI - uI.grad(rhoI)/rhoI with
    // qI = rhoW rhoA/rhoI^2 dU.grad(alphaI). Otherwise zero, as before.
    volScalarField divTargetI
    (
        IOobject("divTargetI", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("divTargetI", dimless/dimTime, 0.0)
    );

    if (incident.momentumBlend())
    {
        const volScalarField qI
        (
            "qI",
            rho1.value()*rho2.value()/sqr(rhoI)*(incident.dUc() & fvc::grad(alphaI))
        );
        divTargetI = qI - (incident.UInc() & fvc::grad(rhoI))/rhoI;
    }

    volScalarField psiP
    (
        IOobject("psiProj", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("psiProj", dimArea/dimTime, 0.0),
        psiTypes
    );
    // Solver settings of psiProj, also inside the outer loop where the final
    // iteration would select psiProjFinal: the projection is not a final solve
    solve
    (
        fvm::laplacian(psiP) == fvc::div(phiInc) - divTargetI,
        mesh.solverDict("psiProj")
    );
    phiInc -= fvc::snGrad(psiP)*mesh.magSf();
    phi = phiComp + phiInc - meshPhi;
    phiIRel = phiInc - meshPhi;
}
