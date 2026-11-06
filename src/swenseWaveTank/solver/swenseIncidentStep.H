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
    // Target of div(phiInc) (theory Sec. 10): kinematic identity div(uI) = q,
    // q = qI (momentum blend) or qv (velocity blend).
    volScalarField divTargetI
    (
        IOobject("divTargetI", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("divTargetI", dimless/dimTime, 0.0)
    );

    if (incident.momentumBlend())
    {
        // Kinematic target qI (eq. kin): div(uI) = qI, no steady mass term
        divTargetI = rho1*rho2/sqr(rhoI)*(incident.dUc() & fvc::grad(alphaI));
    }
    else
    {
        // Velocity blend: div(uI) = qv = dU.grad(alphaI) (eq. qv)
        divTargetI = incident.dUc() & fvc::grad(alphaI);
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
