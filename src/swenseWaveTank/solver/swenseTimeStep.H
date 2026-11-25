// Courant numbers and adaptive time step (as waveIsoFoam: CourantNo.H,
// alphaCourantNo.H and VoF/setDeltaT.H), called before runTime++.
//
// The Courant numbers use the total flux relative to the mesh,
// phi = phiComp + phiInc - meshPhi (the flux that advects alpha1), of the
// previous step, and are always logged in the format of the plain solvers:
//   Courant Number mean: ... max: ...
//   Interface Courant Number mean: ... max: ...
// controlDict (read every step):
//   adjustTimeStep  default no: the time step stays as set (unchanged runs)
//   maxCo           default 1
//   maxAlphaCo      default maxCo
//   maxDeltaT       default GREAT
// With adjustTimeStep yes the step follows setDeltaT.H of the VoF solvers:
// reduction immediate, increase damped (at most 1.2 per step).
{
    const dictionary& cd = runTime.controlDict();
    const bool adjustTimeStepOn = cd.getOrDefault<bool>("adjustTimeStep", false);
    const scalar maxCoS = cd.getOrDefault<scalar>("maxCo", 1);
    const scalar maxAlphaCoS = cd.getOrDefault<scalar>("maxAlphaCo", maxCoS);
    const scalar maxDeltaTS = cd.getOrDefault<scalar>("maxDeltaT", GREAT);

    scalar CoNumS = 0;
    scalar meanCoNumS = 0;
    scalar alphaCoNumS = 0;
    scalar meanAlphaCoNumS = 0;

    if (mesh.nInternalFaces())
    {
        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );
        const scalarField sumPhiI
        (
            mixture.nearInterface()().primitiveField()*sumPhi
        );
        const scalar dt = runTime.deltaTValue();

        CoNumS = 0.5*gMax(sumPhi/mesh.V().field())*dt;
        meanCoNumS = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*dt;
        alphaCoNumS = 0.5*gMax(sumPhiI/mesh.V().field())*dt;
        meanAlphaCoNumS = 0.5*(gSum(sumPhiI)/gSum(mesh.V().field()))*dt;
    }

    Info<< "Courant Number mean: " << meanCoNumS
        << " max: " << CoNumS << endl;
    Info<< "Interface Courant Number mean: " << meanAlphaCoNumS
        << " max: " << alphaCoNumS << endl;

    if (adjustTimeStepOn)
    {
        const scalar maxDeltaTFact =
            Foam::min
            (
                maxCoS/(CoNumS + SMALL),
                maxAlphaCoS/(alphaCoNumS + SMALL)
            );

        const scalar deltaTFact =
            Foam::min(Foam::min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        runTime.setDeltaT
        (
            Foam::min(deltaTFact*runTime.deltaTValue(), maxDeltaTS)
        );

        Info<< "deltaT = " << runTime.deltaTValue() << endl;
    }
}
