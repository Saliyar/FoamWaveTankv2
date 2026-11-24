// alphaHExact on: log max|alphaH - blend weight| over the cells (the blend
// target of alphaHExact off) at each step
if (swenseIncident::alphaHExact(mesh.thisDb()))
{
    const scalarField dB
    (
        mag
        (
            alphaH.primitiveField()
          - incident.blendAlpha(runTime.value()).primitiveField()
        )
    );
    Info<< "SWENSE alphaHExact: max|exact - blend| = " << gMax(dB) << endl;
}
