// Zone-edge distance log (waveProperties logZoneEdge, default off; needs
// plainBodyZone). Every step: the smallest distance from a body wall face
// centre to a zone-edge face centre (internal face with one zone cell), in
// cells: distance / (V/|Sf|) of the zone cell of that edge face (its length
// normal to the face). The zone is a cellZone, so it moves with the mesh.
// The first step below 1 cell is marked BELOW1. Logging only.
{
    const labelList& zc = mesh.cellZones()[plainBodyZoneID];
    boolList inZ(mesh.nCells(), false);
    for (const label c : zc)
    {
        inZ[c] = true;
    }

    scalar dMin = GREAT, dMinCells = GREAT;
    point pBody(Zero), pEdge(Zero);

    for (label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        const label o = mesh.owner()[f];
        const label n = mesh.neighbour()[f];

        if (inZ[o] == inZ[n])
        {
            continue;
        }

        const label cz = inZ[o] ? o : n;
        const scalar len = mesh.V()[cz]/mesh.magSf()[f];
        const point& cf = mesh.Cf()[f];

        for (const label patchi : bodyPatches)
        {
            const vectorField& Cb = mesh.Cf().boundaryField()[patchi];

            forAll(Cb, i)
            {
                const scalar d = mag(Cb[i] - cf);

                if (d/len < dMinCells)
                {
                    dMinCells = d/len;
                    dMin = d;
                    pBody = Cb[i];
                    pEdge = cf;
                }
            }
        }
    }

    reduce(dMinCells, minOp<scalar>());

    Info<< "SWENSE zoneEdge: t " << runTime.value()
        << " minDistCells " << dMinCells << " dist " << dMin
        << " body " << pBody << " edge " << pEdge << endl;

    if (dMinCells < 1 && !zoneEdgeBelow1)
    {
        zoneEdgeBelow1 = true;
        Info<< "SWENSE zoneEdge BELOW1 t " << runTime.value()
            << " minDistCells " << dMinCells << endl;
    }
}
