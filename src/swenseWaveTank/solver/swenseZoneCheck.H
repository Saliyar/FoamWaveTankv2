// Relaxation zones keep their t = 0 weights (swenseRelaxation.C), so the mesh
// must not move inside them (spec item 8): a cell with a relaxation weight
// above zero may not have a point displaced by more than 1e-12 m since t = 0.
{
    const pointField& ptsNow = mesh.points();

    scalar dZone = 0.0;
    forAll(w, c)
    {
        if (w[c] > 0)
        {
            const labelList cp = mesh.cellPoints(c);
            forAll(cp, i)
            {
                dZone = max(dZone, mag(ptsNow[cp[i]] - points0[cp[i]]));
            }
        }
    }

    if (dZone > 1e-12)
    {
        FatalErrorInFunction
            << "The mesh moves inside a relaxation zone (max point displacement "
            << dZone << " m since t = 0). The zone weights are fixed at t = 0."
            << exit(FatalError);
    }
}
