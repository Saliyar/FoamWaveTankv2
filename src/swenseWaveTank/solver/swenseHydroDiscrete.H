// Discrete incident hydrostatics (incidentHydrostatic discrete, see
// swenseSetup.H). Included from swenseUcEqn.H with the current rhoI:
//   pDynD   = p_side - rho_side g (eta - z)
//   pHrgh   = pHydI_rgh, marched down each column from the top cell
//             (analytic rho_side g (eta - z) - rhoI gh there) with
//             pHrgh(lower) = pHrgh(upper) - gh_f (rhoI(lower) - rhoI(upper))
//   pIncD   = pHrgh + rhoI gh + pDynD
// Boundary faces: the same relation from the owner cell,
//   pHrgh_b = pHrgh_c - gh_b (rhoI_b - rhoI_c)    (gh_b = ghf on the patch)
// so the face residual G_f of a boundary face has no hydrostatic part either.
// Moving mesh: the columns are topological (built once at set-up); each
// step the march uses the current gh_f of each face and the current eta and
// side pressure at the moved points. A column whose top face is a wall (the
// cells under a body) starts from the analytic value in its top cell.
{
    incident.sidePressure(pHydA, pDynD);

    scalarField& P = pHrgh.primitiveFieldRef();
    const scalarField& rI = rhoI.primitiveField();
    const scalarField& ghC = gh.primitiveField();

    forAll(colOrder, k)
    {
        const label c = colOrder[k];
        const label a = colAbove[c];

        if (a < 0)
        {
            P[c] = pHydA[c] - rI[c]*ghC[c];
        }
        else
        {
            P[c] = P[a] - ghf[colFace[c]]*(rI[c] - rI[a]);
        }
    }

    pHrgh.correctBoundaryConditions();

    forAll(pHrgh.boundaryField(), patchi)
    {
        fvPatchScalarField& pb = pHrgh.boundaryFieldRef()[patchi];

        if (pb.coupled() || pb.empty())
        {
            continue;
        }

        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const scalarField& rb = rhoI.boundaryField()[patchi];
        const scalarField& gb = ghf.boundaryField()[patchi];

        forAll(pb, i)
        {
            pb[i] = P[fc[i]] - gb[i]*(rb[i] - rI[fc[i]]);
        }
    }

    pIncD == pHrgh + rhoI*gh + pDynD;
}
