// Exact incident reset (controlDict alphaIExactReset, default off; set-up and
// checks in swenseSetup.H). Called in swenseAlphaRelax.H after the relaxation,
// in place of the level correction, before the densities and
// swenseIncidentStep.H. In every cell outside plainBodyZone:
//   aX = swenseIncident::exactAlpha(t)  (the same call as targetAlpha, the
//        relaxation target with waveProperties alphaHExact on)
//   dA = aX - alphaI;  alphaI = aX;  alpha1 = alpha1 + dA
// alpha1 + dA is written as aX + (alpha1 - alphaI), so alpha1 = aX exactly
// where alpha1 = alphaI. alpha1 is clipped to [0, 1]; the clipped volume is
// logged. Log per step: net and abs reset volume (sum dA V, sum |dA| V),
// max|dA| with its position, split into free zone (w = 0) and relaxation
// zones (w > 0), and running totals.
{
    const volScalarField alphaX(incident.exactAlpha(runTime.value()));
    const scalarField& aX = alphaX.primitiveField();

    const scalarField& V = mesh.V();
    const scalarField& wC = w.primitiveField();
    const vectorField& Cc = mesh.C().primitiveField();
    scalarField& aI = alphaI.primitiveFieldRef();
    scalarField& a1 = alpha1.primitiveFieldRef();

    scalar netV = 0, absV = 0, clipV = 0;
    scalar freeNet = 0, freeAbs = 0, relaxNet = 0, relaxAbs = 0;
    scalar dMax = -1;
    vector xMax(Zero);

    forAll(aI, c)
    {
        if (exactResetSkip.size() && exactResetSkip[c])
        {
            continue;
        }

        const scalar dA = aX[c] - aI[c];
        const scalar a1n = aX[c] + (a1[c] - aI[c]);
        const scalar a1c = min(max(a1n, scalar(0)), scalar(1));

        aI[c] = aX[c];
        a1[c] = a1c;

        const scalar dV = dA*V[c];
        clipV += mag(a1n - a1c)*V[c];
        netV += dV;
        absV += mag(dV);
        if (wC[c] > 0)
        {
            relaxNet += dV;
            relaxAbs += mag(dV);
        }
        else
        {
            freeNet += dV;
            freeAbs += mag(dV);
        }
        if (mag(dA) > dMax)
        {
            dMax = mag(dA);
            xMax = Cc[c];
        }
    }

    reduce(netV, sumOp<scalar>());
    reduce(absV, sumOp<scalar>());
    reduce(clipV, sumOp<scalar>());
    reduce(freeNet, sumOp<scalar>());
    reduce(freeAbs, sumOp<scalar>());
    reduce(relaxNet, sumOp<scalar>());
    reduce(relaxAbs, sumOp<scalar>());
    const scalar dLoc = dMax;
    reduce(dMax, maxOp<scalar>());
    scalar xM = (dLoc == dMax) ? xMax.x() : -GREAT;
    scalar zM = (dLoc == dMax) ? xMax.z() : -GREAT;
    reduce(xM, maxOp<scalar>());
    reduce(zM, maxOp<scalar>());

    exactResetTotal += absV;
    exactResetClipTotal += clipV;

    alphaI.correctBoundaryConditions();
    alpha1.correctBoundaryConditions();

    Info<< "SWENSE alphaIExactReset: net = " << netV
        << " m3, abs = " << absV
        << " m3, max|dA| = " << dMax << " at x = " << xM << " z = " << zM
        << ", clipped alpha1 = " << clipV
        << " m3; totals abs " << exactResetTotal
        << " clipped " << exactResetClipTotal << " m3" << nl
        << "SWENSE alphaIExactReset zones: free (w = 0) net = " << freeNet
        << " abs = " << freeAbs << " m3, relax (w > 0) net = " << relaxNet
        << " abs = " << relaxAbs << " m3" << endl;
}
