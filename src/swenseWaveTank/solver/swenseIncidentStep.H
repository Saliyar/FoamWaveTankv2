// Incident field of the new time level (momentum_swense_Rev1, Section 12.1,
// steps 1b to 1d), for the wave models of the last incident.evaluate() and
// the current alphaI. Included once at set-up (level 0) and once per time
// step, after the transport and relaxation of alphaI.
{
    // 1b, 1c: momentum blend in the cells, raw face velocity on the faces
    incident.blend(alphaI);

    // 1/rhoI_f with rhoI_f = a_f rho_w + (1 - a_f) rho_a, a_f the linear face
    // value of alphaI (rhoI is linear in alphaI)
    const surfaceScalarField rRhoIf("rRhoIf", 1.0/linearInterpolate(rhoI));

    // Raw incident flux (step 1c)
    const surfaceScalarField phiRaw("phiRaw", incident.UfInc() & mesh.Sf());

    // 1d: density-weighted projection, eq. (31):
    //   sum_f (1/rhoI_f) snGrad(psi)_f |Sf| = sum_f phiRaw_f,
    // psi = 0 on the atmosphere patch, zero gradient elsewhere. Target
    // divergence zero. Solver settings: psiProj in fvSolution.
    label psiIters = 0;

    for (label nonOrth = 0; nonOrth <= pimple.nNonOrthCorr(); ++nonOrth)
    {
        fvScalarMatrix psiEqn
        (
            fvm::laplacian(rRhoIf, psiProj) == fvc::div(phiRaw)
        );

        psiIters += psiEqn.solve(mesh.solverDict("psiProj")).nIterations();

        if (nonOrth == pimple.nNonOrthCorr())
        {
            // Flux correction (1/rhoI_f) snGrad(psi) |Sf|
            const surfaceScalarField dPhiI("dPhiI", psiEqn.flux());

            phiInc = phiRaw - dPhiI;

            // Cell velocity: uI - (1/rhoI) R[snGrad(psi) |Sf|]
            incident.UIncRW() -= fvc::reconstruct(dPhiI/rRhoIf)/rhoI;
        }
    }

    phiIRel = phiInc - meshPhi;

    // Log: iterations and the size of the flux correction per unit area,
    // split by the face fraction a_f (internal faces)
    {
        const surfaceScalarField afLog
        (
            linearInterpolate(min(max(alphaI, scalar(0)), scalar(1)))
        );
        const surfaceScalarField dUf(mag(phiRaw - phiInc)/mesh.magSf());

        scalar dWet = 0.0, dMid = 0.0, dDry = 0.0;

        forAll(dUf, facei)
        {
            if (afLog[facei] > 0.1)
            {
                dWet = max(dWet, dUf[facei]);
            }
            else if (afLog[facei] < 0.01)
            {
                dDry = max(dDry, dUf[facei]);
            }
            else
            {
                dMid = max(dMid, dUf[facei]);
            }
        }

        reduce(dWet, maxOp<scalar>());
        reduce(dMid, maxOp<scalar>());
        reduce(dDry, maxOp<scalar>());

        Info<< "SWENSE proj: psi iterations = " << psiIters
            << "  max|phiRaw - phiInc|/|Sf|: a_f > 0.1: " << dWet
            << "  0.01 <= a_f <= 0.1: " << dMid
            << "  a_f < 0.01: " << dDry
            << "  max|div(phiRaw)| = " << max(mag(fvc::div(phiRaw))).value()
            << "  max|div(phiInc)| = " << max(mag(fvc::div(phiInc))).value()
            << endl;
    }
}
