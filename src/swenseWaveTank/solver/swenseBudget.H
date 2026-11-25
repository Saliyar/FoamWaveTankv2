// Logging only (switch budgetTerms in controlDict). No equation, field or
// tolerance is changed. Terms of the Uc equation on the final fields of the
// step, as the solver writes them (swenseUcEqn.H), per unit mass:
//
//   T1  ddt(rho, Uc)
//   T2  div(rhoPhi, Uc)
//   T3  ddt(rho, UInc) - ddt(rhoI, UInc)            (= ddt(Drho UInc) form)
//   T4  div(rhoPhi, UInc) - div(rhoPhiI, UInc)      (incident convection difference)
//   T5  reconstruct((ghf snGrad(Drho) + snGrad(pc)) magSf)   (pressure and gravity)
//   T6  viscous: divDevRhoReff(rho, Uc) acting on Uc, plus the incident
//       viscous difference with a minus sign (as in incTerm)
//   R   = T1 + ... + T6, the momentum residual of the solver's own discrete
//         equation (LHS + incTerm - RHS), on the final fields
//
// Set S: class C cells (alpha1 < 0.01 and alphaI < 0.01) with 0 <= z - eta <= 0.2.
// Eta of a cell is the Grid2Grid elevation of the surface row nearest to it.
{
    static bool budgetInit = true;
    static labelList budgetRow;

    if (budgetInit)
    {
        budgetInit = false;

        const pointField& rows = incident.etaPts();
        const vectorField& Cc = mesh.C().primitiveField();

        budgetRow.setSize(mesh.nCells(), 0);

        forAll(Cc, c)
        {
            scalar d2Min = GREAT;

            forAll(rows, r)
            {
                const scalar d2 =
                    sqr(rows[r].x() - Cc[c].x()) + sqr(rows[r].y() - Cc[c].y());

                if (d2 < d2Min)
                {
                    d2Min = d2;
                    budgetRow[c] = r;
                }
            }
        }
    }

    const volVectorField& UInc = incident.UInc();

    const scalar tB = runTime.value();
    const vectorField& CcB = mesh.C().primitiveField();
    const scalarField& a1B = alpha1.primitiveField();
    const scalarField& aIB = alphaI.primitiveField();
    const scalarField& rhoB = rho.primitiveField();

    // Stop rule first: max|Uc| above 3 or any NaN ends the run
    {
        const scalar maxUcB = max(mag(Uc)).value();

        if (maxUcB > 3.0 || maxUcB != maxUcB)
        {
            Info<< "BUDGET STOP t = " << tB << " max|Uc| = " << maxUcB
                << " (limit 3 or NaN)" << endl;
            runTime.writeAndEnd();
        }
    }

    // Fluxes and incident fields at the write times (full fields for item 4)
    if (runTime.writeTime())
    {
        rhoPhi.write();
        rhoPhiI.write();
        phiInc.write();
        phi.write();
        Drho.write();
        UInc.write();
    }

    // Fluxes and fields of the momentum terms
    const volVectorField T1("budT1", fvc::ddt(rho, Uc));
    const volVectorField T2("budT2", fvc::div(rhoPhi, Uc));

    const fvVectorMatrix McB(turbulence->divDevRhoReff(rho, Uc));
    const volVectorField T6c("budT6c", McB & Uc);

    const word convNameB("div(rhoPhi,Uc)");

    const volVectorField T3
    (
        "budT3",
        fvc::ddt(rho, UInc) - fvc::ddt(rhoI, UInc)
    );

    const volVectorField T4
    (
        "budT4",
        fvc::div(rhoPhi, UInc, convNameB) - fvc::div(rhoPhiI, UInc, convNameB)
    );

    const volVectorField T5
    (
        "budT5",
        fvc::reconstruct
        (
            (ghf*fvc::snGrad(Drho) + fvc::snGrad(pc))*mesh.magSf()
        )
    );

    // Incident viscous difference, the same expressions as swenseUcEqn.H
    auto muOfB = [&](const volScalarField& s) -> tmp<volScalarField>
    {
        const volScalarField lim(min(max(s, scalar(0)), scalar(1)));

        return
            lim*rho1*mixture.nuModel1().nu()
          + (scalar(1) - lim)*rho2*mixture.nuModel2().nu();
    };

    const volTensorField gradUIB("budGradUI", fvc::grad(UInc));

    auto viscUIB = [&](const volScalarField& mu) -> tmp<volVectorField>
    {
        return
            fvc::laplacian(mu, UInc, "laplacian(mu,UInc)")
          + fvc::div(mu*dev2(T(gradUIB)), "div(mu*dev2(T(grad(UInc))))");
    };

    const volVectorField T6i
    (
        "budT6i",
        -(viscUIB(muOfB(alpha1)) - viscUIB(muOfB(alphaI)))
    );

    const volVectorField R("budR", T1 + T2 + T3 + T4 + T5 + T6c + T6i);

    const volScalarField gradAfter("budGradAfter", mag(fvc::grad(UInc)));

    // Set S and the per-step statistics
    const label nC = mesh.nCells();
    const scalarField& gA = gradAfter.primitiveField();

    const char* names[8] = {"T1", "T2", "T3", "T4", "T5", "T6c", "T6i", "R"};

    std::vector<const volVectorField*> fieldsB
    {
        &T1, &T2, &T3, &T4, &T5, &T6c, &T6i, &R
    };

    scalar maxUcS = 0.0, maxGradAfter = 0.0, maxGradRaw = 0.0;
    label nS = 0, cellMax = -1;

    scalarList mxT(8, 0.0), ssT(8, 0.0);

    for (label c = 0; c < nC; ++c)
    {
        const scalar zRel = CcB[c].z() - incident.etaRow(budgetRow[c], tB);

        if (a1B[c] < 0.01 && aIB[c] < 0.01 && zRel >= 0.0 && zRel <= 0.2)
        {
            ++nS;

            const scalar ucm = mag(Uc.primitiveField()[c]);

            if (ucm > maxUcS)
            {
                maxUcS = ucm;
                cellMax = c;
            }

            maxGradAfter = max(maxGradAfter, gA[c]);
            maxGradRaw = max(maxGradRaw, budgetRawGradMag[c]);

            for (label k = 0; k < 8; ++k)
            {
                const scalar v =
                    mag((*fieldsB[k]).primitiveField()[c])/rhoB[c];

                mxT[k] = max(mxT[k], v);
                ssT[k] += v*v;
            }
        }
    }

    reduce(nS, sumOp<label>());
    reduce(maxUcS, maxOp<scalar>());
    reduce(maxGradAfter, maxOp<scalar>());
    reduce(maxGradRaw, maxOp<scalar>());

    for (label k = 0; k < 8; ++k)
    {
        reduce(mxT[k], maxOp<scalar>());
        reduce(ssT[k], sumOp<scalar>());
    }

    if (nS > 0)
    {
        Info<< "BUDGET_STATS t " << tB << " nS " << nS
            << " maxUc " << maxUcS
            << " maxGradAfter " << maxGradAfter
            << " maxGradRaw " << maxGradRaw
            << " ratio " << (maxGradRaw > 0 ? maxGradAfter/maxGradRaw : 0.0);

        for (label k = 0; k < 8; ++k)
        {
            Info<< " " << names[k] << "_max " << mxT[k]
                << " " << names[k] << "_rms "
                << Foam::sqrt(ssT[k]/nS);
        }

        Info<< endl;
    }

    if (cellMax >= 0)
    {
        const vector Uk = Uc.primitiveField()[cellMax];
        const scalar rk = rhoB[cellMax];

        Info<< "BUDGET_CELL t " << tB << " cell " << cellMax
            << " x " << CcB[cellMax].x() << " z " << CcB[cellMax].z()
            << " Uc " << Uk
            << " T1 " << T1.primitiveField()[cellMax]/rk
            << " T2 " << T2.primitiveField()[cellMax]/rk
            << " T3 " << T3.primitiveField()[cellMax]/rk
            << " T4 " << T4.primitiveField()[cellMax]/rk
            << " T5 " << T5.primitiveField()[cellMax]/rk
            << " T6c " << T6c.primitiveField()[cellMax]/rk
            << " T6i " << T6i.primitiveField()[cellMax]/rk
            << " R " << R.primitiveField()[cellMax]/rk
            << " dK_T1 " << (Uk & (T1.primitiveField()[cellMax]/rk))
            << " dK_T2 " << (Uk & (T2.primitiveField()[cellMax]/rk))
            << " dK_T3 " << (Uk & (T3.primitiveField()[cellMax]/rk))
            << " dK_T4 " << (Uk & (T4.primitiveField()[cellMax]/rk))
            << " dK_T5 " << (Uk & (T5.primitiveField()[cellMax]/rk))
            << " dK_T6c " << (Uk & (T6c.primitiveField()[cellMax]/rk))
            << " dK_T6i " << (Uk & (T6i.primitiveField()[cellMax]/rk))
            << " dK_R " << (Uk & (R.primitiveField()[cellMax]/rk))
            << endl;
    }
}
