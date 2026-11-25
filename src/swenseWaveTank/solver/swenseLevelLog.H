// Log only (no effect on the solution). Three lines per time step.
//
// SWENSE level: level of alphaI and of alpha1 against the wave-model surface
// elevation (momentum_swense_Rev1, G3). Per column of band cells over one
// (x, y), band = cell centres with |z| < levelBandHalfHeight (controlDict,
// default 0.29 m):
//   level = z_bandBottom + sum_band(alpha V)/A,   Delta = level - eta
// The cells below the band are taken as full of water, so a column that is
// cut by a body below the band is handled. Reported: the largest |Delta| of
// alphaI over the columns outside the relaxation zones (w = 0), and Delta of
// alphaI (I) and alpha1 (1) at the gauges of the surfaceElevation function
// object of controlDict (the column nearest to each gauge), when present.
//
// SWENSE band: |Uc| in the same band, with the relaxation factor (1 - w) of
// the end of the step, split into water (alpha1 >= 0.5) and air cells:
// maximum and volume mean; and the largest |pc| in the band.
//
// SWENSE body: volume flux through the patch named body, when present:
// sum of phiInc, of phiComp, of their sum, and the sum of |phiInc + phiComp|.
//
// Serial log: on several processors each column must lie on one processor
// and the gauge values are those found by the master.
{
    static bool levInit = true;
    static labelList levColOfCell;      // -1 outside the band
    static scalarField levZBot, levWMax;
    static labelList levRow;
    static wordList levGaugeName;
    static labelList levGaugeCol;
    static scalar levHalf = 0.29;

    if (levInit)
    {
        levInit = false;

        levHalf = runTime.controlDict().lookupOrDefault<scalar>
        (
            "levelBandHalfHeight",
            0.29
        );

        const vectorField& Cc = mesh.C().primitiveField();
        const pointField& mpts = mesh.points();

        HashTable<label, string> colOf;
        DynamicList<scalar> zBot, wMax, cx, cy;

        levColOfCell.setSize(mesh.nCells(), -1);

        forAll(Cc, c)
        {
            if (mag(Cc[c].z()) >= levHalf)
            {
                continue;
            }

            const labelList cp(mesh.cellPoints()[c]);
            scalar zb = GREAT;
            forAll(cp, i)
            {
                zb = min(zb, mpts[cp[i]].z());
            }

            const string key
            (
                Foam::name(label(std::llround(Cc[c].x()*1e6))) + "_"
              + Foam::name(label(std::llround(Cc[c].y()*1e6)))
            );

            if (!colOf.found(key))
            {
                colOf.insert(key, zBot.size());
                zBot.append(zb);
                wMax.append(0.0);
                cx.append(Cc[c].x());
                cy.append(Cc[c].y());
            }

            const label col = colOf[key];
            levColOfCell[c] = col;
            zBot[col] = min(zBot[col], zb);
            wMax[col] = max(wMax[col], w[c]);
        }

        levZBot = scalarField(zBot);
        levWMax = scalarField(wMax);

        // Elevation row of each column: the row with the same (x, y)
        const pointField& rows = incident.etaPts();
        levRow.setSize(levZBot.size(), -1);

        forAll(levRow, col)
        {
            scalar d2Min = GREAT;
            forAll(rows, r)
            {
                const scalar d2 =
                    sqr(rows[r].x() - cx[col]) + sqr(rows[r].y() - cy[col]);

                if (d2 < d2Min)
                {
                    d2Min = d2;
                    levRow[col] = r;
                }
            }
        }

        // Gauges of the surfaceElevation function object (controlDict)
        if
        (
            runTime.controlDict().found("functions")
         && runTime.controlDict().subDict("functions").found("surfaceElevation")
         && runTime.controlDict().subDict("functions")
                .subDict("surfaceElevation").found("sets")
        )
        {
            const dictionary& sets =
                runTime.controlDict().subDict("functions")
                    .subDict("surfaceElevation").subDict("sets");

            const wordList names(sets.toc());

            forAll(names, i)
            {
                const point start(sets.subDict(names[i]).lookup("start"));

                label best = -1;
                scalar d2Min = GREAT;
                forAll(levZBot, col)
                {
                    const scalar d2 = sqr(cx[col] - start.x());
                    if (d2 < d2Min)
                    {
                        d2Min = d2;
                        best = col;
                    }
                }

                if (best >= 0)
                {
                    levGaugeName.append(names[i]);
                    levGaugeCol.append(best);
                }
            }
        }

        Info<< "SWENSE level: " << levZBot.size() << " band columns, band"
            << " |z| < " << levHalf << " m" << endl;
    }

    const scalarField& Vc = mesh.V();
    const vectorField& Cc = mesh.C().primitiveField();
    const pointField& mpts = mesh.points();

    // Column heights of alphaI and alpha1 in the band: sum(alpha dz), with
    // dz = cell height and A = V/dz the footprint of the cell
    scalarField colI(levZBot.size(), 0.0);
    scalarField col1(levZBot.size(), 0.0);

    // Band statistics of |Uc| (1 - w)
    scalar ucMaxW = 0.0, ucMaxA = 0.0, pcMaxBand = 0.0;
    scalar ucSumW = 0.0, ucSumA = 0.0, volW = 0.0, volA = 0.0;

    forAll(levColOfCell, c)
    {
        const label col = levColOfCell[c];

        if (col < 0)
        {
            continue;
        }

        const labelList cp(mesh.cellPoints()[c]);
        scalar zb = GREAT, zt = -GREAT;
        forAll(cp, i)
        {
            zb = min(zb, mpts[cp[i]].z());
            zt = max(zt, mpts[cp[i]].z());
        }

        colI[col] += alphaI[c]*(zt - zb);
        col1[col] += alpha1[c]*(zt - zb);

        const scalar uc = mag(Uc[c])*(1.0 - w[c]);

        if (alpha1[c] >= 0.5)
        {
            ucMaxW = max(ucMaxW, uc);
            ucSumW += uc*Vc[c];
            volW += Vc[c];
        }
        else
        {
            ucMaxA = max(ucMaxA, uc);
            ucSumA += uc*Vc[c];
            volA += Vc[c];
        }

        pcMaxBand = max(pcMaxBand, mag(pc[c]));
    }

    scalar levMaxFree = 0.0;
    scalarField dI(levZBot.size(), 0.0), d1(levZBot.size(), 0.0);

    forAll(dI, col)
    {
        const scalar eta = incident.etaRow(levRow[col], runTime.value());

        dI[col] = levZBot[col] + colI[col] - eta;
        d1[col] = levZBot[col] + col1[col] - eta;

        if (levWMax[col] < SMALL)
        {
            levMaxFree = max(levMaxFree, mag(dI[col]));
        }
    }

    reduce(levMaxFree, maxOp<scalar>());
    reduce(ucMaxW, maxOp<scalar>());
    reduce(ucMaxA, maxOp<scalar>());
    reduce(pcMaxBand, maxOp<scalar>());
    reduce(ucSumW, sumOp<scalar>());
    reduce(ucSumA, sumOp<scalar>());
    reduce(volW, sumOp<scalar>());
    reduce(volA, sumOp<scalar>());

    Info<< "SWENSE level: max|levelI - eta| (w = 0 columns) = " << levMaxFree;
    forAll(levGaugeCol, i)
    {
        Info<< "  " << levGaugeName[i] << " I = " << dI[levGaugeCol[i]]
            << " 1 = " << d1[levGaugeCol[i]]
            << " eta = "
            << incident.etaRow(levRow[levGaugeCol[i]], runTime.value());
    }
    Info<< endl;

    Info<< "SWENSE band: |z| < " << levHalf
        << "  max|Uc| water = " << ucMaxW
        << "  air = " << ucMaxA
        << "  mean|Uc| water = " << ucSumW/max(volW, VSMALL)
        << "  air = " << ucSumA/max(volA, VSMALL)
        << "  all = " << (ucSumW + ucSumA)/max(volW + volA, VSMALL)
        << "  max|pc| = " << pcMaxBand << endl;

    const label bodyPatch = mesh.boundaryMesh().findPatchID("body");

    if (bodyPatch >= 0)
    {
        const scalarField& pI = phiInc.boundaryField()[bodyPatch];
        const scalarField& pC = phiComp.boundaryField()[bodyPatch];

        Info<< "SWENSE body: sum(phiInc) = " << gSum(pI)
            << "  sum(phiComp) = " << gSum(pC)
            << "  sum(phiInc + phiComp) = " << gSum(pI + pC)
            << "  sum|phiInc + phiComp| = " << gSum(mag(pI + pC))
            << "  max|phiInc + phiComp|/|Sf| = "
            << gMax
               (
                   mag(pI + pC)
                  /mesh.magSf().boundaryField()[bodyPatch]
               )
            << endl;
    }
}
