// Log only: level of the incident volume fraction against the wave-model
// surface elevation (momentum_swense_Rev1, G3). Per column of cells over one
// (x, y):
//   level = z_bottom + sum(alphaI V)/A,   Delta = level - eta
// with A the column footprint and eta the Grid2Grid elevation of the column.
// Reported: the largest |Delta| over the columns outside the relaxation
// zones (w = 0) and Delta at the gauges of the surfaceElevation function
// object of controlDict (the column nearest to each gauge), when present.
// Serial log: on several processors each column must lie on one processor
// and the gauge values are those found by the master.
{
    static bool levInit = true;
    static labelList levColOfCell;
    static scalarField levZBot, levArea, levWMax;
    static labelList levRow;
    static wordList levGaugeName;
    static labelList levGaugeCol;

    if (levInit)
    {
        levInit = false;

        const vectorField& Cc = mesh.C().primitiveField();
        const pointField& mpts = mesh.points();
        const scalarField& Vc = mesh.V();

        HashTable<label, string> colOf;
        DynamicList<scalar> zBot, area, wMax, cx, cy;

        levColOfCell.setSize(mesh.nCells());

        forAll(Cc, c)
        {
            const labelList cp(mesh.cellPoints()[c]);
            scalar zb = GREAT, zt = -GREAT;
            forAll(cp, i)
            {
                zb = min(zb, mpts[cp[i]].z());
                zt = max(zt, mpts[cp[i]].z());
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
                area.append(Vc[c]/(zt - zb));
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
        levArea = scalarField(area);
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
    }

    scalarField colVol(levZBot.size(), 0.0);
    {
        const scalarField& Vc = mesh.V();
        forAll(levColOfCell, c)
        {
            colVol[levColOfCell[c]] += alphaI[c]*Vc[c];
        }
    }

    scalarField levDelta(levZBot.size(), 0.0);
    scalar levMaxFree = 0.0;

    forAll(levDelta, col)
    {
        levDelta[col] =
            levZBot[col] + colVol[col]/levArea[col]
          - incident.etaRow(levRow[col], runTime.value());

        if (levWMax[col] < SMALL)
        {
            levMaxFree = max(levMaxFree, mag(levDelta[col]));
        }
    }

    reduce(levMaxFree, maxOp<scalar>());

    Info<< "SWENSE level: max|level - eta| (w = 0 columns) = " << levMaxFree;
    forAll(levGaugeCol, i)
    {
        Info<< "  " << levGaugeName[i] << " = " << levDelta[levGaugeCol[i]];
    }
    Info<< endl;
}
