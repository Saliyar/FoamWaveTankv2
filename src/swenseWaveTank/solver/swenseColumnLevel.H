// Level correction of the incident surface (SWENSE, document Rev1 16.3.3).
//
// Per column (vertical stack of cells over one (x, y)) of the cells whose
// bottom face is at or above levelCorrectionZMin (controlDict; default -GREAT,
// i.e. every cell, the old behaviour):
//   Delta = z_bottom + sum(alphaI dz) - eta_H,  eta_H the Grid2Grid eta of the
//   surface row of the column at the current time,
//   alphaIL = alphaI shifted by s = |Delta|/dz_c toward the neighbour above
//   (Delta > 0) or below (Delta < 0); top and bottom cells use themselves,
//   inc = alphaIL - alphaI, added to alphaI and to alpha1.
// alpha1 is clipped to [0, 1] after the increment; alphaI is not clipped.
// Cells below levelCorrectionZMin are not corrected: the formula assumes
// alphaI = 1 there, which is checked every step.
// Optional skip box (levelCorrectionSkipMin/Max in controlDict): columns that
// overlap it, or that lack cells (box cut, refinement), are not corrected.
// Parallel: every column must lie on one processor (checked at set-up). The
// column sums are then local; the global numbers (max|Delta|, clipped volume,
// counts) are reduced, and logs are written by the master.
// The relaxation densities are computed after this, in swenseAlphaRelax.H.

static bool colInit = true;
static labelList colOfCell;     // -1 for cells below levelCorrectionZMin
static labelList nbUp;
static labelList nbDown;
static scalarField dzCell;
static scalarField zBotCol;
static labelList rowOfCol;
static scalarField colX, colY;
static label nCol = 0;
static scalar zLevelMin = -GREAT;
static labelList cellsBelow;
static scalar clipTotal = 0.0;
static stringList colKeyList;
static bool colLogOn = false;
static bool colLogDone = false;

if (colInit)
{
    colInit = false;

    colLogOn = runTime.controlDict().lookupOrDefault<bool>
    (
        "levelCorrectionColumnLog",
        false
    );

    zLevelMin = runTime.controlDict().lookupOrDefault<scalar>
    (
        "levelCorrectionZMin",
        -GREAT
    );

    const vectorField& Cc = mesh.C().primitiveField();
    const pointField& pts = mesh.points();

    // Cell bottoms and tops, cell heights
    dzCell.setSize(mesh.nCells());
    scalarField zBotCell(mesh.nCells());
    scalarField zTopCell(mesh.nCells());
    forAll(dzCell, c)
    {
        scalar zmin = GREAT;
        scalar zmax = -GREAT;
        const labelList cp = mesh.cellPoints(c);
        forAll(cp, i)
        {
            zmin = min(zmin, pts[cp[i]].z());
            zmax = max(zmax, pts[cp[i]].z());
        }
        dzCell[c] = zmax - zmin;
        zBotCell[c] = zmin;
        zTopCell[c] = zmax;
    }

    // Columns: active cells (bottom at or above zLevelMin) with the same (x, y).
    // Keys and cell counts are those of this processor.
    HashTable<label, string> colMap;
    DynamicList<point> colXY;
    DynamicList<string> colKeyL;
    DynamicList<label> colCountL;
    colOfCell.setSize(mesh.nCells(), -1);
    DynamicList<label> below;
    forAll(Cc, c)
    {
        if (zBotCell[c] < zLevelMin - 1e-9)
        {
            below.append(c);
            continue;
        }

        const string key
        (
            name(lround(Cc[c].x()*1e6)) + "_" + name(lround(Cc[c].y()*1e6))
        );
        if (!colMap.found(key))
        {
            colMap.insert(key, colMap.size());
            colXY.append(Cc[c]);
            colKeyL.append(key);
            colCountL.append(0);
        }
        colOfCell[c] = colMap[key];
        colCountL[colOfCell[c]]++;
    }
    nCol = colMap.size();
    cellsBelow = below;

    // Every column must lie on one processor. Each rank sends its keys and
    // cell counts to all ranks, and every rank makes the same check, so all
    // ranks stop together. In serial the list holds only this rank: no split.
    {
        List<stringList> keysAll(Pstream::nProcs());
        List<labelList> cellsAll(Pstream::nProcs());
        keysAll[Pstream::myProcNo()] = stringList(colKeyL);
        cellsAll[Pstream::myProcNo()] = labelList(colCountL);
        Pstream::allGatherList(keysAll);
        Pstream::allGatherList(cellsAll);

        HashTable<label, string> procOf;
        DynamicList<string> splitKeys;
        forAll(keysAll, proc)
        {
            forAll(keysAll[proc], i)
            {
                const string& key = keysAll[proc][i];
                if (!procOf.found(key))
                {
                    procOf.insert(key, proc);
                }
                else if (procOf[key] != proc && !splitKeys.found(key))
                {
                    splitKeys.append(key);
                }
            }
        }

        if (splitKeys.size())
        {
            const label nShow = min(label(5), label(splitKeys.size()));
            string shown;
            for (label i = 0; i < nShow; i++)
            {
                shown += splitKeys[i] + " ";
            }

            FatalErrorInFunction
                << "columns split across processors, use method manual with "
                << "whole columns per rank, or simple (N 1 1) with N dividing "
                << "the number of columns (here " << procOf.size() << ")" << nl
                << "split columns: " << splitKeys.size()
                << ", first " << nShow << " keys (x_y, micrometres): "
                << shown << exit(FatalError);
        }
    }

    // Skip box (optional: controlDict levelCorrectionSkipMin and
    // levelCorrectionSkipMax). A column is skipped when one of its cells
    // overlaps the box, or when it has fewer cells than the full column (the
    // box cut and the refinement around it remove cells from some columns).
    // Skipped cells keep colOfCell -1: no increment and no clipping there.
    // Without the skip box every column is corrected, as before.
    if (runTime.controlDict().found("levelCorrectionSkipMin"))
    {
        const vector bMin
        (
            runTime.controlDict().lookupOrDefault<vector>
            (
                "levelCorrectionSkipMin",
                vector::zero
            )
        );
        const vector bMax
        (
            runTime.controlDict().lookupOrDefault<vector>
            (
                "levelCorrectionSkipMax",
                vector::zero
            )
        );

        labelList nCellsRaw(nCol, 0);
        labelList skipBox(nCol, 0);
        forAll(colOfCell, c)
        {
            const label col = colOfCell[c];
            if (col < 0) continue;

            nCellsRaw[col]++;

            scalar xa = GREAT, xb = -GREAT, ya = GREAT, yb = -GREAT;
            const labelList cp = mesh.cellPoints(c);
            forAll(cp, i)
            {
                xa = min(xa, pts[cp[i]].x());
                xb = max(xb, pts[cp[i]].x());
                ya = min(ya, pts[cp[i]].y());
                yb = max(yb, pts[cp[i]].y());
            }
            if
            (
                xb > bMin.x() && xa < bMax.x()
             && yb > bMin.y() && ya < bMax.y()
             && zTopCell[c] > bMin.z() && zBotCell[c] < bMax.z()
            )
            {
                skipBox[col] = 1;
            }
        }

        // Full column height: the largest count over all processors
        label nFull = 0;
        forAll(nCellsRaw, col)
        {
            nFull = max(nFull, nCellsRaw[col]);
        }
        reduce(nFull, maxOp<label>());

        labelList newCol(nCol, -1);
        DynamicList<point> colXYkeep;
        DynamicList<string> colKeyKeep;
        label nKeep = 0;
        label nBoxSkip = 0;
        label nIncSkip = 0;
        forAll(newCol, col)
        {
            if (skipBox[col])
            {
                nBoxSkip++;
                continue;
            }
            if (nCellsRaw[col] != nFull)
            {
                nIncSkip++;
                continue;
            }
            newCol[col] = nKeep++;
            colXYkeep.append(colXY[col]);
            colKeyKeep.append(colKeyL[col]);
        }
        forAll(colOfCell, c)
        {
            if (colOfCell[c] >= 0)
            {
                colOfCell[c] = newCol[colOfCell[c]];
            }
        }
        colXY = colXYkeep;
        colKeyL = colKeyKeep;
        nCol = nKeep;

        reduce(nBoxSkip, sumOp<label>());
        reduce(nIncSkip, sumOp<label>());

        Info<< "level correction: skip box min " << bMin << " max " << bMax
            << ", skipped columns " << nBoxSkip << " (overlap the box) and "
            << nIncSkip << " (incomplete, " << nFull << " cells expected), "
            << "corrected columns " << returnReduce(nCol, sumOp<label>())
            << endl;
    }

    colKeyList = stringList(colKeyL);

    colX.setSize(nCol);
    colY.setSize(nCol);
    forAll(colXY, col)
    {
        colX[col] = colXY[col].x();
        colY[col] = colXY[col].y();
    }

    // Column bottoms, cell counts and tops; the same count and a top at the
    // top of the mesh are required for every column
    zBotCol.setSize(nCol, GREAT);
    scalarField zTopCol(nCol, -GREAT);
    labelList nCellsCol(nCol, 0);
    forAll(colOfCell, c)
    {
        const label col = colOfCell[c];
        if (col < 0) continue;
        zBotCol[col] = min(zBotCol[col], zBotCell[c]);
        zTopCol[col] = max(zTopCol[col], zTopCell[c]);
        nCellsCol[col]++;
    }

    // Reference count: the largest over all columns of all processors
    label nCellsRef = 0;
    forAll(nCellsCol, col)
    {
        nCellsRef = max(nCellsRef, nCellsCol[col]);
    }
    reduce(nCellsRef, maxOp<label>());

    const scalar zTopMesh = mesh.bounds().max().z();
    forAll(nCellsCol, col)
    {
        if (nCellsCol[col] != nCellsRef)
        {
            FatalErrorInFunction
                << "column " << col << " at x = " << colX[col]
                << " has " << nCellsCol[col] << " cells above zLevelMin, "
                << "the largest column has " << nCellsRef
                << exit(FatalError);
        }
        if (mag(zTopCol[col] - zTopMesh) > 1e-9)
        {
            FatalErrorInFunction
                << "column " << col << " at x = " << colX[col]
                << " does not reach the top of the mesh (top " << zTopCol[col]
                << ", mesh top " << zTopMesh << ")" << exit(FatalError);
        }
    }

    // Cells of each column sorted by height: neighbours above and below
    List<DynamicList<label>> colCells(nCol);
    forAll(colOfCell, c)
    {
        if (colOfCell[c] >= 0)
        {
            colCells[colOfCell[c]].append(c);
        }
    }
    nbUp.setSize(mesh.nCells(), -1);
    nbDown.setSize(mesh.nCells(), -1);
    forAll(colCells, col)
    {
        const labelList cc(colCells[col]);
        scalarField zc(cc.size());
        forAll(cc, i)
        {
            zc[i] = Cc[cc[i]].z();
        }
        labelList order;
        sortedOrder(zc, order);
        forAll(order, i)
        {
            const label c = cc[order[i]];
            nbDown[c] = (i > 0) ? cc[order[i - 1]] : c;
            nbUp[c] = (i + 1 < order.size()) ? cc[order[i + 1]] : c;
        }
    }

    // Surface row of each column: the row with the same (x, y)
    const pointField& etaPts = incident.etaPts();
    rowOfCol.setSize(nCol);
    forAll(rowOfCol, col)
    {
        scalar d2Min = GREAT;
        label rMin = -1;
        forAll(etaPts, r)
        {
            const scalar d2 =
                sqr(etaPts[r].x() - colX[col]) + sqr(etaPts[r].y() - colY[col]);
            if (d2 < d2Min)
            {
                d2Min = d2;
                rMin = r;
            }
        }
        if (Foam::sqrt(d2Min) > 1e-6)
        {
            FatalErrorInFunction
                << "column " << col << " has no surface row (distance "
                << Foam::sqrt(d2Min) << ")" << exit(FatalError);
        }
        rowOfCol[col] = rMin;
    }

    Info<< "level correction: levelCorrectionZMin = " << zLevelMin
        << ", columns = " << returnReduce(nCol, sumOp<label>())
        << ", cells per column = " << nCellsRef
        << ", cells below = " << returnReduce(cellsBelow.size(), sumOp<label>())
        << endl;
}

// Cells below zLevelMin must stay fully water (the column formula assumes it).
// The test is global, so a processor without such cells still takes part.
if (returnReduce(cellsBelow.size(), sumOp<label>()))
{
    const scalarField& aB = alphaI.primitiveField();
    scalar maxBelow = -GREAT;
    forAll(cellsBelow, i)
    {
        maxBelow = max(maxBelow, 1.0 - aB[cellsBelow[i]]);
    }
    reduce(maxBelow, maxOp<scalar>());
    Info<< "level correction: max(1 - alphaI) below zLevelMin = " << maxBelow
        << endl;
    if (maxBelow > 1e-6)
    {
        FatalErrorInFunction
            << "alphaI is not 1 below levelCorrectionZMin: max(1 - alphaI) = "
            << maxBelow << exit(FatalError);
    }
}

const scalarField aI0(alphaI.primitiveField());
const scalar tNow = mesh.time().value();

// Delta per column, before the correction
scalarField colSum(nCol, 0.0);
forAll(colOfCell, c)
{
    if (colOfCell[c] >= 0)
    {
        colSum[colOfCell[c]] += aI0[c]*dzCell[c];
    }
}
scalarField Delta(nCol, 0.0);
forAll(Delta, col)
{
    Delta[col] = zBotCol[col] + colSum[col] - incident.etaRow(rowOfCol[col], tNow);
}
scalar dMax = 0.0;
scalar xDMax = 0.0;
forAll(Delta, col)
{
    if (mag(Delta[col]) > dMax)
    {
        dMax = mag(Delta[col]);
        xDMax = colX[col];
    }
}

// Moving mesh (spec item 7): a column is skipped, as the skip box does, when
// a cell of it has moved since t = 0 by more than 1e-3 of its own dz. On a
// static mesh no column is flagged, so the correction is unchanged.
labelList colMoved(nCol, 0);
if (mesh.moving())
{
    const pointField& ptsNow = mesh.points();

    forAll(colOfCell, c)
    {
        const label col = colOfCell[c];
        if (col < 0) continue;

        const labelList cp = mesh.cellPoints(c);
        forAll(cp, i)
        {
            if (mag(ptsNow[cp[i]] - points0[cp[i]]) > 1e-3*dzCell[c])
            {
                colMoved[col] = 1;
                break;
            }
        }
    }

    Info<< "SWENSE level correction: skipped moved columns = "
        << sum(colMoved) << " of " << nCol << endl;
}

// Global maximum, and the x of a column that reaches it (largest x on a tie
// between processors; serial keeps the first column)
{
    const scalar dMaxLocal = dMax;
    reduce(dMax, maxOp<scalar>());
    xDMax = (dMaxLocal == dMax) ? xDMax : -GREAT;
    reduce(xDMax, maxOp<scalar>());
}

// Per-column Delta for comparing decompositions (controlDict
// levelCorrectionColumnLog, off by default). First correction call only.
if (colLogOn && !colLogDone)
{
    colLogDone = true;

    List<stringList> keysAll(Pstream::nProcs());
    List<scalarList> deltaAll(Pstream::nProcs());
    keysAll[Pstream::myProcNo()] = colKeyList;
    deltaAll[Pstream::myProcNo()] = scalarList(Delta);
    Pstream::allGatherList(keysAll);
    Pstream::allGatherList(deltaAll);

    if (Pstream::master())
    {
        OFstream os(fileName("columnDelta.dat"));
        os.precision(17);
        os << "# column key Delta at t = " << runTime.timeName() << nl;
        forAll(keysAll, proc)
        {
            forAll(keysAll[proc], i)
            {
                os << keysAll[proc][i] << " " << deltaAll[proc][i] << nl;
            }
        }
    }
}

// Shifted alphaI and the increment
scalarField inc(mesh.nCells(), 0.0);
forAll(colOfCell, c)
{
    const label col = colOfCell[c];
    if (col < 0 || colMoved[col] == 1) continue;

    const scalar D = Delta[col];
    const scalar s = mag(D)/dzCell[c];
    scalar aLc = aI0[c];
    if (D > 0)
    {
        aLc = aI0[c] + s*(aI0[nbUp[c]] - aI0[c]);
    }
    else if (D < 0)
    {
        aLc = aI0[c] + s*(aI0[nbDown[c]] - aI0[c]);
    }
    inc[c] = aLc - aI0[c];
}

alphaI.primitiveFieldRef() += inc;

// alpha1 takes the same increment, clipped to [0, 1]; the clipped volume is
// reported (alphaI is not clipped)
const scalarField a1Pre(alpha1.primitiveField() + inc);
const scalarField a1Clip(min(max(a1Pre, scalar(0)), scalar(1)));
const scalarField clipCell(mag(a1Pre - a1Clip)*mesh.V());
const scalar clipVol = gSum(clipCell);
clipTotal += clipVol;
alpha1.primitiveFieldRef() = a1Clip;

alphaI.correctBoundaryConditions();
alpha1.correctBoundaryConditions();

// Delta per column after the correction (should be round-off)
const scalarField& aIn = alphaI.primitiveField();
scalarField dPost(nCol, 0.0);
forAll(colOfCell, c)
{
    if (colOfCell[c] >= 0)
    {
        dPost[colOfCell[c]] += aIn[c]*dzCell[c];
    }
}
scalar dPostMax = 0.0;
forAll(dPost, col)
{
    dPostMax = max
    (
        dPostMax,
        mag(zBotCol[col] + dPost[col] - incident.etaRow(rowOfCol[col], tNow))
    );
}
reduce(dPostMax, maxOp<scalar>());

Info<< "level correction: max|Delta| before = " << dMax << " m at x = " << xDMax
    << ", after = " << dPostMax << " m" << endl;
Info<< "level correction: clipped alpha1 volume = " << clipVol
    << " m3, running total = " << clipTotal << " m3" << endl;
