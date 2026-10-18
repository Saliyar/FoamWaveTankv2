// Level correction of the incident surface (SWENSE, document Rev1 16.3.3).
//
// Per column (vertical stack of cells over one (x, y)):
//   Delta = z_bottom + sum(alphaI dz) - eta_H,  eta_H the Grid2Grid eta of the
//   surface row of the column at the current time,
//   alphaIL = alphaI shifted by s = |Delta|/dz_c toward the neighbour above
//   (Delta > 0) or below (Delta < 0); top and bottom cells use themselves,
//   inc = alphaIL - alphaI, added to alphaI and to alpha1.
// The relaxation densities are computed after this, in swenseAlphaRelax.H.

static bool colInit = true;
static labelList colOfCell;
static labelList nbUp;
static labelList nbDown;
static scalarField dzCell;
static scalarField zBotCol;
static labelList rowOfCol;
static scalarField colX, colY;
static label nCol = 0;

if (colInit)
{
    colInit = false;

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "column level correction needs every column on one processor; "
            << "parallel runs are not supported"
            << exit(FatalError);
    }

    const vectorField& Cc = mesh.C().primitiveField();
    const pointField& pts = mesh.points();

    // Columns: cells with the same (x, y)
    HashTable<label, string> colMap;
    DynamicList<point> colXY;
    colOfCell.setSize(mesh.nCells());
    forAll(Cc, c)
    {
        const string key
        (
            name(lround(Cc[c].x()*1e6)) + "_" + name(lround(Cc[c].y()*1e6))
        );
        if (!colMap.found(key))
        {
            colMap.insert(key, colMap.size());
            colXY.append(Cc[c]);
        }
        colOfCell[c] = colMap[key];
    }
    nCol = colMap.size();

    colX.setSize(nCol);
    colY.setSize(nCol);
    forAll(colXY, col)
    {
        colX[col] = colXY[col].x();
        colY[col] = colXY[col].y();
    }

    // Cell heights and column bottoms
    dzCell.setSize(mesh.nCells());
    zBotCol.setSize(nCol, GREAT);
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
        zBotCol[colOfCell[c]] = min(zBotCol[colOfCell[c]], zmin);
    }

    // Cells of each column sorted by height: neighbours above and below
    List<DynamicList<label>> colCells(nCol);
    forAll(colOfCell, c)
    {
        colCells[colOfCell[c]].append(c);
    }
    nbUp.setSize(mesh.nCells());
    nbDown.setSize(mesh.nCells());
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
}

const scalarField aI0(alphaI.primitiveField());
const scalar tNow = mesh.time().value();

// Delta per column, before the correction
scalarField colSum(nCol, 0.0);
forAll(colOfCell, c)
{
    colSum[colOfCell[c]] += aI0[c]*dzCell[c];
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

// Shifted alphaI and the increment
scalarField inc(mesh.nCells(), 0.0);
forAll(colOfCell, c)
{
    const label col = colOfCell[c];
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
alpha1.primitiveFieldRef() += inc;
alphaI.correctBoundaryConditions();
alpha1.correctBoundaryConditions();

// Delta per column after the correction (should be round-off)
const scalarField& aIn = alphaI.primitiveField();
scalarField dPost(nCol, 0.0);
forAll(colOfCell, c)
{
    dPost[colOfCell[c]] += aIn[c]*dzCell[c];
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

Info<< "level correction: max|Delta| before = " << dMax << " m at x = " << xDMax
    << ", after = " << dPostMax << " m" << endl;
