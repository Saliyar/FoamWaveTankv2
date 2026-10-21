#!/usr/bin/env python3
"""Whole-column decomposition for the SWENSE column level correction.

Reads the cell centres written by postProcess -func writeCellCentres (the C
file of a time directory), groups cells into columns (same x and y, rounded
to 1e-6), sorts the columns by x and gives each rank a block of consecutive
columns. Every cell takes the rank of its column, so no column is split.

Writes an OpenFOAM labelList for decomposePar with method manual:

    python3 swenseColumnDecomposition.py <time>/C <nRanks> <outFile>

Then in system/decomposeParDict:

    numberOfSubdomains <nRanks>;
    method manual;
    manualCoeffs { dataFile "cellDecomposition"; }

with the output file in constant/ (manual reads it from there).
"""

import sys
from collections import defaultdict


def readCellCentres(path):
    """Cell centres from an OpenFOAM vectorField file (ascii)."""
    with open(path) as f:
        lines = f.read().split("\n")

    start = next(i for i, l in enumerate(lines) if l.strip() == "(")

    pts = []
    for l in lines[start + 1:]:
        s = l.strip()
        if s.startswith(")"):
            break
        if s.startswith("("):
            pts.append(tuple(float(v) for v in s.strip("()").split()))
    return pts


def columnBlocks(nCol, nRanks):
    """Rank boundaries over the sorted column index (sizes differ by <= 1)."""
    return [round(nCol*r/nRanks) for r in range(nRanks + 1)]


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__)

    cFile, nRanks, outFile = argv[1], int(argv[2]), argv[3]

    pts = readCellCentres(cFile)
    n = len(pts)

    key = [(round(p[0]*1e6), round(p[1]*1e6)) for p in pts]

    counts = defaultdict(int)
    for k in key:
        counts[k] += 1

    cols = sorted(set(key))
    colIdx = {k: i for i, k in enumerate(cols)}
    nCol = len(cols)

    if nRanks > nCol:
        sys.exit("more ranks (%d) than columns (%d)" % (nRanks, nCol))

    print("cells %d, columns %d, cells per column min %d max %d"
          % (n, nCol, min(counts.values()), max(counts.values())))

    bounds = columnBlocks(nCol, nRanks)

    def rankOfCol(i):
        for r in range(nRanks):
            if bounds[r] <= i < bounds[r + 1]:
                return r

    rankOfCell = [rankOfCol(colIdx[k]) for k in key]

    print("columns per rank", [bounds[r + 1] - bounds[r] for r in range(nRanks)])
    print("cells per rank", [rankOfCell.count(r) for r in range(nRanks)])

    with open(outFile, "w") as f:
        f.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
                "    class       labelList;\n    object      cellDecomposition;\n}\n\n")
        f.write("%d\n(\n" % n)
        for r in rankOfCell:
            f.write("%d\n" % r)
        f.write(")\n")

    print("wrote", outFile)


if __name__ == "__main__":
    main(sys.argv)
