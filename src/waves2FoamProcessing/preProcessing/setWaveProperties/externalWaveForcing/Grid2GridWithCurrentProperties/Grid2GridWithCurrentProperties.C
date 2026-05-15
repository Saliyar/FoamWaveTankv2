/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Grid2GridWithCurrentProperties.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(Grid2GridWithCurrentProperties, 0);

addToRunTimeSelectionTable
(
    setWaveProperties,
    Grid2GridWithCurrentProperties,
    setWaveProperties
);


// * * * * * * * Constructors * * * * * * * //

Grid2GridWithCurrentProperties::Grid2GridWithCurrentProperties
(
    const Time& rT,
    dictionary& dict,
    vector g,
    bool write
)
:
    setWaveProperties(rT, dict, g, write)
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * Member Functions * * * * * * * //

void Grid2GridWithCurrentProperties::set(Ostream& os)
{
    // Begin block
    writeBeginning(os);

    writeGiven(os, "waveType");
    os << endl;

    // ======================================
    // ✅ ORIGINAL GRID2GRID PARAMETERS
    // ======================================
    scalar Tsoft
    (
        dict_.lookupOrDefault<scalar>("Tsoft", 0.0)
    );

    scalar referenceTime
    (
        dict_.lookupOrDefault<scalar>("referenceTime", 0.0)
    );

    vector referencePosition
    (
        dict_.lookupOrDefault<vector>("referencePosition", vector::zero)
    );

    scalar referenceTheta
    (
        dict_.lookupOrDefault<scalar>("referenceTheta", 0.0)
    );

    writeDerived(os, "Tsoft", Tsoft);
    os << endl;

    writeDerived(os, "referenceTime", referenceTime);
    os << endl;

    writeDerived(os, "referencePosition", referencePosition);
    os << endl;

    writeDerived(os, "referenceTheta", referenceTheta);
    os << endl;

    writeGiven(os, "Grid2GridInputFileName");
    os << endl;


    // ======================================
    // ✅ NEW: CURRENT PARAMETERS
    // ======================================

    // Validate profileType
    if (!dict_.found("profileType"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing mandatory entry 'profileType' in Grid2GridWithCurrent." << nl
            << "Valid options are:" << nl
            << "  - constant" << nl
            << "  - depthVarying" << nl
            << "  - parabolic"
            << exit(FatalIOError);
    }

    writeGiven(os, "profileType");

    // Deep current
    writeGiven(os, "Uc");

    // Surface current (optional for constant, but safe to pass)
    if (dict_.found("UF"))
    {
        writeGiven(os, "UF");
    }

    // Transition depth
    if (dict_.found("ZC"))
    {
        writeGiven(os, "ZC");
    }

    // Current ramp
    scalar currentTsoft
    (
        dict_.lookupOrDefault<scalar>("currentTsoft", 0.0)
    );

    writeDerived(os, "currentTsoft", currentTsoft);

    // Optional parabolic inputs
    if (dict_.found("H"))
    {
        writeGiven(os, "H");
    }

    if (dict_.found("z0"))
    {
        writeGiven(os, "z0");
    }


    // ======================================
    // ✅ RELAXATION ZONE (IMPORTANT)
    // ======================================
    writeRelaxationZone(os);

    // End block
    writeEnding(os);
}

} // End namespace Foam


// ************************************************************************* //