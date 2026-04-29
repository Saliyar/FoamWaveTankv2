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


#include "DepthvaryingCurrentProperties.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(DepthvaryingCurrentProperties, 0);

addToRunTimeSelectionTable
(
    setWaveProperties,
    DepthvaryingCurrentProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

DepthvaryingCurrentProperties::DepthvaryingCurrentProperties
(
    const Time& rT,
    dictionary& dict,
    vector g,
    bool write
)
:
    setWaveProperties(rT, dict, g, write)
{
    Info<< "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void DepthvaryingCurrentProperties::set(Ostream& os)
{
    // Write opening of wave sub-dictionary
    writeBeginning(os);

    // Write required entries
    writeGiven(os, "waveType");
    writeGiven(os, "profileType");
    writeGiven(os, "U");
    writeGiven(os, "UF");
    writeGiven(os, "ZC");
    writeGiven(os, "Tsoft");

    if (dict_.found("localSeaLevel"))
    {
        writeGiven(os, "localSeaLevel");
    }

    
// ---- Validate profileType presence ------------------
    if (!dict_.found("profileType"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing mandatory entry 'profileType' in "
            << "DepthvaryingCurrentProperties." << nl
            << "Valid options are:" << nl
            << "  - constant" << nl
            << "  - depthVarying" << nl
            << "  - parabolic"
            << exit(FatalIOError);
    }


    // Write relaxation zone (inlet / outlet)
    writeRelaxationZone(os);

    // Write closing bracket
    writeEnding(os);
}

} // End namespace Foam

// ************************************************************************* //
