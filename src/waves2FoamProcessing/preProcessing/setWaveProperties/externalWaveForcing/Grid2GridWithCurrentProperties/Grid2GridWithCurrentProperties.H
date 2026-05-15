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

Class
    Foam::Grid2GridProperties

Description


SourceFiles
    Grid2GridProperties.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Grid2GridWithCurrentProperties_H
#define Grid2GridWithCurrentProperties_H

#include "setWaveProperties.H"

namespace Foam
{

class Grid2GridWithCurrentProperties
:
    public setWaveProperties
{
    // No additional private data needed
    // (All handled via dictionary)

public:

    //- Runtime type information
    TypeName("Grid2GridWithCurrentProperties");


    // Constructors
    Grid2GridWithCurrentProperties
    (
        const Time& rT,
        dictionary& dict,
        vector g,
        bool write
    );


    // Destructor
    virtual ~Grid2GridWithCurrentProperties() = default;


    // Member Functions

    //- Set wave properties into dictionary
    virtual void set(Ostream& os);
};

} // End namespace Foam

#endif

// ************************************************************************* //
// ************************************************************************* //