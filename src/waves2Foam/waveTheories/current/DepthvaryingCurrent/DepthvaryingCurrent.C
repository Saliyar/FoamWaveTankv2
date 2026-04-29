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


#include "DepthvaryingCurrent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

defineTypeNameAndDebug(DepthvaryingCurrent, 0);
addToRunTimeSelectionTable(waveTheory, DepthvaryingCurrent, dictionary);

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

DepthvaryingCurrent::DepthvaryingCurrent
(
    const word& subDictName,
    const fvMesh& mesh
)
:
    waveTheory(subDictName, mesh),
    U_(vector(coeffDict_.lookup("U"))),
    UF_(vector(coeffDict_.lookup("UF"))),
    ZC_(readScalar(coeffDict_.lookup("ZC"))),
    Tsoft_(readScalar(coeffDict_.lookup("Tsoft"))),
    localSeaLevel_
    (
        coeffDict_.lookupOrDefault<scalar>("localSeaLevel", seaLevel_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DepthvaryingCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << nl
         << "  Deep current U        = " << U_ << nl
         << "  Surface current UF   = " << UF_ << nl
         << "  Transition depth ZC  = " << ZC_ << nl
         << "  Soft start Tsoft     = " << Tsoft_ << endl;
}


scalar DepthvaryingCurrent::factor(const scalar& time) const
{
    scalar f(1.0);

    if (Tsoft_ > 0.0)
    {
        f = sin(PI_/2.0 * min(Tsoft_, time) / Tsoft_);
    }

    return f;
}


scalar DepthvaryingCurrent::eta
(
    const point&,
    const scalar&
) const
{
    return localSeaLevel_;
}


scalar DepthvaryingCurrent::pExcess
(
    const point&,
    const scalar&
) const
{
    return referencePressure(localSeaLevel_);
}


vector DepthvaryingCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    scalar z = x.z();
    vector UC;

    if (z < ZC_)
    {
        UC = UF_;
    }
    else
    {
        UC =
            (UF_ - U_) / (-ZC_) * z
          + UF_;
    }

    return UC * factor(time);
}

Switch DepthvaryingCurrent::conflictTSoftInitialise() const
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //

