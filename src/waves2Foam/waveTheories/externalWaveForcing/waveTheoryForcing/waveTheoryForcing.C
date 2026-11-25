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

#include "waveTheoryForcing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveTheoryForcing, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    waveTheoryForcing,
    externalWaveForcing
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waveTheoryForcing::waveTheoryForcing
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh),
    theory_(waveTheory::New("waveTheoryForcing", mesh)),
    rhoWater_(theory_->rhoWater())
{
    if (!theory_->implementPressure())
    {
        FatalErrorInFunction
            << "The selected wave theory does not implement the pressure,"
            << " so it cannot be used as waveTheoryForcing"
            << " (swenseIncident needs pExcess)."
            << exit(FatalError);
    }

    const scalar Tsoft =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .subDict("waveTheoryForcingCoeffs")
            .lookupOrDefault<scalar>("Tsoft", -1.0);

    if (mag(Tsoft) > SMALL)
    {
        FatalErrorInFunction
            << "waveTheoryForcing requires Tsoft 0 in waveTheoryForcingCoeffs"
            << " (found " << Tsoft << " or none). A start-up ramp is not a"
            << " solution of the free-surface conditions." << nl
            << "    add 'Tsoft 0;' to waveTheoryForcingCoeffs"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void waveTheoryForcing::step()
{
    // Nothing to be done
}


void waveTheoryForcing::close()
{
    // Nothing to be done
}


scalar waveTheoryForcing::eta
(
    const point& x,
    const scalar& time
) const
{
    return theory_->eta(x, time);
}


scalar waveTheoryForcing::pExcess
(
    const point& x,
    const scalar& time
) const
{
    // Theory returns rho (p/rho + g z); swenseIncident wants p/rho + g z
    return theory_->pExcess(x, time)/rhoWater_;
}


vector waveTheoryForcing::U
(
    const point& x,
    const scalar& time
) const
{
    return theory_->U(x, time);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
