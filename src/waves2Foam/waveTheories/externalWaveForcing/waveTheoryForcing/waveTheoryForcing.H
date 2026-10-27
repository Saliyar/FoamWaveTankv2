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
    Foam::waveTheories::waveTheoryForcing
Description
    External wave forcing that evaluates an analytic wave theory (for example
    streamFunction) in place of Grid2Grid. The theory is selected in
    waveProperties by waveTheoryForcingCoeffs/waveType.

    Pressure convention: swenseIncident uses
        pW = rhoWater (pExcess - g z)
    so pExcess is per unit water density. The theory returns pressure in Pa
    (rho included, hydrostatic included), so pExcess = theory.pExcess / rho.

    Requires explicit "Tsoft 0;" in waveTheoryForcingCoeffs: a start-up ramp
    is not a solution of the free-surface conditions.

SourceFiles
    waveTheoryForcing.C

\*---------------------------------------------------------------------------*/

#ifndef waveTheoryForcing_H
#define waveTheoryForcing_H

#include "externalWaveForcing.H"
#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                       Class waveTheoryForcing Declaration
\*---------------------------------------------------------------------------*/

class waveTheoryForcing
:
    public externalWaveForcing
{
protected:

    // Protected data

        //- The analytic wave theory, built from waveTheoryForcingCoeffs
        autoPtr<waveTheory> theory_;

        //- Water density [kg/m3], copied from the theory
        const scalar rhoWater_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        waveTheoryForcing(const waveTheoryForcing&);

        //- Disallow default bitwise assignment
        void operator=(const waveTheoryForcing&);


public:

    //- Runtime type information
    TypeName("waveTheoryForcing");


    // Constructors

        //- Construct from components
        waveTheoryForcing
        (
            IOobject io,
            Time& rT,
            const fvMesh& mesh
        );


    // Destructor

        virtual ~waveTheoryForcing()
        {}


    // Member Functions

        //- Nothing to be done: the theory is analytic
        virtual void step();

        //- Nothing to be done
        virtual void close();

        //- Free-surface elevation of the theory [m]
        virtual scalar eta
        (
            const point& x,
            const scalar& time
        ) const;

        //- Excess pressure per unit water density [m2/s2]
        virtual scalar pExcess
        (
            const point& x,
            const scalar& time
        ) const;

        //- Pressure is implemented (checked in the constructor)
        virtual bool implementPressure() const
        {
            return true;
        };

        //- Water velocity of the theory [m/s]
        virtual vector U
        (
            const point& x,
            const scalar& time
        ) const;

        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
