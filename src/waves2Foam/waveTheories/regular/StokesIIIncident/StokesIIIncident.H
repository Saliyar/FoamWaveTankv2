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
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::waveTheories::StokesIIIncident

Description
    Second-order Stokes wave as an incident field everywhere in the domain,
    the same wave as the ESI waveModels boundary condition StokesII of the
    plain case.

    Inputs are read from the ESI wave dictionary of the plain case
    (waveProperties subDict named by waveDict, usually the inlet patch):
    waveType (must be StokesII), waveHeight, wavePeriod, waveAngle [deg],
    wavePhase (default 1.5 pi, as ESI). The water depth is seaLevel minus the
    lowest point of that patch; a waterDepthRef entry must agree with it.
    rampTime is not used: the incident has no start-up ramp.

    Formulas, from OpenFOAM-v2506 (GPL):
    - eta: src/waveModels/waveGenerationModels/derived/StokesII/
      StokesIIWaveModel.C:54-74
    - velocity: StokesIIWaveModel.C:79-105, z measured from the bottom
    - wave length: .../derived/StokesI/StokesIWaveModel.C:74-89
    - patch frame (x = -mean patch normal, z = -g) and the phase reference:
      src/waveModels/waveModel/waveModel.C:57-120, 402; regularWaveModel.C:80
    ESI has no pressure. The pressure is the second-order Stokes pressure
    (Dean and Dalrymple 1991), the same expression as
    regular/stokesSecond/stokesSecond.C pExcess in this library.

    waveTheoryForcingCoeffs
    {
        waveType    StokesIIIncident;
        waveDict    inlet;
        Tsoft       0;
    }

SourceFiles
    StokesIIIncident.C

\*---------------------------------------------------------------------------*/

#ifndef StokesIIIncident_H
#define StokesIIIncident_H

#include "waveTheory.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                       Class StokesIIIncident Declaration
\*---------------------------------------------------------------------------*/

class StokesIIIncident
:
    public waveTheory
{
protected:

    // Protected data

        //- Wave height H [m], period T [s], angle [rad], phase [rad]
        scalar H_;
        scalar T_;
        scalar angle_;
        scalar phase_;

        //- Patch frame: global = Rlg & local, local = Rgl & global
        tensor Rlg_;
        tensor Rgl_;

        //- Lowest local z of the wave patch [m]
        scalar zBottom_;

        //- Depth h [m], wave length L [m], wave number K [1/m],
        //  angular frequency omega [1/s]
        scalar h_;
        scalar L_;
        scalar K_;
        scalar omega_;

        //- Factor on the first-order velocity (switch fluxConsistentU)
        scalar fU_;


    // Protected member functions

        //- No ramp
        scalar factor(const scalar&) const
        {
            return 1.0;
        }

        //- Phase K_x x + K_y y - omega t + phase in the patch frame
        scalar phaseTot(const point& xl, const scalar time) const;


public:

    //- Runtime type information
    TypeName("StokesIIIncident");


    // Constructors

        StokesIIIncident
        (
            const word& subDictName,
            const fvMesh& mesh_
        );


    // Destructor

        ~StokesIIIncident()
        {}


    // Member Functions

        virtual Switch conflictTSoftInitialise() const
        {
            return false;
        }

        void printCoeffs();

        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;

        //- Dynamic pressure rho (p/rho + g z) [Pa], as stokesSecond
        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;

        bool implementPressure() const
        {
            return true;
        }

        virtual scalar phaseSpeed() const
        {
            return omega_/K_;
        }

        virtual scalar waveNumber() const
        {
            return K_;
        }

        virtual vector U
        (
            const point&,
            const scalar&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
