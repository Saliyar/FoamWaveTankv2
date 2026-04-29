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
    Foam::waveTheories::DepthvaryingCurrent


Description
    Depth-dependent staged current wave theory.

    Implements a piecewise current profile where the velocity is constant
    below a transition depth ZC and varies linearly above ZC up to the
    free surface. A soft-start ramp is applied to avoid numerical transients.


current
{
    type        DepthvaryingCurrent;

    U           (0.5 0 0);   // deep current
    UF          (1.0 0 0);   // surface current
    ZC          -2.0;        // transition depth (negative)
    Tsoft       5;
    localSeaLevel 0;
}

SourceFiles
    DepthvaryingCurrent.C

Author
   Sithik Aliyar, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/


#ifndef DepthvaryingCurrent_H
#define DepthvaryingCurrent_H

#include "waveTheory.H"

namespace Foam
{
namespace waveTheories
{

class DepthvaryingCurrent
:
    public waveTheory
{
    // Private data

        //- Deep current velocity (below ZC)
        vector U_;

        //- Surface current velocity
        vector UF_;

        //- Transition depth (negative)
        scalar ZC_;

        //- Soft-start time
        scalar Tsoft_;

        //- Local still water level
        scalar localSeaLevel_;

        //- Current profile type
        word profileType_;

        //- Parabolic profile parameters
        scalar H_;
        scalar z0_;


public:

    //- Runtime type information
    TypeName("DepthvaryingCurrent");

    // Constructors

        DepthvaryingCurrent
        (
            const word& subDictName,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~DepthvaryingCurrent() = default;


    // Member Functions

        //- Print coefficients
        virtual void printCoeffs();

        //- Soft-start factor
        scalar factor(const scalar& time) const;

        //- Free-surface elevation
        virtual scalar eta
        (
            const point& x,
            const scalar& time
        ) const;

        //- Excess pressure
        virtual scalar pExcess
        (
            const point& x,
            const scalar& time
        ) const;

        //- Velocity field (USED BY RELAXATION ZONES)
        virtual vector U
        (
            const point& x,
            const scalar& time
        ) const;

        //- TSoft conflict
        virtual Switch conflictTSoftInitialise() const;
};

} // End namespace waveTheories
} // End namespace Foam

#endif

