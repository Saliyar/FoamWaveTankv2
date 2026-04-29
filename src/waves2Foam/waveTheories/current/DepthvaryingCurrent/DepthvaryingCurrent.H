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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{
/*---------------------------------------------------------------------------*\
                           Class DepthvaryingCurrent Declaration
\*---------------------------------------------------------------------------*/

class DepthvaryingCurrent
: 
    public waveTheory
{
    //- Deep current velocity
    vector U_;

    //- Surface current velocity
    vector UF_;

    //- Transition depth
    scalar ZC_;

    //- Soft start time
    scalar Tsoft_;

    //- Local sea level
    scalar localSeaLevel_;

    //- Ramp factor
    scalar factor(const scalar& time) const;

public:
    TypeName("DepthvaryingCurrent");

    DepthvaryingCurrent
    (
        const word& subDictName,
        const fvMesh& mesh
    );

    virtual void printCoeffs();

    virtual scalar eta
    (
        const point& x,
        const scalar& time
    ) const;

    virtual vector U
    (
        const point& x,
        const scalar& time
    ) const;

    virtual scalar pExcess
    (
        const point& x,
        const scalar& time
    ) const;

    virtual Switch conflictTSoftInitialise() const;
};

} // namespace waveTheories
} // namespace Foam

#endif
