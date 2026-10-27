/*---------------------------------------------------------------------------*\
Class
    Foam::swenseSurface

Description
    Cut-cell and cut-face geometry against the free surface of a wave theory.

    Uses the convexPolyhedral machinery of waves2Foam with the signed distance
        d(x) = z - eta(x, y, t)
    so that d < 0 is water (below the free surface) and d > 0 is air.
    The wetted volume of a cell is the volume of its negative part, and the
    wetted area of a face is the area of its negative part.

SourceFiles
    swenseSurface.C

\*---------------------------------------------------------------------------*/

#ifndef swenseSurface_H
#define swenseSurface_H

#include "convexPolyhedral.H"
#include "externalWaveForcing.H"
#include "swenseRamp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swenseSurface Declaration
\*---------------------------------------------------------------------------*/

class swenseSurface
:
    public convexPolyhedral
{
    // Private data

        //- External (Grid2Grid) wave forcing that gives eta(x, y, t)
        const waveTheories::externalWaveForcing& ext_;

        //- Time at which the surface is evaluated
        scalar t_;

        //- Start-up ramp time (rampTsoft in waveProperties, 0 = none)
        const scalar Tsoft_;

    // Protected member functions (overrides of convexPolyhedral)

        virtual void signedPointToSurfaceDistance
        (
            const pointField& pp,
            scalarField& sd
        );

        virtual scalar signedPointToSurfaceDistance
        (
            const point& pp
        ) const;


public:

    // Constructors

        //- Construct from mesh, external wave forcing and time
        swenseSurface
        (
            const fvMesh& mesh,
            const waveTheories::externalWaveForcing& ext,
            const scalar t
        );


    // Member Functions

        //- Set the time at which the surface is evaluated
        void setTime(const scalar t);

        //- Wetted (water-side) volume of cell celli
        scalar wetVolume(const label celli);

        //- Wetted (water-side) area of face facei
        scalar wetArea(const label facei);

        //- Centre of the wetted part of face facei
        vector wetCentre(const label facei);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
