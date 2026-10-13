/*---------------------------------------------------------------------------*\
Class
    Foam::swenseIncident

Description
    Incident (wave) fields for the SWENSE solver.

    Each cell and face is classified by its centre: water if z < eta(x, y, t),
    air otherwise. Water values come from the Grid2Grid external wave forcing
    (waves2Foam), air values from the Air2Air solver.

    Wave frame (Grid2Grid convention, from Grid2GridCoeffs):
        x_w = R(theta) (x - referencePosition),  R(theta) = [c -s; s c]
    The Air2Air time is the Grid2Grid time t + referenceTime.

    Water:
        UInc = U_grid2grid(x, t)
        pInc = rho_w (pExcess(x, t) - g_mag z_w)        [hydrostatic included]
    Air:
        UInc = R^T U_air(x_w, y_w, z - seaLevel)
        pInc = p_air(x_w, y_w, z - seaLevel)

    Cell alphaInc is the wetted volume fraction (swenseSurface).

    Interface correction: on internal faces whose two neighbouring cells have
    different classification, pfInc is the water-side pressure at the free
    surface above the face centre.

SourceFiles
    swenseIncident.C

\*---------------------------------------------------------------------------*/

#ifndef swenseIncident_H
#define swenseIncident_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "autoPtr.H"
#include "externalWaveForcing.H"
#include "Air2AirWrapper.H"
#include "swenseSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swenseIncident Declaration
\*---------------------------------------------------------------------------*/

class swenseIncident
{
    // Private data

        const fvMesh& mesh_;

        //- Grid2Grid external wave forcing (registered in the mesh database)
        const waveTheories::externalWaveForcing& ext_;

        //- Air2Air solver (optional, required for air values)
        autoPtr<Air2AirWrapper> air_;

        //- Wetted-geometry helper
        swenseSurface surf_;

        //- Wave frame from Grid2GridCoeffs (same frame as the Grid2Grid solver)
        const vector refPos_;
        const scalar theta_;
        const scalar refTime_;

        //- Mean water level and water density
        const scalar seaLevel_;
        const scalar rhoWater_;

        //- Gravity magnitude
        const scalar magG_;

        //- Incident fields
        volVectorField UInc_;
        volScalarField pInc_;
        volScalarField alphaInc_;
        surfaceVectorField UfInc_;
        surfaceScalarField pfInc_;

        //- Cell classification from the last update: 1 water, 0 air
        scalarField cellWater_;

    // Private member functions

        //- Rotation from the mesh frame to the wave frame (about z)
        vector toWave(const vector& x) const;

        //- Rotation from the wave frame back to the mesh frame (about z)
        vector toMesh(const vector& u) const;

        //- Water velocity and pressure (mesh frame)
        vector waterU(const point& x, const scalar t) const;
        scalar waterP(const point& x, const scalar t) const;

        //- Air velocity and pressure (mesh frame)
        vector airU(const point& x) const;
        scalar airP(const point& x) const;

        //- Water-side pressure at the free surface above x
        scalar interfaceP(const point& x, const scalar t) const;

        //- Disallow default bitwise copy
        swenseIncident(const swenseIncident&);
        void operator=(const swenseIncident&);


public:

    // Constructors

        //- Construct from mesh, the registered external wave forcing and the
        //  water density. Reads Grid2GridCoeffs and the optional Air2Air block
        //  of farfieldCoeffs from waveProperties.
        swenseIncident
        (
            const fvMesh& mesh,
            const waveTheories::externalWaveForcing& ext,
            const dimensionedScalar& rhoWater
        );


    // Member Functions

        //- Update all incident fields for simulation time t
        void update(const scalar t);

        //- Accessors
        const volVectorField& UInc() const { return UInc_; }
        const volScalarField& pInc() const { return pInc_; }
        const volScalarField& alphaInc() const { return alphaInc_; }
        const surfaceVectorField& UfInc() const { return UfInc_; }
        const surfaceScalarField& pfInc() const { return pfInc_; }
        const scalarField& cellWater() const { return cellWater_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
