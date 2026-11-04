/*---------------------------------------------------------------------------*\
Class
    Foam::swenseIncident

Description
    Blended incident (wave) fields for the SWENSE solver.

    Water and air incident fields are blended with the incident volume
    fraction alphaI (document swense_blended_incident, Section 6):

        alphaI = clamp(0.5 - (z - eta)/delta, 0, 1)
        UInc   = alphaI UW + (1 - alphaI) UA
        pInc   = alphaI pW + (1 - alphaI) pA
        rhoI   = alphaI rhoW + (1 - alphaI) rhoA

    The blending law is only a starting weight: alphaI is then advected with
    the incident velocity by the solver (Section 14.2), so the blend weight
    returned by blendWeight() is used only for its initial value and for the
    farfield boundary condition.

    Water (alphaI = 1 region), Grid2Grid, in the wave frame:
        UW = U_grid2grid(x, t)
        pW = rho_w (pExcess(x, t) - g_mag z_w)        [hydrostatic included]
    Air (alphaI = 0 region), Air2Air:
        UA = R^T U_air(x_w, y_w, z - seaLevel)
        pA = p_air(x_w, y_w, z - seaLevel) - rho_a g z_w   [hydrostatic]

    Wave frame (Grid2Grid convention, from Grid2GridCoeffs):
        x_w = R(theta) (x - referencePosition),  R(theta) = [c -s; s c]
    The Air2Air time is the Grid2Grid time t + referenceTime.

    The water or air value is only evaluated where its weight is non-zero, so
    Air2Air is never called deep in the water and Grid2Grid is never called
    deep in the air.

    Blending thickness delta is read from waveProperties (blendThickness,
    default 0.024 m).

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
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "Grid2Grid.H"
#include "swenseRamp.H"

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

        //- Wave frame from Grid2GridCoeffs (same frame as the Grid2Grid solver)
        const vector refPos_;
        const scalar theta_;
        const scalar refTime_;

        //- Mean water level and water density
        const scalar seaLevel_;
        const scalar rhoWater_;

        //- Start-up ramp time (rampTsoft, 0 = none) and the simulation time
        //  of the last update (the Air2Air time is shifted by referenceTime)
        const scalar Tsoft_;
        scalar tNow_;

        //- Start-up ramp factor f(t) of all incident values
        scalar ramp(const scalar t) const;

        //- Blending thickness delta (m)
        const scalar delta_;

        //- Momentum-weighted blend and linear absolute pressure blend
        //  (waveProperties momentumBlend, default false: current scheme)
        const bool momentumBlend_;

        //- Guard distance d_g = max(delta, 3 dz), dz = cell height at the
        //  interface (document Rev4, 16.3 rule 5)
        scalar guard_;

        // Point tables (see buildPointTables)
        pointField pts_;
        label offCb_, offCf_, offCfb_;
        labelList etaIdx_;          // per flat point: row of its distinct (x, y)
        pointField etaPts_;         // per row: (x, y, 0)
        scalarField etaVal_;        // per row: eta at time etaTime_
        scalar etaTime_;
        labelList airIx_, airIy_, airIz_;   // per flat point: Air2Air table rows

        //- [option A] Grid2Grid wrapper (null unless the external forcing is Grid2Grid)
        const waveTheories::Grid2Grid* g2g_;

        //- [option A] Grid2Grid values at all flat points for the current step
        bool fastEval_;
        scalarField etaPt_, pExPt_;
        // Scratch fix: per-row surface pressures (water, air) at (x, y, eta_r),
        // and one flat point per row to read eta_r
        scalarField pWs_, pAs_;
        labelList rowPt_;
        vectorField UPt_;

        //- Air density (Air2Air pressure is linear in it; hydrostatic -rho_a g z)
        const scalar rhoAir_;

        //- Gravity magnitude
        const scalar magG_;

        //- Blended incident fields
        volVectorField UInc_;
        volScalarField pInc_;
        surfaceVectorField UfInc_;
        surfaceScalarField pfInc_;
        surfaceScalarField psiWf_;

    // Private member functions

        //- Rotation from the mesh frame to the wave frame (about z)
        vector toWave(const vector& x) const;

        //- Rotation from the wave frame back to the mesh frame (about z)
        vector toMesh(const vector& u) const;

        //- Water velocity and pressure (mesh frame)
        vector waterU(const point& x, const scalar t) const;
        scalar waterP(const point& x, const scalar t) const;

        //- psiWf value at flat point ip (scratch)
        scalar psiWValue(const label ip, const scalar t) const;

        //- Air velocity (mesh frame) and pressure at flat point ip, from the
        //  Air2Air tables (one evaluation for both)
        vector airUP(const label ip, scalar& p) const;

        //- Weights of the water and air parts at x for the advected alphaI
        //  (document Rev3 16.3.5): W = clamp((alphaI - eps)/(1 - 2 eps), 0, 1),
        //  water part only where W > 0 and z - eta <= delta, air part only
        //  where W < 1 and z - eta >= -delta; a part that is not evaluated has
        //  weight exactly 0 and the other part weight exactly 1
        void weights
        (
            const point& x,
            const scalar eta,
            const scalar alphaI,
            scalar& wWater,
            scalar& wAir
        ) const;

        //- Grid2Grid eta at flat point ip: the per-step value when it is
        //  current for time t, otherwise evaluated here
        scalar etaAt(const label ip, const scalar t) const;

        //- Blended velocity u and pressure p at flat point ip for the advected
        //  alphaI (weights evaluated once for both)
        void blendUP
        (
            const label ip,
            const scalar t,
            const scalar alphaI,
            vector& u,
            scalar& p
        ) const;

        //- Point tables of the mesh: flat order cells, boundary cells, faces,
        //  boundary faces (the order of the update loops), distinct (x, y) for
        //  eta and distinct x, y, z for the Air2Air tables
        void buildPointTables();

        //- Disallow default bitwise copy
        swenseIncident(const swenseIncident&);
        void operator=(const swenseIncident&);


public:

    // Static Member Functions

        //- Blend weight of the water incident field for a height zRel above
        //  the free surface (zRel = z - eta): 1 deep in water, 0 in air
        static scalar blendWeight(const scalar zRel, const scalar delta);

        //- Blending thickness from waveProperties (default 0.024 m)
        static scalar blendThickness(const objectRegistry& db);

        //- Cell height at the interface: largest spacing of the cell-centre
        //  z levels with |z| < 0.1 m
        static scalar interfaceCellHeight(const fvMesh& mesh);


    // Constructors

        //- Construct from mesh, the registered external wave forcing, and the
        //  water and air densities. Reads Grid2GridCoeffs and the optional
        //  Air2Air block of farfieldCoeffs from waveProperties. The Air2Air
        //  rhoAir must equal the air density.
        swenseIncident
        (
            const fvMesh& mesh,
            const waveTheories::externalWaveForcing& ext,
            const dimensionedScalar& rhoWater,
            const dimensionedScalar& rhoAir
        );


    // Member Functions

        //- HOS target alphaH = clamp(0.5 - (z - eta)/delta, 0, 1) at cell
        //  centres for time t (relaxation target, initial alphaI)
        volScalarField blendAlpha(const scalar t) const;

        //- Update all blended incident fields for time t, blended with the
        //  advected incident fraction alphaI (cells, then faces with alphaIf
        //  = linear interpolation of alphaI)
        void update(const scalar t, const volScalarField& alphaI);

        //- Rebuild the evaluation points from the current cell and face
        //  centres after a mesh motion (call when mesh.moving()). Only the
        //  point lists and their Grid2Grid basis change; the Grid2Grid and
        //  Air2Air objects are not modified.
        void movePoints();

        //- Accessors
        const volVectorField& UInc() const { return UInc_; }
        //  Non-const access for the isoAdvector of alphaI, which holds a
        //  reference to UInc: only to make UInc relative to the mesh around
        //  its advection step, and restore it after
        volVectorField& UIncRW() { return UInc_; }
        const volScalarField& pInc() const { return pInc_; }
        const surfaceVectorField& UfInc() const { return UfInc_; }
        const surfaceScalarField& pfInc() const { return pfInc_; }

        //- Water pressure relative to the row surface value, per m2/s2 (scratch)
        const surfaceScalarField& psiWf() const { return psiWf_; }

        //- Surface rows (x, y, 0) and Grid2Grid eta of row r at time t
        const pointField& etaPts() const { return etaPts_; }
        scalar etaRow(const label r, const scalar t) const { return etaAt(rowPt_[r], t); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
