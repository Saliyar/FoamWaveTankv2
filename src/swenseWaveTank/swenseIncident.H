/*---------------------------------------------------------------------------*\
Class
    Foam::swenseIncident

Description
    Blended incident (wave) fields for the SWENSE solver.

    Momentum blend of the water and air incident fields with the incident
    volume fraction alphaI (momentum_swense_Rev1, Definition 1 and Section
    12.1), gamma = alphaI limited to [0, 1]:

        rhoI = gamma rhoW + (1 - gamma) rhoA
        mI   = gamma rhoW UW + (1 - gamma) rhoA UA
        UInc = mI/rhoI
        pInc = gamma pW + (1 - gamma) pA              [output only]

    Faces: a_f = linear face value of gamma, W_f = a_f rhoW/rhoI_f,
        UfInc = W_f UW + (1 - W_f) UA                 [raw, not projected]

    alphaI is advected by the solver. blendWeight() gives its initial value,
    the relaxation target and the farfield boundary value:
        clamp(0.5 - (z - eta)/delta, 0, 1)

    Water (alphaI = 1 region), Grid2Grid, in the wave frame:
        UW = U_grid2grid(x, t)
        pW = rho_w (pExcess(x, t) - g_mag z_w)        [hydrostatic included]
    Air (alphaI = 0 region), Air2Air:
        UA = R^T U_air(x_w, y_w, z - seaLevel)
        pA = p_air(x_w, y_w, z - seaLevel) - rho_a g z_w   [hydrostatic]

    Wave frame (Grid2Grid convention, from Grid2GridCoeffs):
        x_w = R(theta) (x - referencePosition),  R(theta) = [c -s; s c]
    The Air2Air time is the Grid2Grid time t + referenceTime.

    The water value is evaluated only where gamma > 0 (a_f > 0 on faces) and
    the air value only where gamma < 1 (a_f < 1).

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
#include "swenseHeightRampForcing.H"
#include "PtrList.H"

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

        //- Air velocity handed to the blend (airFill in waveProperties):
        //  air2air (default, unchanged), water, or layer (smoothstep between
        //  the water and Air2Air velocity over zeta = z - eta in
        //  [matchStart, matchEnd] m)
        const word airFill_;
        const scalar matchStart_;
        const scalar matchEnd_;

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
        //- One flat point per distinct (x, y) row, to read eta of the row
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

        //- [rampMode height] the height-ramped forcing (null otherwise).
        //  Air2Air (source theory) per level of the forcing: a level is
        //  built from its theory when it becomes active and freed when it is
        //  no longer one of the two neighbours; air values are interpolated
        //  with the weights of the water side (levels airI0_, airI1_, weights
        //  airW0_, airW1_ of the last evaluate(); index < 0: weight 0)
        const swenseHeightRampForcing* hr_;
        dictionary a2aDict_;
        PtrList<Air2AirWrapper> airLev_;
        label airI0_, airI1_;
        scalar airW0_, airW1_;

    // Private member functions

        //- [rampMode height] build the Air2Air of level i from its theory
        //  and its tables for the current points
        void makeAirLevel(const label i);

        //- Air2Air tables of one wrapper for the current points
        void airTables(Air2AirWrapper& a);

        //- Rotation from the mesh frame to the wave frame (about z)
        vector toWave(const vector& x) const;

        //- Rotation from the wave frame back to the mesh frame (about z)
        vector toMesh(const vector& u) const;

        //- Water velocity and pressure (mesh frame)
        vector waterU(const point& x, const scalar t) const;
        scalar waterP(const point& x, const scalar t) const;

        //- Air velocity (mesh frame) and pressure at flat point ip, from the
        //  Air2Air tables (one evaluation for both)
        vector airUP(const label ip, scalar& p) const;

        //- Air velocity handed to the blend at flat point ip, from the Air2Air
        //  velocity UAir of airUP (returned unchanged for airFill air2air)
        vector airVelocity
        (
            const label ip,
            const scalar t,
            const vector& UAir
        ) const;

        //- Water velocity and absolute water pressure at flat point ip (the
        //  per-step Grid2Grid values when they are current)
        vector waterUAt(const label ip, const scalar t) const;
        scalar waterPAt(const label ip, const scalar t) const;

        //- Log only: count a water value used above eta + d_g or an air
        //  value used below eta - d_g
        void countDrift
        (
            const label ip,
            const scalar t,
            const bool useWater,
            const bool useAir
        ) const;

        //- Momentum blend at a cell-value point ip for the water fraction
        //  gamma (spec step 1b): u = mI/rhoI. p is the output pressure
        //  gamma pW + (1 - gamma) pA. Water only where gamma > 0, air only
        //  where gamma < 1
        void blendCell
        (
            const label ip,
            const scalar t,
            const scalar gamma,
            vector& u,
            scalar& p
        ) const;

        //- Raw face velocity at a face point ip for the face fraction af
        //  (spec step 1c): W_f uW + (1 - W_f) uA
        vector blendFace
        (
            const label ip,
            const scalar t,
            const scalar af
        ) const;

        //- Grid2Grid eta at flat point ip: the per-step value when it is
        //  current for time t, otherwise evaluated here
        scalar etaAt(const label ip, const scalar t) const;

        //- Point tables of the mesh: flat order cells, boundary cells, faces,
        //  boundary faces (the order of the update loops), distinct (x, y) for
        //  eta and distinct x, y, z for the Air2Air tables
        void buildPointTables();

        //- True if a cell or face centre differs from the stored point list
        bool centresMoved() const;

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

        //- Evaluate the water and air wave models for time t (no blending).
        //  Call once per time step before blendAlpha() and blend()
        void evaluate(const scalar t);

        //- Momentum blend of the incident fields at the time of the last
        //  evaluate(), with the incident volume fraction alphaI (spec steps
        //  1b and 1c): UInc = mI/rhoI in the cells, the raw face velocity
        //  UfInc with the linear face value of alphaI, and pInc for output
        void blend(const volScalarField& alphaI);

        //- Exact wetted volume fraction below the ramped incident surface at
        //  time t: cells cut on their current points (swenseSurface, the
        //  routine of incidentResetInZone), boundary faces by the wetted face
        //  area. Cells whose points are all below (above) the surface are
        //  1 (0) without a cut.
        tmp<volScalarField> exactAlpha(const scalar t) const;

        //- Relaxation and initial target: exactAlpha(t) with waveProperties
        //  alphaHExact on, otherwise blendAlpha(t)
        tmp<volScalarField> targetAlpha(const scalar t) const;

        //- waveProperties alphaHExact (default off)
        static bool alphaHExact(const objectRegistry& db);

        //- Cell blend of the listed cells for the fractions gamma (clipped to
        //  [0, 1] here), at the time of the last evaluate(): u = mI/rhoI and
        //  p = gamma pW + (1 - gamma) pA, as blend() gives in those cells.
        //  Nothing is stored (incidentResetInZone).
        void blendCells
        (
            const labelUList& cells,
            const scalarField& gamma,
            vectorField& u,
            scalarField& p
        ) const;

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
        volScalarField& pIncRW() { return pInc_; }
        const volScalarField& pInc() const { return pInc_; }
        const surfaceVectorField& UfInc() const { return UfInc_; }

        //- Side pressure split (incidentHydrostatic discrete), at the time
        //  of the last evaluate(), cells and boundary faces. Side: water
        //  where z_w < eta, air above (eta of the local HOS surface):
        //    pHyd = rho_side g (eta - z_w)            (analytic hydrostatic)
        //    pDyn = p_side - pHyd                     (p_side = waterP or the
        //                                              Air2Air pressure)
        void sidePressure(volScalarField& pHyd, volScalarField& pDyn) const;

        //- Surface rows (x, y, 0) and Grid2Grid eta of row r at time t
        const pointField& etaPts() const { return etaPts_; }
        scalar etaRow(const label r, const scalar t) const { return etaAt(rowPt_[r], t); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
