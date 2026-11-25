/*---------------------------------------------------------------------------*\
Class
    Foam::swenseHeightRampForcing

Description
    Height-ramped analytic incident wave (waveProperties rampMode height).

    At start-up one wave theory is built per height H_i = H i/n, i = 1..n
    (n = 20): the theory dictionary waveTheoryForcingCoeffs is copied with
    height H_i, setWaveProperties New(...)->set() writes the derived
    coefficients into an OStringStream, the text is parsed back and the
    theory is built from it. setWaveProperties runs on the master only, one
    height at a time (streamFunction runs fenton4Foam in the case directory),
    and the text is broadcast.

    At time t the target height is Hr = H f(t), f the swenseRamp law with
    rampTsoft. eta, U and pExcess are evaluated from the two theories whose
    heights bracket Hr and interpolated linearly in H. Below H_1 the lower
    neighbour is still water (all values 0). For t >= rampTsoft only H_n is
    evaluated.

    This is the one place that supplies the ramped incident. With rampMode
    height swenseRamp::Tsoft() is 0, so swenseIncident, swenseSurface and the
    swenseWaveAlpha and swenseWaveVelocity conditions read these values
    without a factor. The Air2Air air side (source theory) is built per
    level from level(i) and interpolated with the same weights
    (swenseIncident).

    The base class (waveTheoryForcing) keeps the full-height theory of the
    file coefficients; phaseSpeed() and waveNumber() are taken from it.

SourceFiles
    swenseHeightRampForcing.C

\*---------------------------------------------------------------------------*/

#ifndef swenseHeightRampForcing_H
#define swenseHeightRampForcing_H

#include "waveTheoryForcing.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class swenseHeightRampForcing Declaration
\*---------------------------------------------------------------------------*/

class swenseHeightRampForcing
:
    public waveTheories::waveTheoryForcing
{
    // Private data

        //- Full wave height [m] (waveTheoryForcingCoeffs height)
        scalar H_;

        //- Ramp time (rampTsoft) [s]
        scalar Tramp_;

        //- Heights H_i and their theories, i = 1..n stored at 0..n-1
        scalarList Hlev_;
        PtrList<waveTheories::waveTheory> lev_;


    // Private Member Functions

        swenseHeightRampForcing(const swenseHeightRampForcing&) = delete;
        void operator=(const swenseHeightRampForcing&) = delete;


public:

    //- Runtime type information
    TypeName("swenseHeightRampForcing");


    // Constructors

        //- Construct as the registered "externalWaveForcing"
        swenseHeightRampForcing
        (
            IOobject io,
            Time& rT,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~swenseHeightRampForcing() = default;


    // Member Functions

        //- Number of levels n
        label nLevels() const { return lev_.size(); }

        //- Height of level i (0..n-1) [m]
        scalar levelHeight(const label i) const { return Hlev_[i]; }

        //- Theory of level i (0..n-1), not ramped
        const waveTheories::waveTheory& level(const label i) const
        {
            return lev_[i];
        }

        //- Full height [m]
        scalar fullHeight() const { return H_; }

        //- Ramped height Hr = H f(t) [m]
        scalar rampedHeight(const scalar t) const;

        //- Interpolation at time t: value = w0 v(i0) + w1 v(i1); a level
        //  index < 0 has weight 0 (still water or not needed)
        void weights
        (
            const scalar t,
            label& i0,
            scalar& w0,
            label& i1,
            scalar& w1
        ) const;

        //- Ramped free-surface elevation [m]
        virtual scalar eta(const point& x, const scalar& time) const;

        //- Ramped excess pressure per unit water density [m2/s2]
        virtual scalar pExcess(const point& x, const scalar& time) const;

        //- Ramped water velocity [m/s]
        virtual vector U(const point& x, const scalar& time) const;


    // Static set-up functions (called by the solver; kept in this library so
    // that the solver needs no type information of the forcing classes)

        //- The registered external wave forcing: rampMode factor (default)
        //  the forcing of externalForcing (as createExternalWaveForcing.H),
        //  rampMode height a swenseHeightRampForcing (externalForcing must
        //  be waveTheoryForcing)
        static autoPtr<waveTheories::externalWaveForcing> newForcing
        (
            Time& rT,
            const fvMesh& mesh
        );

        //- Start-up consistency of an analytic incident wave at the full
        //  height (swenseIncidentResidual.H): stop above incidentResidualMax
        //  of waveProperties (default 0.005). Grid2Grid is not checked
        static void residualCheck
        (
            const waveTheories::externalWaveForcing& ext,
            const fvMesh& mesh
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
