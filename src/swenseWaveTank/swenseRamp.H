/*---------------------------------------------------------------------------*\
Description
    Start-up ramp of the SWENSE incident fields. One law, one key, used by
    swenseIncident (water eta, U, pExcess; Air2Air U and Pd), swenseSurface
    and the swenseWaveAlpha and swenseWaveVelocity boundary conditions.

    rampTsoft in waveProperties (default 0 = no ramp):
        f(t) = sin(pi/2 min(t, Tsoft)/Tsoft)       (waves2Foam Grid2Grid law)
    Values are multiplied by f(t). Time derivatives would be f' q + f q_t;
    the SWENSE solver consumes none of them.

    The Grid2Grid soft start (Grid2GridCoeffs Tsoft) and the theory soft
    start (waveTheoryForcingCoeffs Tsoft) must be 0, so the ramp is applied
    once.

SourceFiles
    swenseRamp.H (header only)

\*---------------------------------------------------------------------------*/

#ifndef swenseRamp_H
#define swenseRamp_H

#include "IOdictionary.H"
#include "objectRegistry.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace swenseRamp
{

//- Ramp time from the waveProperties of this database (0 = no ramp)
inline scalar Tsoft(const objectRegistry& db)
{
    return db.lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<scalar>("rampTsoft", 0.0);
}

//- Ramp factor at time t; exactly 1 when there is no ramp
inline scalar factor(const scalar Tsoft, const scalar t)
{
    if (Tsoft <= 0)
    {
        return 1.0;
    }

    return Foam::sin(constant::mathematical::pi/2*min(t, Tsoft)/Tsoft);
}

}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
