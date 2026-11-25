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

    rampMode in waveProperties (default factor):
        factor  the values are multiplied by f(t) at each site (above)
        height  the wave height is ramped instead, Hr = H f(t), by the
                external forcing swenseHeightRampForcing, which is the one
                place that supplies the ramped incident. Tsoft() then returns
                0, so no site multiplies by f; rampTime() gives the ramp time.

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

//- rampMode of the waveProperties of this database: true for height
inline bool heightMode(const objectRegistry& db)
{
    const word mode
    (
        db.lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<word>("rampMode", word("factor"))
    );

    if (mode != "factor" && mode != "height")
    {
        FatalErrorInFunction
            << "rampMode must be factor or height, not " << mode
            << exit(FatalError);
    }

    return mode == "height";
}

//- Ramp time rampTsoft of the waveProperties (0 = no ramp), in either mode
inline scalar rampTime(const objectRegistry& db)
{
    return db.lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<scalar>("rampTsoft", 0.0);
}

//- Ramp time of the factor applied at the sites: rampTsoft with rampMode
//  factor, 0 (no factor anywhere) with rampMode height
inline scalar Tsoft(const objectRegistry& db)
{
    return heightMode(db) ? 0.0 : rampTime(db);
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
