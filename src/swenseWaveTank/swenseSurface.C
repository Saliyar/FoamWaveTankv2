/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseSurface.H"

// * * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::swenseSurface::swenseSurface
(
    const fvMesh& mesh,
    const waveTheories::externalWaveForcing& ext,
    const scalar t
)
:
    convexPolyhedral(mesh, true),
    ext_(ext),
    t_(t)
{}


// * * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

void Foam::swenseSurface::setTime(const scalar t)
{
    t_ = t;
}


void Foam::swenseSurface::signedPointToSurfaceDistance
(
    const pointField& pp,
    scalarField& sd
)
{
    forAll(pp, pointi)
    {
        sd[pointi] = signedPointToSurfaceDistance(pp[pointi]);
    }
}


Foam::scalar Foam::swenseSurface::signedPointToSurfaceDistance
(
    const point& pp
) const
{
    return pp.z() - ext_.eta(pp, t_);
}


Foam::scalar Foam::swenseSurface::wetVolume(const label celli)
{
    return dividePolyhedral(celli, point::zero, point::one).magNeg();
}


Foam::scalar Foam::swenseSurface::wetArea(const label facei)
{
    return divideFace(facei).negMag();
}


Foam::vector Foam::swenseSurface::wetCentre(const label facei)
{
    return divideFace(facei).negCentre();
}

// ************************************************************************* //
