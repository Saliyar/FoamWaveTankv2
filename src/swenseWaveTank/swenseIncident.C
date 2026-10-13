/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseIncident.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * //

Foam::swenseIncident::swenseIncident
(
    const fvMesh& mesh,
    const waveTheories::externalWaveForcing& ext,
    const dimensionedScalar& rhoWater
)
:
    mesh_(mesh),
    ext_(ext),
    air_(),
    surf_(mesh, ext, 0.0),
    refPos_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .subDict("Grid2GridCoeffs")
            .lookupOrDefault<vector>("referencePosition", vector::zero)
    ),
    theta_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .subDict("Grid2GridCoeffs")
            .lookupOrDefault<scalar>("referenceTheta", 0.0)*M_PI/180.0
    ),
    refTime_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .subDict("Grid2GridCoeffs")
            .lookupOrDefault<scalar>("referenceTime", 0.0)
    ),
    seaLevel_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<scalar>("seaLevel", 0.0)
    ),
    rhoWater_(rhoWater.value()),
    magG_
    (
        mag(meshObjects::gravity::New(mesh.thisDb().time()).value())
    ),
    UInc_
    (
        IOobject("UInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("UInc", dimVelocity, vector::zero)
    ),
    pInc_
    (
        IOobject("pInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pInc", dimPressure, 0.0)
    ),
    alphaInc_
    (
        IOobject("alphaInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("alphaInc", dimless, 0.0)
    ),
    UfInc_
    (
        IOobject("UfInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("UfInc", dimVelocity, vector::zero)
    ),
    pfInc_
    (
        IOobject("pfInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pfInc", dimPressure, 0.0)
    ),
    cellWater_(mesh.nCells(), 0.0)
{
    const dictionary& wp =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties");

    // The SWENSE path needs the dynamic pressure without the soft start
    if (wp.subDict("Grid2GridCoeffs").lookupOrDefault<scalar>("Tsoft", 0.0) > SMALL)
    {
        FatalErrorInFunction
            << "swenseWaveFoam requires Grid2GridCoeffs Tsoft = 0"
            << exit(FatalError);
    }

    // Air2Air block (optional), in the farfield coefficients
    if (wp.subDict("farfieldCoeffs").found("Air2Air"))
    {
        const dictionary& a2a = wp.subDict("farfieldCoeffs").subDict("Air2Air");
        fileName modesFile;
        a2a.lookup("modesFile") >> modesFile;
        air_.reset(new Air2AirWrapper(modesFile, a2a));
    }
}


// * * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

Foam::vector Foam::swenseIncident::toWave(const vector& x) const
{
    // R(theta) = [c -s; s c], the Grid2Grid convention
    const scalar c = Foam::cos(theta_);
    const scalar s = Foam::sin(theta_);

    return vector(c*x.x() - s*x.y(), s*x.x() + c*x.y(), x.z());
}


Foam::vector Foam::swenseIncident::toMesh(const vector& u) const
{
    // R(theta)^T
    const scalar c = Foam::cos(theta_);
    const scalar s = Foam::sin(theta_);

    return vector(c*u.x() + s*u.y(), -s*u.x() + c*u.y(), u.z());
}


Foam::vector Foam::swenseIncident::waterU
(
    const point& x,
    const scalar t
) const
{
    return ext_.U(x, t);
}


Foam::scalar Foam::swenseIncident::waterP
(
    const point& x,
    const scalar t
) const
{
    // pExcess is the dynamic pressure per unit water density; the
    // hydrostatic term uses the wave-frame height (rotation keeps z)
    const scalar zw = (x - refPos_).z();

    return rhoWater_*(ext_.pExcess(x, t) - magG_*zw);
}


Foam::vector Foam::swenseIncident::airU(const point& x) const
{
    if (!air_.valid())
    {
        FatalErrorInFunction
            << "Air point found but no Air2Air block is configured"
            << exit(FatalError);
    }

    const vector xw = toWave(x - refPos_);

    return toMesh(air_->velocity(xw.x(), xw.y(), x.z() - seaLevel_));
}


Foam::scalar Foam::swenseIncident::airP(const point& x) const
{
    if (!air_.valid())
    {
        FatalErrorInFunction
            << "Air point found but no Air2Air block is configured"
            << exit(FatalError);
    }

    const vector xw = toWave(x - refPos_);

    return air_->pressure(xw.x(), xw.y(), x.z() - seaLevel_);
}


Foam::scalar Foam::swenseIncident::interfaceP
(
    const point& x,
    const scalar t
) const
{
    const point xs(x.x(), x.y(), ext_.eta(x, t));

    return waterP(xs, t);
}


void Foam::swenseIncident::update(const scalar t)
{
    // Grid2Grid time is t + referenceTime; Air2Air uses the same time
    const scalar tAir = t + refTime_;

    if (air_.valid())
    {
        air_->update(tAir);
    }

    surf_.setTime(t);

    const scalarField& V = mesh_.V();

    // Cells: classification by centre, values, wetted volume fraction
    const volVectorField& C = mesh_.C();

    forAll(C, celli)
    {
        const point& x = C[celli];
        const bool water = (x.z() - ext_.eta(x, t)) < 0.0;

        cellWater_[celli] = water ? 1.0 : 0.0;

        if (water)
        {
            UInc_[celli] = waterU(x, t);
            pInc_[celli] = waterP(x, t);
        }
        else
        {
            UInc_[celli] = airU(x);
            pInc_[celli] = airP(x);
        }

        alphaInc_[celli] = surf_.wetVolume(celli)/V[celli];
    }

    // Boundary patches of the cell fields: same classification, patch centres
    forAll(C.boundaryField(), patchi)
    {
        const vectorField& Cb = C.boundaryField()[patchi];

        forAll(Cb, facei)
        {
            const point& x = Cb[facei];
            const bool water = (x.z() - ext_.eta(x, t)) < 0.0;

            UInc_.boundaryFieldRef()[patchi][facei] =
                water ? waterU(x, t) : airU(x);
            pInc_.boundaryFieldRef()[patchi][facei] =
                water ? waterP(x, t) : airP(x);
            alphaInc_.boundaryFieldRef()[patchi][facei] = water ? 1.0 : 0.0;
        }
    }

    // Internal faces: classification by face centre, interface correction
    // where the two neighbouring cells differ
    const vectorField& Cfi = mesh_.Cf().internalField();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    forAll(Cfi, facei)
    {
        const point& x = Cfi[facei];
        const bool water = (x.z() - ext_.eta(x, t)) < 0.0;

        UfInc_.primitiveFieldRef()[facei] = water ? waterU(x, t) : airU(x);
        pfInc_.primitiveFieldRef()[facei] = water ? waterP(x, t) : airP(x);

        if (cellWater_[own[facei]] != cellWater_[nei[facei]])
        {
            pfInc_.primitiveFieldRef()[facei] = interfaceP(x, t);
        }
    }

    // Boundary faces: classification by face centre
    const surfaceVectorField::Boundary& Cfb = mesh_.Cf().boundaryField();

    forAll(Cfb, patchi)
    {
        forAll(Cfb[patchi], facei)
        {
            const point& x = Cfb[patchi][facei];
            const bool water = (x.z() - ext_.eta(x, t)) < 0.0;

            UfInc_.boundaryFieldRef()[patchi][facei] =
                water ? waterU(x, t) : airU(x);
            pfInc_.boundaryFieldRef()[patchi][facei] =
                water ? waterP(x, t) : airP(x);
        }
    }
}

// ************************************************************************* //
