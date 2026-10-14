/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseIncident.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * Static Member Functions  * * * * * * * * * //

Foam::scalar Foam::swenseIncident::blendWeight
(
    const scalar zRel,
    const scalar delta
)
{
    return min(max(0.5 - zRel/delta, scalar(0)), scalar(1));
}


Foam::scalar Foam::swenseIncident::interfaceCellHeight(const fvMesh& mesh)
{
    List<scalar> zs;

    forAll(mesh.C(), celli)
    {
        if (mag(mesh.C()[celli].z()) < 0.1)
        {
            zs.append(mesh.C()[celli].z());
        }
    }

    Foam::sort(zs);

    scalar dz = 0.0;

    for (label i = 1; i < zs.size(); ++i)
    {
        if (zs[i] - zs[i-1] > SMALL)
        {
            dz = max(dz, zs[i] - zs[i-1]);
        }
    }

    return dz;
}


Foam::scalar Foam::swenseIncident::blendThickness(const objectRegistry& db)
{
    return db.lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<scalar>("blendThickness", 0.024);
}


// * * * * * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * //

Foam::swenseIncident::swenseIncident
(
    const fvMesh& mesh,
    const waveTheories::externalWaveForcing& ext,
    const dimensionedScalar& rhoWater,
    const dimensionedScalar& rhoAir
)
:
    mesh_(mesh),
    ext_(ext),
    air_(),
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
    delta_(blendThickness(mesh.thisDb())),
    guard_(0.0),
    rhoAir_(rhoAir.value()),
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
    )
{
    const dictionary& wp =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties");

    // Guard distance of the interface blend (document Rev4, rule 5)
    const scalar dzInterface = interfaceCellHeight(mesh);
    guard_ = max(delta_, 3.0*dzInterface);

    Info<< "SWENSE guard: delta = " << delta_ << " m, dz at interface = "
        << dzInterface << " m, guard d_g = " << guard_ << " m" << endl;

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

        // The Air2Air pressure is -rho_air (psi_t + ...): its density must be
        // the air density of the two-phase model
        const scalar rhoAir2Air = a2a.lookupOrDefault<scalar>("rhoAir", 1.225);

        if (mag(rhoAir2Air - rhoAir_) > SMALL*rhoAir_)
        {
            FatalErrorInFunction
                << "Air2Air rhoAir (" << rhoAir2Air << ") must equal the air "
                << "density rho2 (" << rhoAir_ << ")" << nl
                << "    set rhoAir in farfieldCoeffs/Air2Air of waveProperties"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
    // No Air2Air block: still air (the calm-water test of document 16.4)
    if (!air_.valid())
    {
        return vector::zero;
    }

    const vector xw = toWave(x - refPos_);

    return toMesh(air_->velocity(xw.x(), xw.y(), x.z() - seaLevel_));
}


Foam::scalar Foam::swenseIncident::airP(const point& x) const
{
    // Full air pressure (document eq. 5.2): Air2Air part plus hydrostatic
    // -rho_a g z_w, with the same height as the water hydrostatic term.
    // Without an Air2Air block the air is at rest: only the hydrostatic part
    const scalar zw = (x - refPos_).z();

    if (!air_.valid())
    {
        return -rhoAir_*magG_*zw;
    }

    const vector xw = toWave(x - refPos_);

    return air_->pressure(xw.x(), xw.y(), x.z() - seaLevel_)
         - rhoAir_*magG_*zw;
}


void Foam::swenseIncident::weights
(
    const point& x,
    const scalar t,
    const scalar alphaI,
    scalar& wWater,
    scalar& wAir
) const
{
    const scalar eps = 1e-6;
    const scalar W = min(max((alphaI - eps)/(1.0 - 2.0*eps), scalar(0)), scalar(1));

    const scalar zRel = x.z() - ext_.eta(x, t);

    const bool evalWater = (W > 0) && (zRel <= guard_);
    const bool evalAir = (W < 1) && (zRel >= -guard_);

    if (!evalWater && !evalAir)
    {
        FatalErrorInFunction
            << "alphaI = " << alphaI << " at " << x
            << " (z - eta = " << zRel << ") is outside both incident fields:"
            << " alphaI has drifted from eta" << exit(FatalError);
    }

    // Exactly 0 for a part that is not evaluated, exactly 1 for the other
    wWater = evalWater ? (evalAir ? W : scalar(1)) : scalar(0);
    wAir = evalAir ? (evalWater ? 1.0 - W : scalar(1)) : scalar(0);

    if (W > 0 && W < 1 && (zRel > guard_ || zRel < -guard_))
    {
        WarningInFunction
            << "W = " << W << " strictly between 0 and 1 at " << x
            << " with |z - eta| = " << mag(zRel) << " > d_g: alphaI has drifted"
            << endl;
    }
}


Foam::vector Foam::swenseIncident::blendU
(
    const point& x,
    const scalar t,
    const scalar alphaI
) const
{
    scalar wWater, wAir;
    weights(x, t, alphaI, wWater, wAir);

    vector u(Zero);

    if (wWater > 0)
    {
        u += wWater*waterU(x, t);
    }

    if (wAir > 0)
    {
        u += wAir*airU(x);
    }

    return u;
}


Foam::scalar Foam::swenseIncident::blendP
(
    const point& x,
    const scalar t,
    const scalar alphaI
) const
{
    scalar wWater, wAir;
    weights(x, t, alphaI, wWater, wAir);

    scalar p = 0.0;

    if (wWater > 0)
    {
        p += wWater*waterP(x, t);
    }

    if (wAir > 0)
    {
        p += wAir*airP(x);
    }

    return p;
}


Foam::volScalarField Foam::swenseIncident::blendAlpha(const scalar t) const
{
    volScalarField alpha
    (
        IOobject("alphaBlendInit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("alphaBlendInit", dimless, 0.0)
    );

    const volVectorField& C = mesh_.C();

    forAll(C, celli)
    {
        const point& x = C[celli];

        alpha[celli] = blendWeight(x.z() - ext_.eta(x, t), delta_);
    }

    // Boundary faces: the same target at the face centres
    const surfaceVectorField::Boundary& Cfb = mesh_.Cf().boundaryField();

    forAll(Cfb, patchi)
    {
        forAll(Cfb[patchi], facei)
        {
            const point& x = Cfb[patchi][facei];

            alpha.boundaryFieldRef()[patchi][facei] =
                blendWeight(x.z() - ext_.eta(x, t), delta_);
        }
    }

    return alpha;
}


void Foam::swenseIncident::update
(
    const scalar t,
    const volScalarField& alphaI
)
{
    // Grid2Grid time is t + referenceTime; Air2Air uses the same time
    const scalar tAir = t + refTime_;

    if (air_.valid())
    {
        air_->update(tAir);
    }

    // Cells: blend with the advected alphaI (clamped to [0, 1])
    const volVectorField& C = mesh_.C();

    forAll(C, celli)
    {
        const point& x = C[celli];
        UInc_[celli] = blendU(x, t, alphaI[celli]);
        pInc_[celli] = blendP(x, t, alphaI[celli]);
    }

    // Boundary patches of the cell fields: same blend, patch centres
    forAll(C.boundaryField(), patchi)
    {
        const vectorField& Cb = C.boundaryField()[patchi];
        const scalarField& ab = alphaI.boundaryField()[patchi];

        forAll(Cb, facei)
        {
            UInc_.boundaryFieldRef()[patchi][facei] = blendU(Cb[facei], t, ab[facei]);
            pInc_.boundaryFieldRef()[patchi][facei] = blendP(Cb[facei], t, ab[facei]);
        }
    }

    // Faces: alphaIf is the linear interpolation of alphaI
    const surfaceScalarField alphaIf(fvc::interpolate(alphaI));

    const vectorField& Cfi = mesh_.Cf().internalField();
    const scalarField& aif = alphaIf.internalField();

    forAll(Cfi, facei)
    {
        const point& x = Cfi[facei];
        UfInc_.primitiveFieldRef()[facei] = blendU(x, t, aif[facei]);
        pfInc_.primitiveFieldRef()[facei] = blendP(x, t, aif[facei]);
    }

    const surfaceVectorField::Boundary& Cfb = mesh_.Cf().boundaryField();

    forAll(Cfb, patchi)
    {
        const fvsPatchScalarField& ab = alphaIf.boundaryField()[patchi];

        forAll(Cfb[patchi], facei)
        {
            const point& x = Cfb[patchi][facei];
            UfInc_.boundaryFieldRef()[patchi][facei] = blendU(x, t, ab[facei]);
            pfInc_.boundaryFieldRef()[patchi][facei] = blendP(x, t, ab[facei]);
        }
    }
}

// ************************************************************************* //
