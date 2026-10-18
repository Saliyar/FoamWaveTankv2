/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseIncident.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"
#include "Grid2Grid.H"

#include <cmath>
#include <map>
#include <utility>
#include <vector>

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
    ),
    psiWf_
    (
        IOobject("psiWf", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("psiWf", dimVelocity*dimVelocity, 0.0)
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

    g2g_ = dynamic_cast<const waveTheories::Grid2Grid*>(&ext_);
    fastEval_ = false;

    buildPointTables();

    if (g2g_)
    {
        // [option A] stored spline basis of all points, once
        g2g_->setupPoints(pts_);
        etaPt_.setSize(pts_.size(), 0.0);
        pExPt_.setSize(pts_.size(), 0.0);
        UPt_.setSize(pts_.size(), vector::zero);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swenseIncident::buildPointTables()
{
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();

    // Flat point list in the order of the update loops
    std::vector<point> pts;

    forAll(C, celli)
    {
        pts.push_back(C[celli]);
    }

    offCb_ = label(pts.size());

    forAll(C.boundaryField(), patchi)
    {
        forAll(C.boundaryField()[patchi], facei)
        {
            pts.push_back(C.boundaryField()[patchi][facei]);
        }
    }

    offCf_ = label(pts.size());

    forAll(Cf, facei)
    {
        pts.push_back(Cf[facei]);
    }

    offCfb_ = label(pts.size());

    forAll(Cf.boundaryField(), patchi)
    {
        forAll(Cf.boundaryField()[patchi], facei)
        {
            pts.push_back(Cf.boundaryField()[patchi][facei]);
        }
    }

    const label n = label(pts.size());

    pts_.setSize(n);
    forAll(pts_, ip)
    {
        pts_[ip] = pts[ip];
    }

    // Grid2Grid eta depends on (x, y) only: one row per distinct (x, y)
    std::map<std::pair<long long, long long>, label> rowOf;
    std::vector<point> etaRows;

    etaIdx_.setSize(n);
    forAll(pts_, ip)
    {
        const std::pair<long long, long long> key
        (
            std::llround(pts_[ip].x()*1e14), std::llround(pts_[ip].y()*1e14)
        );

        auto it = rowOf.find(key);

        if (it == rowOf.end())
        {
            it = rowOf.insert(std::make_pair(key, label(etaRows.size()))).first;
            etaRows.push_back(point(pts_[ip].x(), pts_[ip].y(), 0.0));
        }

        etaIdx_[ip] = it->second;
    }

    etaPts_.setSize(etaRows.size());
    forAll(etaPts_, r)
    {
        etaPts_[r] = etaRows[r];
    }

    etaVal_.setSize(etaRows.size(), 0.0);
    etaTime_ = -GREAT;

    Info<< "SWENSE point tables: " << n << " points, " << etaRows.size()
        << " distinct (x, y) for eta" << endl;

    if (air_.valid())
    {
        scalarList xs(n), ys(n), zs(n);

        forAll(pts_, ip)
        {
            const vector xw = toWave(pts_[ip] - refPos_);
            xs[ip] = xw.x();
            ys[ip] = xw.y();
            zs[ip] = pts_[ip].z() - seaLevel_;
        }

        air_->solver().buildTables(xs, ys, zs, airIx_, airIy_, airIz_);
    }
}


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


Foam::vector Foam::swenseIncident::airUP
(
    const label ip,
    scalar& p
) const
{
    // Full air field (document eq. 5.2): Air2Air part plus hydrostatic
    // -rho_a g z_w, from the Air2Air tables. Without an Air2Air block the
    // air is at rest: only the hydrostatic part
    const scalar zw = (pts_[ip] - refPos_).z();

    if (!air_.valid())
    {
        p = -rhoAir_*magG_*zw;
        return vector::zero;
    }

    vector Uw;
    scalar Pd;
    air_->solver().flowTables(airIx_[ip], airIy_[ip], airIz_[ip], Uw, Pd);

    p = Pd - rhoAir_*magG_*zw;

    return toMesh(Uw);
}


Foam::scalar Foam::swenseIncident::etaAt(const label ip, const scalar t) const
{
    if (etaTime_ == t)
    {
        return fastEval_ ? etaPt_[ip] : etaVal_[etaIdx_[ip]];
    }

    return ext_.eta(pts_[ip], t);
}


void Foam::swenseIncident::weights
(
    const point& x,
    const scalar eta,
    const scalar alphaI,
    scalar& wWater,
    scalar& wAir
) const
{
    const scalar eps = 1e-6;
    const scalar W = min(max((alphaI - eps)/(1.0 - 2.0*eps), scalar(0)), scalar(1));

    const scalar zRel = x.z() - eta;

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


void Foam::swenseIncident::blendUP
(
    const label ip,
    const scalar t,
    const scalar alphaI,
    vector& u,
    scalar& p
) const
{
    const point& x = pts_[ip];

    // Weights once per point, used by the velocity only (scratch fix)
    const scalar eta = etaAt(ip, t);
    scalar wWater, wAir;
    weights(x, eta, alphaI, wWater, wAir);

    u = Zero;
    p = 0.0;

    scalar pa = 0.0;
    bool havePa = false;

    if (wWater > 0)
    {
        if (fastEval_)
        {
            u += wWater*UPt_[ip];
        }
        else
        {
            u += wWater*waterU(x, t);
        }
    }

    if (wAir > 0)
    {
        const vector Ua = airUP(ip, pa);
        havePa = true;
        u += wAir*Ua;
    }

    // Pressure (scratch fix): one side only, relative to its value at the
    // surface of this (x, y) row; only the side that is used is evaluated
    const scalar zRel = x.z() - eta;
    const label r = etaIdx_[ip];

    if (zRel <= 0)
    {
        const scalar pW =
            fastEval_
          ? rhoWater_*(pExPt_[ip] - magG_*(x - refPos_).z())
          : waterP(x, t);
        p = pW - pWs_[r];
    }
    else
    {
        if (!havePa)
        {
            airUP(ip, pa);
        }
        p = pa - pAs_[r];
    }
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
        alpha[celli] = blendWeight(C[celli].z() - etaAt(celli, t), delta_);
    }

    // Boundary faces: the same target at the face centres
    const surfaceVectorField::Boundary& Cfb = mesh_.Cf().boundaryField();

    label ip = offCfb_;

    forAll(Cfb, patchi)
    {
        forAll(Cfb[patchi], facei)
        {
            alpha.boundaryFieldRef()[patchi][facei] =
                blendWeight(Cfb[patchi][facei].z() - etaAt(ip, t), delta_);
            ++ip;
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

    if (g2g_)
    {
        // [option A] Grid2Grid eta, U and pExcess at all stored points, once per step
        g2g_->evaluatePoints(t, etaPt_, UPt_, pExPt_);
        fastEval_ = true;
    }
    else
    {
        // Grid2Grid eta once per distinct (x, y) for this time (used by all points)
        fastEval_ = false;
        forAll(etaPts_, r)
        {
            etaVal_[r] = ext_.eta(etaPts_[r], t);
        }
    }
    etaTime_ = t;

    // Scratch fix: surface pressures once per step for every distinct (x, y) row
    if (rowPt_.size() != etaPts_.size())
    {
        rowPt_.setSize(etaPts_.size());
        rowPt_ = -1;
        forAll(etaIdx_, ip)
        {
            if (rowPt_[etaIdx_[ip]] < 0)
            {
                rowPt_[etaIdx_[ip]] = ip;
            }
        }
    }
    pWs_.setSize(etaPts_.size());
    pAs_.setSize(etaPts_.size());
    forAll(etaPts_, r)
    {
        const point xs(etaPts_[r].x(), etaPts_[r].y(), etaAt(rowPt_[r], t));
        const scalar zws = (xs - refPos_).z();

        pWs_[r] = waterP(xs, t);

        if (!air_.valid())
        {
            pAs_[r] = -rhoAir_*magG_*zws;
        }
        else
        {
            const vector xw = toWave(xs - refPos_);
            vector Uw;
            scalar Pd;
            air_->solver().flow(xw.x(), xw.y(), xs.z() - seaLevel_, Uw, Pd);
            pAs_[r] = Pd - rhoAir_*magG_*zws;
        }
    }

    // Cells: blend with the advected alphaI (clamped to [0, 1])
    const volVectorField& C = mesh_.C();

    forAll(C, celli)
    {
        blendUP(celli, t, alphaI[celli], UInc_[celli], pInc_[celli]);
    }

    // Boundary patches of the cell fields: same blend, patch centres
    label ip = offCb_;

    forAll(C.boundaryField(), patchi)
    {
        const scalarField& ab = alphaI.boundaryField()[patchi];

        forAll(ab, facei)
        {
            blendUP
            (
                ip, t, ab[facei],
                UInc_.boundaryFieldRef()[patchi][facei],
                pInc_.boundaryFieldRef()[patchi][facei]
            );
            ++ip;
        }
    }

    // Faces: alphaIf is the linear interpolation of alphaI
    const surfaceScalarField alphaIf(fvc::interpolate(alphaI));

    const scalarField& aif = alphaIf.internalField();

    forAll(aif, facei)
    {
        blendUP
        (
            offCf_ + facei, t, aif[facei],
            UfInc_.primitiveFieldRef()[facei],
            pfInc_.primitiveFieldRef()[facei]
        );
        psiWf_.primitiveFieldRef()[facei] = psiWValue(offCf_ + facei, t);
    }

    const surfaceVectorField::Boundary& Cfb = mesh_.Cf().boundaryField();

    ip = offCfb_;

    forAll(Cfb, patchi)
    {
        const fvsPatchScalarField& ab = alphaIf.boundaryField()[patchi];

        forAll(Cfb[patchi], facei)
        {
            blendUP
            (
                ip, t, ab[facei],
                UfInc_.boundaryFieldRef()[patchi][facei],
                pfInc_.boundaryFieldRef()[patchi][facei]
            );
            psiWf_.boundaryFieldRef()[patchi][facei] = psiWValue(ip, t);
            ++ip;
        }
    }

}

// ************************************************************************* //


// Scratch: psiWf = (pW(face) - pWs[row]) / rho1, water pressure relative to the row surface
Foam::scalar Foam::swenseIncident::psiWValue(const label ip, const scalar t) const
{
    const scalar pw = fastEval_
      ? rhoWater_*(pExPt_[ip] - magG_*(pts_[ip] - refPos_).z())
      : waterP(pts_[ip], t);

    return (pw - pWs_[etaIdx_[ip]])/rhoWater_;
}
