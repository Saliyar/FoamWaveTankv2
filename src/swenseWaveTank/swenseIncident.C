/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseIncident.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"
#include "Grid2Grid.H"
#include "waveTheoryForcing.H"

#include <cmath>
#include <map>
#include <utility>
#include <vector>

// Drift of alphaI from eta (weights()): counted there, reported once per call
// of update() as one log line, instead of one warning per cell
static Foam::label driftCount = 0;
static Foam::scalar driftWorstZ = 0.0;
static Foam::point driftWorstX = Foam::point::zero;

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

    // Heights from every processor: a processor holding part of the mesh
    // would see a smaller gap, and the guard distance must be the same on all
    List<scalarList> zsProc(Pstream::nProcs());
    zsProc[Pstream::myProcNo()] = scalarList(zs);
    Pstream::allGatherList(zsProc);

    DynamicList<scalar> zsAll;
    forAll(zsProc, proc)
    {
        forAll(zsProc[proc], i)
        {
            zsAll.append(zsProc[proc][i]);
        }
    }

    Foam::sort(zsAll);

    scalar dz = 0.0;

    for (label i = 1; i < zsAll.size(); ++i)
    {
        if (zsAll[i] - zsAll[i-1] > SMALL)
        {
            dz = max(dz, zsAll[i] - zsAll[i-1]);
        }
    }

    return dz;
}


Foam::scalar Foam::swenseIncident::blendThickness(const objectRegistry& db)
{
    return db.lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<scalar>("blendThickness", 0.024);
}


// Grid2GridCoeffs of waveProperties for a Grid2Grid forcing. Any other forcing
// (for example waveTheoryForcing) is already in CFD coordinates, so an empty
// dictionary gives the identity wave frame and zero time offset.
static Foam::dictionary grid2GridCoeffsFor
(
    const Foam::objectRegistry& db,
    const Foam::waveTheories::externalWaveForcing& ext
)
{
    if (dynamic_cast<const Foam::waveTheories::Grid2Grid*>(&ext))
    {
        return db.lookupObject<IOdictionary>("waveProperties")
            .subDict("Grid2GridCoeffs");
    }

    return Foam::dictionary();
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
        grid2GridCoeffsFor(mesh.thisDb(), ext)
            .lookupOrDefault<vector>("referencePosition", vector::zero)
    ),
    theta_
    (
        grid2GridCoeffsFor(mesh.thisDb(), ext)
            .lookupOrDefault<scalar>("referenceTheta", 0.0)*M_PI/180.0
    ),
    refTime_
    (
        grid2GridCoeffsFor(mesh.thisDb(), ext)
            .lookupOrDefault<scalar>("referenceTime", 0.0)
    ),
    seaLevel_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<scalar>("seaLevel", 0.0)
    ),
    rhoWater_(rhoWater.value()),
    Tsoft_(swenseRamp::Tsoft(mesh.thisDb())),
    tNow_(0.0),
    delta_(blendThickness(mesh.thisDb())),
    momentumBlend_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<bool>("momentumBlend", false)
    ),
    momentumBlendPressure_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<bool>("momentumBlendPressure", false)
    ),
    guard_(0.0),
    rhoAir_(rhoAir.value()),
    magG_
    (
        mag(meshObjects::gravity::New(mesh.thisDb().time()).value())
    ),
    UInc_
    (
        IOobject("UInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("UInc", dimVelocity, vector::zero)
    ),
    pInc_
    (
        IOobject("pInc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE),
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
    ),
    dUc_
    (
        IOobject("dUc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("dUc", dimVelocity, vector::zero)
    ),
    dPf_
    (
        IOobject("dPf", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("dPf", dimPressure, 0.0)
    ),
    pWc_
    (
        IOobject("pWc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pWc", dimPressure, 0.0)
    ),
    pAc_
    (
        IOobject("pAc", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pAc", dimPressure, 0.0)
    )
{
    const dictionary& wp =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties");

    // Guard distance of the interface blend (document Rev4, rule 5)
    const scalar dzInterface = interfaceCellHeight(mesh);
    guard_ = max(delta_, 3.0*dzInterface);

    Info<< "SWENSE guard: delta = " << delta_ << " m, dz at interface = "
        << dzInterface << " m, guard d_g = " << guard_ << " m" << endl;

    // The Grid2Grid soft start would be applied twice with the SWENSE ramp
    // (rampTsoft in waveProperties, applied once in swenseIncident). The
    // analytic theory checks its own Tsoft in waveTheoryForcing
    if (grid2GridCoeffsFor(mesh.thisDb(), ext_).lookupOrDefault<scalar>("Tsoft", 0.0) > SMALL)
    {
        FatalErrorInFunction
            << "swenseWaveFoam requires Grid2GridCoeffs Tsoft = 0; use rampTsoft"
            << " in waveProperties for the start-up ramp" << exit(FatalError);
    }

    // The Grid2Grid pExcess is used with g z measured from the origin
    if (mag(seaLevel_) > SMALL)
    {
        FatalErrorInFunction
            << "swenseWaveFoam requires seaLevel 0 in waveProperties (found "
            << seaLevel_ << " m)" << exit(FatalError);
    }

    // Air2Air block (optional), in the farfield coefficients
    if (wp.subDict("farfieldCoeffs").found("Air2Air"))
    {
        const dictionary& a2a = wp.subDict("farfieldCoeffs").subDict("Air2Air");
        const word source = a2a.lookupOrDefault<word>("source", word("modesFile"));

        if (source == "theory")
        {
            // Periodic modes from the wave theory (one period, t = 0)
            const waveTheories::waveTheoryForcing* wtf =
                dynamic_cast<const waveTheories::waveTheoryForcing*>(&ext_);

            if (!wtf)
            {
                FatalErrorInFunction
                    << "Air2Air source theory needs externalWaveForcing"
                    << " waveTheoryForcing" << exit(FatalError);
            }

            const scalar c = wtf->phaseSpeed();
            const scalar K = wtf->waveNumber();

            if (c <= 0 || K <= 0)
            {
                FatalErrorInFunction
                    << "Air2Air source theory: the wave theory gives no"
                    << " phase speed or wave number" << exit(FatalError);
            }

            // Period length L = N_L wavelength, N_L = 1
            const scalar Lp = 2*M_PI/K;
            const label M = a2a.lookupOrDefault<label>("nModes", 64);
            const label N = 4*M;

            std::vector<scalar> eta(N);
            for (label i = 0; i < N; ++i)
            {
                eta[i] = wtf->eta(point((i + 0.5)*Lp/N, 0, 0), 0.0);
            }

            // Length scale h = 1/K: the fundamental has k h = 1
            air_.reset(new Air2AirWrapper("", a2a));
            air_->setTheory
            (
                1.0/K, magG_, c, Lp, eta, M,
                a2a.lookupOrDefault<word>("cutoff", word("auto")),
                a2a.lookupOrDefault<scalar>("cutoffC", 0.5)
            );
        }
        else
        {
            fileName modesFile;
            a2a.lookup("modesFile") >> modesFile;
            air_.reset(new Air2AirWrapper(modesFile, a2a));
        }

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
    if (Tsoft_ > 0)
    {
        return ramp(t)*ext_.U(x, t);
    }

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

    if (Tsoft_ > 0)
    {
        return rhoWater_*(ramp(t)*ext_.pExcess(x, t) - magG_*zw);
    }

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

    // Air2Air part ramped; the hydrostatic part is not
    if (Tsoft_ > 0)
    {
        const scalar f = ramp(tNow_);
        Pd *= f;
        Uw *= f;
    }

    p = Pd - rhoAir_*magG_*zw;

    return toMesh(Uw);
}


Foam::scalar Foam::swenseIncident::ramp(const scalar t) const
{
    return swenseRamp::factor(Tsoft_, t);
}


Foam::scalar Foam::swenseIncident::etaAt(const label ip, const scalar t) const
{
    if (etaTime_ == t)
    {
        return fastEval_ ? etaPt_[ip] : etaVal_[etaIdx_[ip]];
    }

    if (Tsoft_ > 0)
    {
        return ramp(t)*ext_.eta(pts_[ip], t);
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
    const scalar a = min(max((alphaI - eps)/(1.0 - 2.0*eps), scalar(0)), scalar(1));

    // Momentum blend (waveProperties momentumBlend true): the water weight is
    // the share of the incident momentum density, W = a rhoW/rhoI, with rhoI
    // blended by the same a. Otherwise the water weight is a itself.
    const scalar rhoI = a*rhoWater_ + (1.0 - a)*rhoAir_;
    const scalar W = momentumBlend_ ? a*rhoWater_/rhoI : a;

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
        ++driftCount;

        if (mag(zRel) > driftWorstZ)
        {
            driftWorstZ = mag(zRel);
            driftWorstX = x;
        }
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

    // Momentum blend (waveProperties momentumBlend true): absolute pressures
    // blended linearly with alphaI, pI = alphaI pW + (1 - alphaI) pA, the form
    // used in the incident residual derivation. Only the sides that are
    // evaluated (weights above) enter; a side with zero weight is not used.
    // E1: the absolute blend is the only pressure convention (no momentumBlendPressure gate).
    {
        const scalar a =
            min(max((alphaI - 1e-6)/(1.0 - 2e-6), scalar(0)), scalar(1));

        const bool useW = wWater > 0;
        const bool useA = wAir > 0;

        scalar pWabs = 0.0, pAabs = 0.0;

        if (useW)
        {
            pWabs =
                fastEval_
              ? rhoWater_*(pExPt_[ip] - magG_*(x - refPos_).z())
              : waterP(x, t);
        }

        if (useA)
        {
            if (!havePa)
            {
                airUP(ip, pa);
            }
            pAabs = pa;
        }

        if (useW && useA)
        {
            p = a*pWabs + (1.0 - a)*pAabs;
        }
        else if (useW)
        {
            p = pWabs;
        }
        else
        {
            p = pAabs;
        }

        return;
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


void Foam::swenseIncident::movePoints()
{
    // The point lists follow the current cell and face centres. The
    // Grid2Grid basis is stored for these points, so it is set up again only
    // when a point has moved (a zero motion keeps the stored basis, as static)
    const pointField oldPts(pts_);

    buildPointTables();

    if (pts_ == oldPts)
    {
        return;
    }

    rowPt_.clear();

    if (g2g_)
    {
        g2g_->setupPoints(pts_);
        etaPt_.setSize(pts_.size(), 0.0);
        pExPt_.setSize(pts_.size(), 0.0);
        UPt_.setSize(pts_.size(), vector::zero);
    }
}


void Foam::swenseIncident::fieldsAt
(
    const label ip,
    const scalar t,
    vector& dU,
    scalar& pW,
    scalar& pa
) const
{
    const point& x = pts_[ip];

    const vector Ua = airUP(ip, pa);

    const vector Uw = fastEval_ ? UPt_[ip] : waterU(x, t);
    dU = Uw - Ua;

    // Same pW as blendUP(): the fastEval_ path when it is set
    pW =
        fastEval_
      ? rhoWater_*(pExPt_[ip] - magG_*(x - refPos_).z())
      : waterP(x, t);
}


void Foam::swenseIncident::dUdP
(
    const label ip,
    const scalar t,
    vector& dU,
    scalar& dP
) const
{
    scalar pW = 0.0;
    scalar pa = 0.0;
    fieldsAt(ip, t, dU, pW, pa);
    dP = pW - pa;
}


void Foam::swenseIncident::update
(
    const scalar t,
    const volScalarField& alphaI
)
{
    // One line for the drift events counted since the previous call
    Info<< "SWENSE drift: " << driftCount
        << " alphaI-drift events since the last update, worst |z - eta| = "
        << driftWorstZ << " m at " << driftWorstX << endl;
    driftCount = 0;
    driftWorstZ = 0.0;
    driftWorstX = Foam::point::zero;

    // Grid2Grid time is t + referenceTime; Air2Air uses the same time
    const scalar tAir = t + refTime_;
    tNow_ = t;

    if (air_.valid())
    {
        air_->update(tAir);
    }

    if (g2g_)
    {
        // [option A] Grid2Grid eta, U and pExcess at all stored points, once per step
        g2g_->evaluatePoints(t, etaPt_, UPt_, pExPt_);

        // Start-up ramp on the Grid2Grid values (eta, U, pExcess)
        if (Tsoft_ > 0)
        {
            const scalar f = ramp(t);
            etaPt_ *= f;
            UPt_ *= f;
            pExPt_ *= f;
        }
        fastEval_ = true;
    }
    else
    {
        // Grid2Grid eta once per distinct (x, y) for this time (used by all points)
        fastEval_ = false;
        forAll(etaPts_, r)
        {
            etaVal_[r] = ext_.eta(etaPts_[r], t);

            if (Tsoft_ > 0)
            {
                etaVal_[r] *= ramp(t);
            }
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
            if (Tsoft_ > 0)
            {
                Pd *= ramp(t);
            }

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

    // Momentum blend: dU = uW - uA at the cell centres and dP = pW - pA at
    // the faces, at the same points as blendUP (not filled when the switch is off)
    // Cell-centre dU, pW and pa, always (the velocity-blend F^v needs them too)
    {
        vector dU;
        scalar pW, pa;

        forAll(dUc_, celli)
        {
            fieldsAt(celli, t, dU, pW, pa);
            dUc_[celli] = dU;
            pWc_[celli] = pW;
            pAc_[celli] = pa;
        }

        label ipc = offCb_;

        forAll(dUc_.boundaryField(), patchi)
        {
            forAll(dUc_.boundaryField()[patchi], facei)
            {
                fieldsAt(ipc, t, dU, pW, pa);
                dUc_.boundaryFieldRef()[patchi][facei] = dU;
                pWc_.boundaryFieldRef()[patchi][facei] = pW;
                pAc_.boundaryFieldRef()[patchi][facei] = pa;
                ++ipc;
            }
        }
    }

    // dP = pW - pA at the faces (same points as blendUP), always: used by the
    // momentum blend and by the velocity-blend F^v
    {
        vector dU;
        scalar dP;

        forAll(dPf_.primitiveField(), facei)
        {
            dUdP(offCf_ + facei, t, dU, dP);
            dPf_.primitiveFieldRef()[facei] = dP;
        }

        label ipf = offCfb_;

        forAll(dPf_.boundaryField(), patchi)
        {
            forAll(dPf_.boundaryField()[patchi], facei)
            {
                dUdP(ipf, t, dU, dP);
                dPf_.boundaryFieldRef()[patchi][facei] = dP;
                ++ipf;
            }
        }
    }

}

// ************************************************************************* //


// psiWf = pW(face)/rho_w (absolute, diagnostic only since the absolute-pressure change)
Foam::scalar Foam::swenseIncident::psiWValue(const label ip, const scalar t) const
{
    const scalar pw = fastEval_
      ? rhoWater_*(pExPt_[ip] - magG_*(pts_[ip] - refPos_).z())
      : waterP(pts_[ip], t);

    // E2: absolute water psi, no surface subtraction. pWs_ and pAs_ are kept
    // only for diagnostics; the force uses pfInc/rhoI (swenseUcEqn.H), not this value.
    return pw/rhoWater_;
}
