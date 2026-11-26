/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseIncident.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"
#include "Grid2Grid.H"
#include "waveTheoryForcing.H"
#include "swenseSurface.H"

#include <cmath>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

// Wave-model values used outside the guard distance (countDrift()): counted
// there, reported by blend() as one log line
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
    airFill_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<word>("airFill", word("air2air"))
    ),
    matchStart_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<scalar>("matchStart", 0.0)
    ),
    matchEnd_
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<scalar>("matchEnd", 1.0)
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
    hr_(dynamic_cast<const swenseHeightRampForcing*>(&ext)),
    a2aDict_(),
    airLev_(),
    airI0_(-1),
    airI1_(-1),
    airW0_(0),
    airW1_(0)
{
    const dictionary& wp =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties");

    // Guard distance of the interface blend (document Rev4, rule 5)
    const scalar dzInterface = interfaceCellHeight(mesh);
    guard_ = max(delta_, 3.0*dzInterface);

    Info<< "SWENSE guard: delta = " << delta_ << " m, dz at interface = "
        << dzInterface << " m, guard d_g = " << guard_ << " m" << endl;

    // Air velocity mode of the blend (airFill); layer limits in metres
    if (airFill_ != "air2air" && airFill_ != "water" && airFill_ != "layer")
    {
        FatalIOErrorInFunction(wp)
            << "airFill must be air2air, water or layer, not " << airFill_
            << exit(FatalIOError);
    }

    if (airFill_ == "layer" && !(matchEnd_ > matchStart_))
    {
        FatalIOErrorInFunction(wp)
            << "layer airFill needs matchEnd > matchStart, got "
            << matchStart_ << " and " << matchEnd_
            << exit(FatalIOError);
    }

    Info<< "SWENSE airFill: " << airFill_ << ", matchStart = " << matchStart_
        << " m, matchEnd = " << matchEnd_ << " m" << endl;

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

        if (source == "theory" && hr_)
        {
            // rampMode height: one Air2Air per level of the height ramp,
            // built from that level's theory when it becomes active
            // (evaluate()), interpolated with the water-side weights
            a2aDict_ = a2a;
            airLev_.setSize(hr_->nLevels());

            Info<< "SWENSE Air2Air source theory, rampMode height: "
                << hr_->nLevels() << " levels, each from its own theory"
                << " (eta of that height), interpolated in H as the water"
                << endl;
        }
        else if (hr_)
        {
            FatalErrorInFunction
                << "rampMode height needs Air2Air source theory (found "
                << source << ")" << exit(FatalError);
        }
        else if (source == "theory")
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

    // Grid2Grid eta depends on (x, y) only: one row per distinct (x, y).
    // Hash map: the same rows in the same order (first appearance) as an
    // ordered map, without the O(n log n) tree
    struct keyHash
    {
        size_t operator()(const std::pair<long long, long long>& k) const
        {
            return std::hash<long long>()(k.first)*1000003u
                ^ std::hash<long long>()(k.second);
        }
    };
    std::unordered_map<std::pair<long long, long long>, label, keyHash> rowOf;
    rowOf.reserve(2*n);
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
        airTables(air_());
    }

    forAll(airLev_, i)
    {
        if (airLev_.set(i))
        {
            airTables(airLev_[i]);
        }
    }
}


void Foam::swenseIncident::airTables(Air2AirWrapper& a)
{
    const label n = pts_.size();
    scalarList xs(n), ys(n), zs(n);

    forAll(pts_, ip)
    {
        const vector xw = toWave(pts_[ip] - refPos_);
        xs[ip] = xw.x();
        ys[ip] = xw.y();
        zs[ip] = pts_[ip].z() - seaLevel_;
    }

    // The table rows depend on the points only: airIx_, airIy_, airIz_ are
    // the same for every wrapper
    a.solver().buildTables(xs, ys, zs, airIx_, airIy_, airIz_);
}


void Foam::swenseIncident::makeAirLevel(const label i)
{
    // As the source theory set-up of the constructor, with the theory of
    // level i (full values at its own height, not ramped)
    const waveTheories::waveTheory& th = hr_->level(i);
    const scalar c = th.phaseSpeed();
    const scalar K = th.waveNumber();

    if (c <= 0 || K <= 0)
    {
        FatalErrorInFunction
            << "Air2Air source theory: the wave theory gives no"
            << " phase speed or wave number" << exit(FatalError);
    }

    const scalar Lp = 2*M_PI/K;
    const label M = a2aDict_.lookupOrDefault<label>("nModes", 64);
    const label N = 4*M;

    std::vector<scalar> eta(N);
    for (label j = 0; j < N; ++j)
    {
        eta[j] = th.eta(point((j + 0.5)*Lp/N, 0, 0), 0.0);
    }

    Info<< "SWENSE Air2Air level " << i + 1 << " (H " << hr_->levelHeight(i)
        << " m)" << endl;

    airLev_.set(i, new Air2AirWrapper("", a2aDict_));
    airLev_[i].setTheory
    (
        1.0/K, magG_, c, Lp, eta, M,
        a2aDict_.lookupOrDefault<word>("cutoff", word("auto")),
        a2aDict_.lookupOrDefault<scalar>("cutoffC", 0.5)
    );

    airTables(airLev_[i]);
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

    if (airLev_.size())
    {
        // rampMode height: Air2Air of the two levels, weights of the water
        // side; no factor
        vector Uw = vector::zero;
        scalar Pd = 0;

        const label il[2] = {airI0_, airI1_};
        const scalar wl[2] = {airW0_, airW1_};

        for (label k = 0; k < 2; ++k)
        {
            if (il[k] >= 0)
            {
                vector Ul;
                scalar Pl;
                airLev_[il[k]].solver().flowTables
                (
                    airIx_[ip], airIy_[ip], airIz_[ip], Ul, Pl
                );
                Uw += wl[k]*Ul;
                Pd += wl[k]*Pl;
            }
        }

        p = Pd - rhoAir_*magG_*zw;

        return toMesh(Uw);
    }

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


Foam::vector Foam::swenseIncident::airVelocity
(
    const label ip,
    const scalar t,
    const vector& UAir
) const
{
    if (airFill_ == "air2air")
    {
        return UAir;
    }

    if (airFill_ == "water")
    {
        return waterUAt(ip, t);
    }

    // layer: smoothstep s(t) between the water velocity below matchStart
    // and the Air2Air velocity above matchEnd, zeta = z - eta
    const scalar zeta = pts_[ip].z() - etaAt(ip, t);
    const scalar tt =
        min(max((zeta - matchStart_)/(matchEnd_ - matchStart_), scalar(0)), scalar(1));
    const scalar s = 3.0*tt*tt - 2.0*tt*tt*tt;

    return (1.0 - s)*waterUAt(ip, t) + s*UAir;
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


Foam::vector Foam::swenseIncident::waterUAt
(
    const label ip,
    const scalar t
) const
{
    return fastEval_ ? UPt_[ip] : waterU(pts_[ip], t);
}


Foam::scalar Foam::swenseIncident::waterPAt
(
    const label ip,
    const scalar t
) const
{
    return
        fastEval_
      ? rhoWater_*(pExPt_[ip] - magG_*(pts_[ip] - refPos_).z())
      : waterP(pts_[ip], t);
}


void Foam::swenseIncident::countDrift
(
    const label ip,
    const scalar t,
    const bool useWater,
    const bool useAir
) const
{
    // Log only: a water value used above eta + d_g or an air value used
    // below eta - d_g means that alphaI has left the surface of the wave model
    const scalar zRel = pts_[ip].z() - etaAt(ip, t);

    if ((useWater && zRel > guard_) || (useAir && zRel < -guard_))
    {
        ++driftCount;

        if (mag(zRel) > driftWorstZ)
        {
            driftWorstZ = mag(zRel);
            driftWorstX = pts_[ip];
        }
    }
}


void Foam::swenseIncident::blendCell
(
    const label ip,
    const scalar t,
    const scalar gamma,
    vector& u,
    scalar& p
) const
{
    // Momentum blend (spec Definition 1, step 1b):
    //   rhoI = gamma rho_w + (1 - gamma) rho_a
    //   mI   = gamma rho_w uW + (1 - gamma) rho_a uA,   uI = mI/rhoI
    // The water field is used only where gamma > 0 and the air field only
    // where gamma < 1 (M2). The pressure is for output only.
    const bool useWater = gamma > 0;
    const bool useAir = gamma < 1;

    vector mI(Zero);
    p = 0.0;

    if (useWater)
    {
        mI += gamma*rhoWater_*waterUAt(ip, t);
        p += gamma*waterPAt(ip, t);
    }

    if (useAir)
    {
        scalar pa = 0.0;
        const vector Ua = airVelocity(ip, t, airUP(ip, pa));
        mI += (1.0 - gamma)*rhoAir_*Ua;
        p += (1.0 - gamma)*pa;
    }

    u = mI/(gamma*rhoWater_ + (1.0 - gamma)*rhoAir_);

    countDrift(ip, t, useWater, useAir);
}


Foam::vector Foam::swenseIncident::blendFace
(
    const label ip,
    const scalar t,
    const scalar af
) const
{
    // Raw face velocity (spec step 1c): mass-fraction weight of the face
    //   W_f = a_f rho_w/(a_f rho_w + (1 - a_f) rho_a)
    //   u_f = W_f uW(x_f) + (1 - W_f) uA(x_f)
    // uW only where a_f > 0, uA only where a_f < 1
    const bool useWater = af > 0;
    const bool useAir = af < 1;

    const scalar Wf = af*rhoWater_/(af*rhoWater_ + (1.0 - af)*rhoAir_);

    vector u(Zero);

    if (useWater)
    {
        u += Wf*waterUAt(ip, t);
    }

    if (useAir)
    {
        scalar pa = 0.0;
        u += (1.0 - Wf)*airVelocity(ip, t, airUP(ip, pa));
    }

    countDrift(ip, t, useWater, useAir);

    return u;
}



void Foam::swenseIncident::sidePressure
(
    volScalarField& pHyd,
    volScalarField& pDyn
) const
{
    const scalar t = tNow_;

    auto split = [&](const label ip, scalar& ph, scalar& pd)
    {
        const scalar zw = (pts_[ip] - refPos_).z();
        const scalar eta = etaAt(ip, t);

        scalar ps = 0.0;
        scalar rs = 0.0;

        if (zw < eta)
        {
            ps = waterPAt(ip, t);
            rs = rhoWater_;
        }
        else
        {
            airUP(ip, ps);
            rs = rhoAir_;
        }

        ph = rs*magG_*(eta - zw);
        pd = ps - ph;
    };

    scalarField& ph = pHyd.primitiveFieldRef();
    scalarField& pd = pDyn.primitiveFieldRef();

    forAll(ph, celli)
    {
        split(celli, ph[celli], pd[celli]);
    }

    label ip = offCb_;

    forAll(pHyd.boundaryField(), patchi)
    {
        const label nf = mesh_.C().boundaryField()[patchi].size();

        for (label facei = 0; facei < nf; ++facei)
        {
            scalar a = 0.0, b = 0.0;
            split(ip, a, b);

            if (!pHyd.boundaryField()[patchi].coupled())
            {
                pHyd.boundaryFieldRef()[patchi][facei] = a;
                pDyn.boundaryFieldRef()[patchi][facei] = b;
            }
            ++ip;
        }
    }

    pHyd.correctBoundaryConditions();
    pDyn.correctBoundaryConditions();
}


bool Foam::swenseIncident::alphaHExact(const objectRegistry& db)
{
    return
        db.foundObject<IOdictionary>("waveProperties")
     && db.lookupObject<IOdictionary>("waveProperties")
            .getOrDefault<Switch>("alphaHExact", false);
}


Foam::tmp<Foam::volScalarField>
Foam::swenseIncident::exactAlpha(const scalar t) const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject("alphaExact", mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("alphaExact", dimless, 0.0)
        )
    );
    volScalarField& alpha = talpha.ref();

    swenseSurface surf(mesh_, ext_, t);

    // Signed distance d = z - eta at every mesh point, as the cut uses it
    const pointField& P = mesh_.points();
    scalarField d(P.size());
    {
        const scalar Tsoft = swenseRamp::Tsoft(mesh_.thisDb());
        const scalar f = Tsoft > 0 ? swenseRamp::factor(Tsoft, t) : 1.0;
        forAll(P, pi)
        {
            d[pi] = P[pi].z() - f*ext_.eta(P[pi], t);
        }
    }

    const labelListList& cP = mesh_.cellPoints();
    const scalarField& V = mesh_.V();

    forAll(alpha, c)
    {
        scalar dMin = GREAT, dMax = -GREAT;
        for (const label pi : cP[c])
        {
            dMin = min(dMin, d[pi]);
            dMax = max(dMax, d[pi]);
        }
        if (dMax <= 0)
        {
            alpha[c] = 1.0;
        }
        else if (dMin >= 0)
        {
            alpha[c] = 0.0;
        }
        else
        {
            alpha[c] = min(max(surf.wetVolume(c)/V[c], scalar(0)), scalar(1));
        }
    }

    // Boundary faces: wetted face fraction
    const faceList& F = mesh_.faces();
    forAll(alpha.boundaryField(), patchi)
    {
        fvPatchScalarField& ab = alpha.boundaryFieldRef()[patchi];
        const label start = mesh_.boundary()[patchi].start();
        forAll(ab, i)
        {
            const label facei = start + i;
            scalar dMin = GREAT, dMax = -GREAT;
            for (const label pi : F[facei])
            {
                dMin = min(dMin, d[pi]);
                dMax = max(dMax, d[pi]);
            }
            if (dMax <= 0)
            {
                ab[i] = 1.0;
            }
            else if (dMin >= 0)
            {
                ab[i] = 0.0;
            }
            else
            {
                const scalar A = F[facei].mag(P);
                ab[i] = A > VSMALL
                  ? min(max(surf.wetArea(facei)/A, scalar(0)), scalar(1))
                  : 0.0;
            }
        }
    }

    return talpha;
}


Foam::tmp<Foam::volScalarField>
Foam::swenseIncident::targetAlpha(const scalar t) const
{
    if (alphaHExact(mesh_.thisDb()))
    {
        return exactAlpha(t);
    }
    return tmp<volScalarField>(new volScalarField(blendAlpha(t)));
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


bool Foam::swenseIncident::centresMoved() const
{
    // Same order as buildPointTables: cells, boundary cell values, faces,
    // boundary faces
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();

    label ip = 0;
    const label n = pts_.size();

    forAll(C, celli)
    {
        if (ip >= n || C[celli] != pts_[ip++]) return true;
    }
    forAll(C.boundaryField(), patchi)
    {
        forAll(C.boundaryField()[patchi], facei)
        {
            if (ip >= n || C.boundaryField()[patchi][facei] != pts_[ip++]) return true;
        }
    }
    forAll(Cf, facei)
    {
        if (ip >= n || Cf[facei] != pts_[ip++]) return true;
    }
    forAll(Cf.boundaryField(), patchi)
    {
        forAll(Cf.boundaryField()[patchi], facei)
        {
            if (ip >= n || Cf.boundaryField()[patchi][facei] != pts_[ip++]) return true;
        }
    }

    return ip != n;
}


void Foam::swenseIncident::movePoints()
{
    // The point lists follow the current cell and face centres. The
    // Grid2Grid basis is stored for these points, so it is set up again only
    // when a point has moved (a zero motion keeps the stored basis, as static).
    // No centre moved: the tables are already those of the current points.
    if (!centresMoved())
    {
        return;
    }

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


void Foam::swenseIncident::evaluate(const scalar t)
{
    // Grid2Grid time is t + referenceTime; Air2Air uses the same time
    const scalar tAir = t + refTime_;
    tNow_ = t;

    if (air_.valid())
    {
        air_->update(tAir);
    }

    if (airLev_.size())
    {
        // The two levels of the water-side interpolation at time t; others
        // are freed (the height ramp only moves up)
        hr_->weights(t, airI0_, airW0_, airI1_, airW1_);

        forAll(airLev_, i)
        {
            const bool active = (i == airI0_ || i == airI1_);

            if (active && !airLev_.set(i))
            {
                makeAirLevel(i);
            }
            else if (!active && airLev_.set(i))
            {
                airLev_.set(i, nullptr);
            }

            if (active)
            {
                airLev_[i].update(tAir);
            }
        }
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

    // One flat point per distinct (x, y) row, to read eta of the row
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
}


void Foam::swenseIncident::blendCells
(
    const labelUList& cells,
    const scalarField& gamma,
    vectorField& u,
    scalarField& p
) const
{
    const scalar t = tNow_;

    u.setSize(cells.size());
    p.setSize(cells.size());

    forAll(cells, i)
    {
        // Flat point index of a cell centre is the cell label (as in blend)
        blendCell
        (
            cells[i],
            t,
            min(max(gamma[i], scalar(0)), scalar(1)),
            u[i],
            p[i]
        );
    }
}


void Foam::swenseIncident::blend(const volScalarField& alphaI)
{
    const scalar t = tNow_;

    driftCount = 0;
    driftWorstZ = 0.0;
    driftWorstX = Foam::point::zero;

    // Incident water fraction gamma = alphaI limited to [0, 1], the same
    // limit as the density function swenseRho
    const volScalarField gamma
    (
        "gammaI",
        min(max(alphaI, scalar(0)), scalar(1))
    );

    // Non-const access first: stores the old-time value of UInc before any
    // cell value is changed
    vectorField& Ui = UInc_.primitiveFieldRef();
    scalarField& pi = pInc_.primitiveFieldRef();

    // Cells (step 1b)
    forAll(Ui, celli)
    {
        blendCell(celli, t, gamma[celli], Ui[celli], pi[celli]);
    }

    // Logging only: largest |water velocity - air velocity handed to the
    // blend| in cells with 0 < alphaI < 1 (no effect on the fields)
    scalar gapMax = 0.0;

    forAll(Ui, celli)
    {
        if (alphaI[celli] > 0 && alphaI[celli] < 1)
        {
            scalar pa = 0.0;
            const vector Ua = airVelocity(celli, t, airUP(celli, pa));
            gapMax = max(gapMax, mag(waterUAt(celli, t) - Ua));
        }
    }

    // Boundary values of the cell fields: same blend at the patch face
    // centres with the boundary value of gamma
    label ip = offCb_;

    forAll(gamma.boundaryField(), patchi)
    {
        const scalarField& gb = gamma.boundaryField()[patchi];

        forAll(gb, facei)
        {
            blendCell
            (
                ip, t, gb[facei],
                UInc_.boundaryFieldRef()[patchi][facei],
                pInc_.boundaryFieldRef()[patchi][facei]
            );
            ++ip;
        }
    }

    // Faces (step 1c): a_f is the linear face value of gamma
    const surfaceScalarField af("afI", linearInterpolate(gamma));

    vectorField& Ufi = UfInc_.primitiveFieldRef();

    forAll(Ufi, facei)
    {
        Ufi[facei] = blendFace(offCf_ + facei, t, af[facei]);
    }

    ip = offCfb_;

    forAll(af.boundaryField(), patchi)
    {
        const fvsPatchScalarField& ab = af.boundaryField()[patchi];

        forAll(ab, facei)
        {
            UfInc_.boundaryFieldRef()[patchi][facei] =
                blendFace(ip, t, ab[facei]);
            ++ip;
        }
    }

    // One line for the points where a wave model was used outside its guard
    Info<< "SWENSE drift: " << driftCount
        << " points used outside the guard, worst |z - eta| = "
        << driftWorstZ << " m at " << driftWorstX << endl;

    Info<< "SWENSE gap: max|water - air handed to blend| in 0<alphaI<1 = "
        << gapMax << endl;
}

// ************************************************************************* //
