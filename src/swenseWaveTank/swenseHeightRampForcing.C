/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseHeightRampForcing.H"
#include "swenseRamp.H"
#include "setWaveProperties.H"
#include "gravityMeshObject.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "Pstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "mathematicalConstants.H"
#include "swenseIncidentResidual.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swenseHeightRampForcing, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swenseHeightRampForcing::swenseHeightRampForcing
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    waveTheories::waveTheoryForcing(io, rT, mesh),
    H_(0),
    Tramp_(swenseRamp::rampTime(mesh.thisDb())),
    Hlev_(),
    lev_()
{
    const label n = 20;

    // The registered waveProperties: the per-level coefficient dictionaries
    // are added to it while the theory is built (waveTheory reads its
    // coefficients from the registered dictionary) and removed after
    IOdictionary& wp = const_cast<IOdictionary&>
    (
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
    );

    const dictionary base(wp.subDict("waveTheoryForcingCoeffs"));
    H_ = base.get<scalar>("height");

    const vector g(meshObjects::gravity::New(rT).value());

    Hlev_.setSize(n);
    lev_.setSize(n);

    Info<< "swenseHeightRampForcing: rampMode height, rampTsoft " << Tramp_
        << " s, " << word(base.lookup("waveType")) << ", H " << H_ << " m, "
        << n << " levels" << endl;

    for (label i = 0; i < n; ++i)
    {
        Hlev_[i] = H_*scalar(i + 1)/scalar(n);

        dictionary d(base);
        d.set("height", Hlev_[i]);
        d.name() = "waveTheoryForcingCoeffs";

        // setWaveProperties on the master only, one height at a time
        string txt;

        if (Pstream::master())
        {
            const bool oldPar = UPstream::parRun(false);

            OStringStream os;
            autoPtr<setWaveProperties> sp
            (
                setWaveProperties::New(rT, d, g, true)
            );
            sp->set(os);
            txt = os.str();

            UPstream::parRun(oldPar);
        }

        Pstream::broadcast(txt);

        IStringStream is(txt);
        const dictionary parsed(is);

        dictionary coeffs;
        bool found = false;

        for (const entry& e : parsed)
        {
            if (e.isDict())
            {
                coeffs = e.dict();
                found = true;
                break;
            }
        }

        if (!found)
        {
            FatalErrorInFunction
                << "setWaveProperties gave no coefficient dictionary for height "
                << Hlev_[i] << nl << txt << exit(FatalError);
        }

        // No soft start inside the theory: the height ramp is the only ramp
        coeffs.set("Tsoft", scalar(0));

        const word key("swenseHeightRamp" + Foam::name(i + 1));

        wp.set(word(key + "Coeffs"), coeffs);
        lev_.set(i, waveTheories::waveTheory::New(key, mesh).ptr());
        wp.remove(word(key + "Coeffs"));

        Info<< "    level " << i + 1 << ": H " << Hlev_[i] << " m, k "
            << lev_[i].waveNumber() << " 1/m, c " << lev_[i].phaseSpeed()
            << " m/s" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::swenseHeightRampForcing::rampedHeight(const scalar t) const
{
    if (Tramp_ <= 0)
    {
        return H_;
    }

    return H_*swenseRamp::factor(Tramp_, t);
}


void Foam::swenseHeightRampForcing::weights
(
    const scalar t,
    label& i0,
    scalar& w0,
    label& i1,
    scalar& w1
) const
{
    const label n = lev_.size();

    // Position of Hr on the level grid: level j (1..n) has height H j/n
    const scalar r = rampedHeight(t)/H_*scalar(n);

    if (r >= scalar(n))
    {
        i0 = n - 1;
        w0 = 1;
        i1 = -1;
        w1 = 0;
        return;
    }

    const label j = label(std::floor(r));
    const scalar frac = r - scalar(j);

    // Lower neighbour level j (index j - 1; j = 0 is still water), upper
    // neighbour level j + 1 (index j)
    i0 = j - 1;
    w0 = 1 - frac;
    i1 = j;
    w1 = frac;

    if (i0 < 0)
    {
        w0 = 0;
    }

    if (w1 == 0)
    {
        i1 = -1;
    }
}


Foam::scalar Foam::swenseHeightRampForcing::eta
(
    const point& x,
    const scalar& time
) const
{
    label i0, i1;
    scalar w0, w1;
    weights(time, i0, w0, i1, w1);

    scalar e = 0;
    if (i0 >= 0) e += w0*lev_[i0].eta(x, time);
    if (i1 >= 0) e += w1*lev_[i1].eta(x, time);

    return e;
}


Foam::scalar Foam::swenseHeightRampForcing::pExcess
(
    const point& x,
    const scalar& time
) const
{
    label i0, i1;
    scalar w0, w1;
    weights(time, i0, w0, i1, w1);

    // As waveTheoryForcing: the theory returns rho (p/rho + g z)
    scalar p = 0;
    if (i0 >= 0) p += w0*lev_[i0].pExcess(x, time);
    if (i1 >= 0) p += w1*lev_[i1].pExcess(x, time);

    return p/rhoWater_;
}


Foam::vector Foam::swenseHeightRampForcing::U
(
    const point& x,
    const scalar& time
) const
{
    label i0, i1;
    scalar w0, w1;
    weights(time, i0, w0, i1, w1);

    vector u = vector::zero;
    if (i0 >= 0) u += w0*lev_[i0].U(x, time);
    if (i1 >= 0) u += w1*lev_[i1].U(x, time);

    return u;
}


// * * * * * * * * * * * * * * Static Set-up Functions * * * * * * * * * * //

Foam::autoPtr<Foam::waveTheories::externalWaveForcing>
Foam::swenseHeightRampForcing::newForcing(Time& rT, const fvMesh& mesh)
{
    const IOobject extIO
    (
        "externalWaveForcing",
        rT.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (swenseRamp::heightMode(mesh.thisDb()))
    {
        const word extType
        (
            mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
                .lookupOrDefault<word>("externalForcing", "emptyExternal")
        );

        if (extType != "waveTheoryForcing")
        {
            FatalErrorInFunction
                << "rampMode height needs externalForcing waveTheoryForcing"
                << " (found " << extType << ")" << exit(FatalError);
        }

        return autoPtr<waveTheories::externalWaveForcing>
        (
            new swenseHeightRampForcing(extIO, rT, mesh)
        );
    }

    Info << "\nReading waveProperties" << endl;

    autoPtr<waveTheories::externalWaveForcing> ext
    (
        waveTheories::externalWaveForcing::New(extIO, rT, mesh)
    );

    // SWENSE takes the incident pressure from the forcing (p = pInc + pc +
    // Drho gh): a forcing whose pExcess is a placeholder 0 would run with
    // the hydrostatic incident pressure only. waveTheoryForcing refuses such
    // theories itself (cnoidalFirst); this catches external forcings
    // (irregularFast)
    if (!ext->implementPressure())
    {
        FatalErrorInFunction
            << "externalForcing " << ext->type() << " does not implement the"
            << " incident pressure (implementPressure false, pExcess = 0)."
            << nl << "    SWENSE needs pExcess for pInc. For a linear"
            << " irregular sea use externalForcing waveTheoryForcing with"
            << " waveType irregular (linear pressure)"
            << exit(FatalError);
    }

    return ext;
}


void Foam::swenseHeightRampForcing::residualCheck
(
    const waveTheories::externalWaveForcing& ext,
    const fvMesh& mesh
)
{
    const scalar resMax =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .lookupOrDefault<scalar>("incidentResidualMax", 0.005);

    const scalar zb = gMin(mesh.points().component(vector::Z)());

    const IOdictionary& wp =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties");

    const swenseHeightRampForcing* hr =
        dynamic_cast<const swenseHeightRampForcing*>(&ext);
    const waveTheories::waveTheoryForcing* wtf =
        dynamic_cast<const waveTheories::waveTheoryForcing*>(&ext);

    // c and K of the checked wave: the theory's own (phaseSpeed(),
    // waveNumber()), else (stokesFirst, irregular) linear dispersion with the
    // period of the coefficients (period, else the peak period Tp) and depth
    auto cK = [&](const scalar c0, const scalar K0, scalar& c, scalar& K)
    {
        if (c0 > 0 && K0 > 0)
        {
            c = c0;
            K = K0;
            return word("own c, K");
        }

        const dictionary& cd = wp.subDict("waveTheoryForcingCoeffs");
        const scalar T =
            cd.found("period") ? cd.get<scalar>("period") : cd.get<scalar>("Tp");
        const scalar d = cd.get<scalar>("depth");
        const scalar gm = mag(meshObjects::gravity::New(mesh.time()).value());
        const scalar om = constant::mathematical::twoPi/T;

        K = om*om/gm;
        for (label it = 0; it < 200; ++it)
        {
            K = om*om/(gm*Foam::tanh(K*d));
        }
        c = om/K;

        return word("linear c, K at T " + Foam::name(T));
    };

    scalar c = 0, K = 0;

    if (hr)
    {
        const waveTheories::waveTheory& th = hr->level(hr->nLevels() - 1);
        const word how(cK(th.phaseSpeed(), th.waveNumber(), c, K));
        swenseIncidentResidual::check
        (
            th, c, K, zb, resMax,
            "height ramp, " + word(th.type()) + ", level H "
          + Foam::name(hr->fullHeight()) + ", " + how
        );
    }
    else if (wtf)
    {
        const word how(cK(wtf->phaseSpeed(), wtf->waveNumber(), c, K));
        swenseIncidentResidual::check
        (
            *wtf, c, K, zb, resMax,
            "waveTheoryForcing, "
          + word(wp.subDict("waveTheoryForcingCoeffs").lookup("waveType"))
          + ", " + how
        );
    }
    else
    {
        Info<< "SWENSE incident residual: not checked (forcing "
            << ext.type() << ")" << endl;
    }

    // Optional: the incident elevation as SWENSE uses it (ramp included) at
    // given points, written once at start-up for startTime..endTime, for
    // comparison with the gauges (diag/s14.py). incidentEtaProbes (points);
    // incidentEtaProbesDt (default 0.01 s)
    if (wp.found("incidentEtaProbes") && Pstream::master())
    {
        const List<point> pr(wp.lookup("incidentEtaProbes"));
        const scalar dtS = wp.lookupOrDefault<scalar>("incidentEtaProbesDt", 0.01);
        const scalar Ts = swenseRamp::Tsoft(mesh.thisDb());
        const Time& rt = mesh.time();

        const fileName dir(rt.globalPath()/"postProcessing"/"incidentEta");
        mkDir(dir);
        OFstream os(dir/"eta.dat");
        os.precision(10);

        os << "# incident eta (ramp included), points";
        for (const point& q : pr)
        {
            os << " " << q;
        }
        os << nl;

        const label nt = label((rt.endTime().value() - rt.startTime().value())/dtS + 0.5);

        for (label i = 0; i <= nt; ++i)
        {
            const scalar t = rt.startTime().value() + i*dtS;
            os << t;
            for (const point& q : pr)
            {
                os << " " << ext.eta(q, t)*swenseRamp::factor(Ts, t);
            }
            os << nl;
        }

        Info<< "SWENSE incident eta at " << pr.size() << " points written to "
            << dir/"eta.dat" << endl;
    }
}


// ************************************************************************* //
