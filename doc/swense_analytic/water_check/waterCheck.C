// Phase 3 water checks for waveTheoryForcing (streamFunction, stokesFifth, stokesSecond)
#include "fvCFD.H"
#include "externalWaveForcing.H"
#include "IOdictionary.H"
#include <cmath>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary waveProperties
    (
        IOobject("waveProperties", runTime.constant(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE)
    );

    autoPtr<waveTheories::externalWaveForcing> extP =
        waveTheories::externalWaveForcing::New
        (
            IOobject("externalWaveForcing", runTime.constant(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            runTime, mesh
        );

    const waveTheories::externalWaveForcing& ext = extP();

    const dictionary& cd = waveProperties.subDict("waveTheoryForcingCoeffs");
    const word type(cd.lookup("waveType"));
    const scalar H = readScalar(cd.lookup("height"));
    const scalar h = readScalar(cd.lookup("depth"));
    const scalar omega = readScalar(cd.lookup("omega"));
    const scalar k = vector(cd.lookup("waveNumber")).x();
    const scalar T = 2*M_PI/omega;
    const scalar a = 0.5*H;
    const scalar g = 9.81;
    const scalar lambda = 2*M_PI/k;

    Info<< "theory " << type << "  H " << H << " h " << h << " T " << T
        << " k " << k << " omega " << omega << " a " << a << nl;

    // ---- 1. divergence at 50 points in the water ----
    const scalar dl = 1e-4;
    scalar maxDiv = 0.0;
    for (label i = 0; i < 50; ++i)
    {
        const scalar x = (i + 0.5)*lambda/50;
        const scalar z = -4.0 + 3.0*((37*i) % 50)/50.0;   // water, below eta
        const scalar t = 0.25*T;   // theory eta is zero for t < Tstart = 0 (central difference needs t > 0)
        const point p(x, 0, z);
        const scalar dudx = (ext.U(point(x + dl, 0, z), t).x()
                              - ext.U(point(x - dl, 0, z), t).x())/(2*dl);
        const scalar dwdz = (ext.U(point(x, 0, z + dl), t).z()
                              - ext.U(point(x, 0, z - dl), t).z())/(2*dl);
        maxDiv = max(maxDiv, mag(dudx + dwdz));
    }
    const scalar scaleKAW = k*a*omega;
    Info<< "CHECK div: max|div U| = " << maxDiv << "  (k a w = " << scaleKAW
        << ", ratio " << maxDiv/scaleKAW << ", rule 1e-6)" << nl;

    // ---- 2. kinematic FSBC at z = eta, 200 points over one period ----
    const scalar tau = T/1e5;
    const scalar dx = 1e-4;
    scalar maxKin = 0.0;
    for (label j = 0; j < 200; ++j)
    {
        const scalar x = (j + 0.5)*lambda/200;
        const scalar t = 0.25*T;   // theory eta is zero for t < Tstart = 0 (central difference needs t > 0)
        const scalar eta = ext.eta(point(x, 0, 0), t);
        const scalar etaT = (ext.eta(point(x, 0, 0), t + tau)
                             - ext.eta(point(x, 0, 0), t - tau))/(2*tau);
        const scalar etaX = (ext.eta(point(x + dx, 0, 0), t)
                             - ext.eta(point(x - dx, 0, 0), t))/(2*dx);
        const point ps(x, 0, eta);
        const vector U = ext.U(ps, t);
        const scalar R = etaT + U.x()*etaX - U.z();
        maxKin = max(maxKin, mag(R));
    }
    Info<< "CHECK kin: max|eta_t + u eta_x - w| at eta = " << maxKin
        << "  (a w = " << a*omega << ", ratio " << maxKin/(a*omega)
        << ", rule 1e-6)" << nl;

    // ---- 3. dynamic FSBC: total pressure p/rho = pExcess - g eta at eta ----
    scalar maxDyn = 0.0;
    for (label j = 0; j < 200; ++j)
    {
        const scalar x = (j + 0.5)*lambda/200;
        const scalar t = 0.25*T;   // theory eta is zero for t < Tstart = 0 (central difference needs t > 0)
        const scalar eta = ext.eta(point(x, 0, 0), t);
        const scalar p = ext.pExcess(point(x, 0, eta), t) - g*eta;
        maxDyn = max(maxDyn, mag(p));
    }
    Info<< "CHECK dyn: max|p/rho| at eta = " << maxDyn
        << "  (g a = " << g*a << ", ratio " << maxDyn/(g*a)
        << ", rule 1e-5)" << nl;

    Info<< "End" << nl;
    return 0;
}
