/*---------------------------------------------------------------------------*\
Description
    Start-up consistency check of an analytic incident wave (diag/s12.py
    method). Over one period (the second, t in [T, 2T)) at x = 0, which is
    one wavelength for a progressive wave along x, at the full height:

        kinematic   R_k = eta_t + u_s eta_x - w_s           at z = eta
        mass flux   R_m = eta_t + dQ/dx,  Q = int_{zb}^{eta} u dz

    both scaled by omega H/2 (H = crest - trough of the samples). Derivatives
    by central differences (dx = 1e-4 L, dt = 1e-4 T), Q by Simpson's rule
    (400 intervals). zb is the lowest mesh point. W is any object with
    eta(point, t) and U(point, t) in the mesh frame, not ramped.

SourceFiles
    swenseIncidentResidual.H (header only)

\*---------------------------------------------------------------------------*/

#ifndef swenseIncidentResidual_H
#define swenseIncidentResidual_H

#include "error.H"
#include "messageStream.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace swenseIncidentResidual
{

//- Max |R_k| and |R_m| over one period, scaled by omega H/2; prints one line
//  and stops above maxAllowed
template<class W>
void check
(
    const W& w,
    const scalar c,
    const scalar K,
    const scalar zb,
    const scalar maxAllowed,
    const string& what
)
{
    const scalar twoPi = constant::mathematical::twoPi;
    const scalar L = twoPi/K;
    const scalar T = L/c;
    const scalar om = twoPi/T;
    const scalar dx = 1e-4*L;
    const scalar dt = 1e-4*T;
    const label n = 256;
    const label nz = 400;

    auto etaAt = [&](const scalar x, const scalar t)
    {
        return w.eta(point(x, 0, 0), t);
    };

    // Q = int_{zb}^{eta} u dz, Simpson
    auto Q = [&](const scalar x, const scalar t)
    {
        const scalar e = etaAt(x, t);
        const scalar dz = (e - zb)/nz;
        scalar s = 0;

        for (label j = 0; j <= nz; ++j)
        {
            const scalar wt = (j == 0 || j == nz) ? 1 : (j % 2 ? 4 : 2);
            s += wt*w.U(point(x, 0, zb + j*dz), t).x();
        }

        return s*dz/3.0;
    };

    scalar Rk = 0, Rm = 0, eMax = -GREAT, eMin = GREAT;

    for (label i = 0; i < n; ++i)
    {
        // Second period: some theories are 0 before t = 0 (Tstart)
        const scalar t = T*(1 + scalar(i)/scalar(n));
        const scalar e = etaAt(0, t);
        const scalar et = (etaAt(0, t + dt) - etaAt(0, t - dt))/(2*dt);
        const scalar ex = (etaAt(dx, t) - etaAt(-dx, t))/(2*dx);
        const vector us = w.U(point(0, 0, e), t);

        Rk = max(Rk, mag(et + us.x()*ex - us.z()));
        Rm = max(Rm, mag(et + (Q(dx, t) - Q(-dx, t))/(2*dx)));
        eMax = max(eMax, e);
        eMin = min(eMin, e);
    }

    const scalar H = eMax - eMin;
    const scalar s = om*H/2;

    Info<< "SWENSE incident residual (" << what.c_str() << "): kinematic "
        << Rk/s << ", mass flux " << Rm/s << " of omega H/2 (H " << H
        << " m, T " << T << " s, L " << L << " m, bottom z " << zb
        << " m), incidentResidualMax " << maxAllowed << endl;

    if (max(Rk, Rm)/s > maxAllowed)
    {
        FatalErrorInFunction
            << "The incident wave (" << what.c_str() << ") does not satisfy"
            << " the free-surface kinematic condition: kinematic residual "
            << Rk/s << ", mass-flux residual " << Rm/s << " of omega H/2,"
            << " above incidentResidualMax " << maxAllowed << nl
            << "    use a wave theory that satisfies it (for example"
            << " streamFunction) or raise incidentResidualMax in"
            << " waveProperties" << exit(FatalError);
    }
}

}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
