// isoAdvector transport of alpha1 (total flux phi) and alphaI (incident flux
// phiInc) with the same isoAdvection settings, in every outer corrector
// (document Rev4, 16.3.2 and 16.7). Both restart from their old-time values.
// Advective form: the source is Sp = div(flux) of the same flux that advects
// the field, Su = 0, as the MULES alpha step of waveFoamSWENSE.
#include "rhofs.H"
// Row 4 fix (Rev 5 audit): both advections use the current total velocity
U == incident.UInc() + Uc;

{
    volScalarField Sp("Sp", fvc::div(phi));
    zeroField Su;

    alpha1 = alpha1.oldTime();
    advector.advect(Sp, Su);
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
}

{
    volScalarField SpI("SpI", fvc::div(phiIRel));
    zeroField Su;

    alphaI = alphaI.oldTime();
    advectorI.advect(SpI, Su);
    rhoPhiI = advectorI.getRhoPhi(rho1f, rho2f);
}
