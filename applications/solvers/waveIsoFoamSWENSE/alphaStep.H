// Volume fractions by the isoAdvector routine, one set of settings
// (momentum_swense_Rev1, steps 1a and 2, M4). Once per time step, before the
// incident field of the new level is built, so the fluxes and velocities
// are those of level n:
//   1a  alphaI with phiIRel (= phiInc of level n) and UInc of level n
//   2   alpha1 with phi = phiInc + phiComp and U = UInc + Uc of level n
// Both calls are advector.advect(Sp, Su) with Sp = div(flux), Su = 0, as
// before. rhoPhiI and rhoPhi are the mass fluxes of the two transports.
#include "rhofs.H"

U == incident.UInc() + Uc;

{
    volScalarField SpI("SpI", fvc::div(phiIRel));
    zeroField Su;

    alphaI = alphaI.oldTime();

    advectorI.advect(SpI, Su);

    rhoPhiI = advectorI.getRhoPhi(rho1f, rho2f);
}

{
    volScalarField Sp("Sp", fvc::div(phi));
    zeroField Su;

    alpha1 = alpha1.oldTime();

    advector.advect(Sp, Su);

    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
}
