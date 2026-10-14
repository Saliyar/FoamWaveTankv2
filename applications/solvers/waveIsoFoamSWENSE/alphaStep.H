// isoAdvector transport of alpha1 (total flux phi) and alphaI (incident flux
// phiInc) with the same isoAdvection settings, in every outer corrector
// (document Rev1, 16.3.2). Both restart from their old-time values.
#include "rhofs.H"

{
    #include "alphaSuSp.H"

    alpha1 = alpha1.oldTime();
    advector.advect(Sp, Su);
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
}

{
    #include "alphaSuSp.H"

    alphaI = alphaI.oldTime();
    advectorI.advect(Sp, Su);
    rhoPhiI = advectorI.getRhoPhi(rho1f, rho2f);
}
