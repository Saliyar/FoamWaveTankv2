// Relaxation of alpha1 and alphaI toward the same HOS target alphaH with the
// same weight (document Rev1, 16.3.3), then the densities rho and rhoI from
// the same function (16.1), the density ratio (16.2) and the incident part
// of the pressure pStar = (rho/rhoI) pI. Applied after the alpha transport
// of each outer corrector; both fields restart from their old-time values
// in every corrector, so the relaxation is not compounded.

alpha1 = (1.0 - w)*alpha1 + w*alphaH;
alphaI = (1.0 - w)*alphaI + w*alphaH;

// Level correction of the incident surface toward the HOS level eta_H
#include "swenseColumnLevel.H"

rho == swenseRho(alpha1, rho1, rho2);
rhoI == swenseRho(alphaI, rho1, rho2);
rhoByRhoI == rho/rhoI;

alpha2 = 1.0 - alpha1;
mixture.correct();

pStar == rhoByRhoI*incident.pInc();
