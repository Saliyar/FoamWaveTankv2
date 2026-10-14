// VoF transport of alpha1 (total flux phi) and alphaI (incident flux phiInc)
// by the same function, in every outer corrector (document Rev1, 16.3.2).
// Both restart from their old-time values inside swenseAdvanceVoF.
swenseAdvanceVoF(alpha1, phi, U, rhoPhi, rho1, rho2);
swenseAdvanceVoF(alphaI, phiInc, UInc, rhoPhiI, rho1, rho2);
