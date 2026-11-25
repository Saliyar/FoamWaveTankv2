// Relaxation of alpha1 and alphaI toward the same wave-model target alphaH
// with the same weight, then the densities rho and rhoI from the same
// function and their difference Drho (momentum_swense_Rev1, Section 12.2).
// The same operations act on both fractions, so alpha1 = alphaI is kept bit
// for bit. The column level correction (controlDict levelCorrection, default
// off) adds the same increment to both, then clips alpha1 only.

alpha1 = (1.0 - w)*alpha1 + w*alphaH;
alphaI = (1.0 - w)*alphaI + w*alphaH;

if (levelCorrection)
{
    #include "swenseColumnLevel.H"
}

rho == swenseRho(alpha1, rho1, rho2);
rhoI == swenseRho(alphaI, rho1, rho2);
Drho == rho - rhoI;

alpha2 = 1.0 - alpha1;
mixture.correct();
