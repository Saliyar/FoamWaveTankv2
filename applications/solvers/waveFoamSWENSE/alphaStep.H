// VoF transport of alpha1 (total flux phi) and alphaI (incident flux phiInc)
// by the same function, in every outer corrector (document Rev1, 16.3.2).
// Both restart from their old-time values inside swenseAdvanceVoF.
// On a moving mesh the fluxes advect relative to the mesh (phi, phiIRel) and
// the source is the divergence of the absolute flux (phi + meshPhi). On a
// static mesh this is the production call, unchanged.
if (mesh.moving())
{
    const volScalarField SpAlpha1("SpAlpha1", fvc::div(phi + meshPhi));
    const volScalarField SpAlphaI("SpAlphaI", fvc::div(phiIRel + meshPhi));

    swenseAdvanceVoF(alpha1, phi, SpAlpha1, U, rhoPhi, rho1, rho2);
    swenseAdvanceVoF(alphaI, phiIRel, SpAlphaI, UInc, rhoPhiI, rho1, rho2);
}
else
{
    swenseAdvanceVoF(alpha1, phi, U, rhoPhi, rho1, rho2);
    swenseAdvanceVoF(alphaI, phiInc, UInc, rhoPhiI, rho1, rho2);
}
