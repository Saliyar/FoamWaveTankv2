// isoAdvector transport of alpha1 (total flux phi) and alphaI (incident flux
// phiInc) with the same isoAdvection settings, in every outer corrector
// (document Rev4, 16.3.2 and 16.7). Both restart from their old-time values.
// Advective form: the source is Sp = div of the ABSOLUTE flux (phi + meshPhi
// for alpha1, phiIRel + meshPhi = phiInc for alphaI), Su = 0, as the MULES
// alpha step of waveFoamSWENSE. The advected velocities are relative to the
// mesh: U - reconstruct(meshPhi) and UInc - reconstruct(meshPhi) on a moving
// mesh, restored after the step, as waveIsoFoam alphaEqn.H does.
#include "rhofs.H"
// Row 4 fix (Rev 5 audit): both advections use the current total velocity
U == incident.UInc() + Uc;

{
    volScalarField Sp("Sp", fvc::div(phi + meshPhi));
    zeroField Su;

    alpha1 = alpha1.oldTime();

    if (mesh.moving())
    {
        U -= fvc::reconstruct(meshPhi);
    }

    advector.advect(Sp, Su);

    if (mesh.moving())
    {
        U += fvc::reconstruct(meshPhi);
    }

    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
}

{
    volScalarField SpI("SpI", fvc::div(phiIRel + meshPhi));
    zeroField Su;

    alphaI = alphaI.oldTime();

    if (mesh.moving())
    {
        incident.UIncRW() -= fvc::reconstruct(meshPhi);
    }

    advectorI.advect(SpI, Su);

    if (mesh.moving())
    {
        incident.UIncRW() += fvc::reconstruct(meshPhi);
    }

    rhoPhiI = advectorI.getRhoPhi(rho1f, rho2f);
}
