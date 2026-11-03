// isoAdvector transport of alpha1 (total flux phi) and alphaI (incident flux
// phiInc) with the same isoAdvection settings, in every outer corrector
// (document Rev4, 16.3.2 and 16.7). Both restart from their old-time values.
// Advective form: the source is Sp = div(flux), Su = 0, as the MULES alpha
// step of waveFoamSWENSE. On a moving mesh the source is the divergence of the
// ABSOLUTE flux (phi + meshPhi for alpha1, phiIRel + meshPhi for alphaI) and
// the advected velocities are relative to the mesh: U - reconstruct(meshPhi)
// and UInc - reconstruct(meshPhi), restored after the step, as waveIsoFoam
// alphaEqn.H does. On a static mesh these are the production expressions.
#include "rhofs.H"
// Row 4 fix (Rev 5 audit): both advections use the current total velocity
U == incident.UInc() + Uc;

{
    volScalarField Sp("Sp", fvc::div(phi));
    zeroField Su;

    if (mesh.moving())
    {
        Sp = fvc::div(phi + meshPhi);
    }

    alpha1 = alpha1.oldTime();

    // isoAdvection::advect() reads alpha1.oldTime() unscaled on a moving mesh
    // (its V0/V scaling only reaches alpha1In_, which is then overwritten), so
    // the old value is scaled to the current cell volume here and restored
    // after the step, as FloatStepper alphaEqn.H does (upstream b118249).
    scalarField alpha1OldValues;

    if (mesh.moving())
    {
        alpha1OldValues = alpha1.oldTime().primitiveField();
        alpha1.oldTime().primitiveFieldRef() *= (mesh.Vsc0()/mesh.Vsc());

        U -= fvc::reconstruct(meshPhi);
    }

    advector.advect(Sp, Su);

    if (mesh.moving())
    {
        alpha1.oldTime().primitiveFieldRef() = alpha1OldValues;

        U += fvc::reconstruct(meshPhi);
    }

    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
}

{
    volScalarField SpI("SpI", fvc::div(phiInc));
    zeroField Su;

    if (mesh.moving())
    {
        SpI = fvc::div(phiIRel + meshPhi);
    }

    alphaI = alphaI.oldTime();

    // Same V0/V scaling of the old value as for alpha1 above
    scalarField alphaIOldValues;

    if (mesh.moving())
    {
        alphaIOldValues = alphaI.oldTime().primitiveField();
        alphaI.oldTime().primitiveFieldRef() *= (mesh.Vsc0()/mesh.Vsc());

        incident.UIncRW() -= fvc::reconstruct(meshPhi);
    }

    advectorI.advect(SpI, Su);

    if (mesh.moving())
    {
        alphaI.oldTime().primitiveFieldRef() = alphaIOldValues;

        incident.UIncRW() += fvc::reconstruct(meshPhi);
    }

    rhoPhiI = advectorI.getRhoPhi(rho1f, rho2f);
}
