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
    volScalarField SpI("SpI", fvc::div(swenseNoRelFlux ? phiInc : phiIRel));
    zeroField Su;

    alphaI = alphaI.oldTime();

    // incidentAlphaWall: body patches get the incident wetted fraction on
    // inflow faces (t^{n+1/2}, current face points) before the advection;
    // upwind(phi).flux(alphaI) takes this patch value on those faces.
    if (incidentAlphaWall)
    {
        label nIn = 0;
        scalar wIn = 0, inFlux = 0, outFlux = 0;

        forAll(bodyPatches, i)
        {
            const label patchi = bodyPatches[i];
            fvPatchScalarField& ab = alphaI.boundaryFieldRef()[patchi];
            ab.evaluate();

            const swenseIncidentAlphaFvPatchScalarField& sb =
                refCast<const swenseIncidentAlphaFvPatchScalarField>(ab);
            nIn += sb.nInflow();
            wIn += sb.wetInflow();

            const scalarField& pf =
                (swenseNoRelFlux ? phiInc : phiIRel).boundaryField()[patchi];
            forAll(pf, facei)
            {
                const scalar q = pf[facei]*ab[facei];
                if (pf[facei] < 0) { inFlux += q; } else { outFlux += q; }
            }
        }

        reduce(nIn, sumOp<label>());
        reduce(wIn, sumOp<scalar>());
        reduce(inFlux, sumOp<scalar>());
        reduce(outFlux, sumOp<scalar>());

        Info<< "SWENSE incidentAlphaWall: inflow faces " << nIn
            << " mean wet " << (nIn > 0 ? wIn/nIn : 0.0)
            << " alphaI flux in " << inFlux << " out " << outFlux
            << " m3/s" << endl;
    }

    advectorI.advect(SpI, Su);

    rhoPhiI = advectorI.getRhoPhi(rho1f, rho2f);
}

{
    // Absolute total flux for switch a (updated every step; unused otherwise)
    phiAbsTot = phi + meshPhi;

    volScalarField Sp("Sp", fvc::div(swenseNoRelFlux ? phiAbsTot : phi));
    zeroField Su;

    alpha1 = alpha1.oldTime();

    advector.advect(Sp, Su);

    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
}
