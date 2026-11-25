// Incident reset in the body zone (waveProperties incidentResetInZone,
// default off). Included after the alphaI step, the relaxation and the
// incident step of the new level (swenseTimeLoop.H), so state A below is the
// unreset incident of level n+1:
//   alphaI_c = V_wet(eta(t^{n+1}))/V_c   in the body-zone cells, the exact
//              wetted fraction below the ramped incident surface on the
//              current cell points (moving mesh)
//   rhoI, Drho                           rebuilt from alphaI
//   uI, pInc (cell blend, mI/rhoI)       rebuilt in the zone cells
//   UInc_c = uI_B - Rc/rhoI_B,  Rc = rhoI_A (uI_A - UInc_A)
//              (the projection correction R[dPhiI/rRhoIf] of state A, kept)
//   Uc  += UInc_A - UInc_B      (zone cells; UInc changes only there)
//   pc  += (pIncUse + Drho gh)_A - (pIncUse + Drho gh)_B   (every cell; p =
//              pIncUse + pc + Drho gh, swensePcEqn.H. With the discrete
//              hydrostatics a change of rhoI changes pIncD below it)
// phiInc, UfInc, phiIRel, phi and rho (from alpha1) are not touched, so the
// total U = UInc + Uc, p = pIncUse + pc + Drho gh, rho and phi are unchanged at the
// reset; only the split moves. The face incident (UfInc, phiInc) stays that
// of state A.
{
    const labelList& zc = mesh.cellZones()[plainBodyZoneID];
    const scalar tNew = runTime.value();

    // Exact wetted cell fraction on the current points
    swenseSurface resetSurf
    (
        mesh,
        mesh.thisDb().lookupObject<waveTheories::externalWaveForcing>
        (
            "externalWaveForcing"
        ),
        tNew
    );

    const scalarField& V = mesh.V();
    scalarField aNew(zc.size());
    forAll(zc, i)
    {
        const label c = zc[i];
        aNew[i] = min(max(resetSurf.wetVolume(c)/V[c], scalar(0)), scalar(1));
    }

    // State A
    const scalarField aA(alphaI.primitiveField(), zc);
    const scalarField rhoIA(rhoI.primitiveField(), zc);
    const vectorField UIncA(incident.UInc().primitiveField(), zc);

    vectorField uIA, uIB;
    scalarField pBA, pBB;
    incident.blendCells(zc, aA, uIA, pBA);
    incident.blendCells(zc, aNew, uIB, pBB);

    // Total-U and total-p check fields (state A)
    const volVectorField UtotA("UtotA", incident.UInc() + Uc);

    tmp<volScalarField> tpIncA;
    if (hydroDiscrete)
    {
        #include "swenseHydroDiscrete.H"
        tpIncA = tmp<volScalarField>(new volScalarField("pIncA", pIncD));
    }
    else
    {
        tpIncA = tmp<volScalarField>(new volScalarField("pIncA", incident.pInc()));
    }
    const volScalarField ptotA("ptotA", tpIncA() + pc + Drho*gh);
    const scalarField pShiftA(ptotA.primitiveField() - pc.primitiveField());

    // Reset alphaI and rebuild rhoI, Drho
    scalarField& aI = alphaI.primitiveFieldRef();
    forAll(zc, i)
    {
        aI[zc[i]] = aNew[i];
    }
    alphaI.correctBoundaryConditions();

    rhoI == swenseRho(alphaI, rho1, rho2);
    Drho == rho - rhoI;

    // Cell incident in the zone (projection correction of state A kept)
    vectorField& UI = incident.UIncRW().primitiveFieldRef();
    scalarField& pI = incident.pIncRW().primitiveFieldRef();
    vectorField& UcI = Uc.primitiveFieldRef();
    const scalarField& rI = rhoI.primitiveField();

    scalar dAlphaMax = 0, dUIMax = 0;
    forAll(zc, i)
    {
        const label c = zc[i];
        const vector Rc = rhoIA[i]*(uIA[i] - UIncA[i]);
        const vector UB = uIB[i] - Rc/rI[c];

        UcI[c] += UIncA[i] - UB;
        UI[c] = UB;
        pI[c] = pBB[i];

        dAlphaMax = max(dAlphaMax, mag(aNew[i] - aA[i]));
        dUIMax = max(dUIMax, mag(UB - UIncA[i]));
    }
    // No boundary update of UInc and Uc: the zone cells touch only body
    // patches, whose values do not depend on the cell values.

    // pc shift by the change of the closure incident pressure
    scalarField dpc(mesh.nCells());
    if (hydroDiscrete)
    {
        #include "swenseHydroDiscrete.H"
        dpc = pShiftA - (pIncD.primitiveField() + Drho.primitiveField()*gh.primitiveField());
    }
    else
    {
        dpc = pShiftA - (incident.pInc().primitiveField() + Drho.primitiveField()*gh.primitiveField());
    }
    pc.primitiveFieldRef() += dpc;

    // Boundary values: shifted with their face cells on the patches that do
    // not fix the value (fixedFluxPressure cannot be evaluated outside the
    // pc equation); fixed-value patches keep their value.
    forAll(pc.boundaryField(), patchi)
    {
        fvPatchScalarField& pb = pc.boundaryFieldRef()[patchi];
        if (pb.coupled())
        {
            pb.evaluate();
        }
        else if (!pb.fixesValue() && pb.size())
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(pb, i)
            {
                pb[i] += dpc[fc[i]];
            }
        }
    }

    // Check: totals before and after the reset
    const volScalarField& pIncB = hydroDiscrete ? pIncD : incident.pInc();
    const scalar dU = gMax(mag(incident.UInc().primitiveField() + Uc.primitiveField() - UtotA.primitiveField()));
    const scalar dp = gMax(mag(pIncB.primitiveField() + pc.primitiveField() + Drho.primitiveField()*gh.primitiveField() - ptotA.primitiveField()));
    const scalar pScale = gMax(mag(ptotA.primitiveField()));

    reduce(dAlphaMax, maxOp<scalar>());
    reduce(dUIMax, maxOp<scalar>());

    Info<< "SWENSE resetInZone: cells " << returnReduce(zc.size(), sumOp<label>())
        << " max|dAlphaI| " << dAlphaMax
        << " max|dUInc| " << dUIMax
        << " check max|dU_tot| " << dU
        << " max|dp_tot| " << dp << " (max|p_tot| " << pScale << ")"
        << endl;
}
