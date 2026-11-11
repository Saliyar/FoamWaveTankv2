// Logging only (switch theoremLog in controlDict, default off). Evaluated after
// the pc solve, with fvc on the final fields and the solver's fluxes.
//
// R_tot: the plain two-phase momentum residual of waveIsoFoam on the TOTAL
//   fields, per unit mass:
//     ddt(rho,U) + div(rhoPhi,U) + divDevRhoReff(rho,U)
//       - reconstruct((-snGrad(p_rgh) - ghf snGrad(rho)) magSf)
//   with p_rgh the total p_rgh (p_rgh = p - rho gh, as the solver stores it).
//
// r_I: the same operator on the INCIDENT fields, per unit mass (rhoI):
//     ddt(rhoI,UInc) + div(rhoPhiI,UInc) + [viscous difference of the incident
//     fields as in swenseUcEqn.H, mu(alphaI)]
//       - reconstruct((-snGrad(pIrgh) - ghf snGrad(rhoI)) magSf)
//   with pIrgh = pInc - rhoI gh: the incident pressure of the relaxation zone
//   pInc (incident.pInc(), cell values, Air2Air and wave model) converted to
//   the rgh form with the incident density. snGrad is OpenFOAM's corrected
//   face gradient of the cell field.
//
// Set S: air cells (alpha.water < 0.5) with 0 <= z - eta <= 0.2, eta of the
// nearest surface row (Grid2Grid). Split: free region |x| < 1.35 m and zone.
{
    static bool thmInit = true;
    static labelList thmRow;

    if (thmInit)
    {
        thmInit = false;

        const pointField& rowsT = incident.etaPts();
        const vectorField& CcT = mesh.C().primitiveField();

        thmRow.setSize(mesh.nCells(), 0);

        forAll(CcT, c)
        {
            scalar d2Min = GREAT;

            forAll(rowsT, r)
            {
                const scalar d2 =
                    sqr(rowsT[r].x() - CcT[c].x()) + sqr(rowsT[r].y() - CcT[c].y());

                if (d2 < d2Min)
                {
                    d2Min = d2;
                    thmRow[c] = r;
                }
            }
        }
    }

    const scalar tT = runTime.value();
    const vectorField& CT = mesh.C().primitiveField();
    const scalarField& a1T = alpha1.primitiveField();
    const scalarField& aIT = alphaI.primitiveField();
    const scalarField& rhoT = rho.primitiveField();
    const scalarField& rhoIT = rhoI.primitiveField();

    // R_tot (total fields, per unit mass later)
    const volVectorField TT1("thmT1", fvc::ddt(rho, U));
    const volVectorField TT2("thmT2", fvc::div(rhoPhi, U));
    const volVectorField TTv("thmTv", turbulence->divDevRhoReff(rho, U) & U);
    const volVectorField TTp
    (
        "thmTp",
        fvc::reconstruct((-fvc::snGrad(p_rgh) - ghf*fvc::snGrad(rho))*mesh.magSf())
    );
    const volVectorField Rtot("thmRtot", TT1 + TT2 + TTv - TTp);

    // r_I (incident fields)
    const volVectorField& UIT = incident.UInc();
    const volScalarField pIrgh("thmPIrgh", incident.pInc() - rhoI*gh);

    auto muOfT = [&](const volScalarField& s) -> tmp<volScalarField>
    {
        const volScalarField lim(min(max(s, scalar(0)), scalar(1)));

        return
            lim*rho1*mixture.nuModel1().nu()
          + (scalar(1) - lim)*rho2*mixture.nuModel2().nu();
    };

    const volTensorField gradUIT("thmGradUI", fvc::grad(UIT));

    auto viscT = [&](const volScalarField& mu) -> tmp<volVectorField>
    {
        return
            fvc::laplacian(mu, UIT, "laplacian(mu,UInc)")
          + fvc::div(mu*dev2(T(gradUIT)), "div(mu*dev2(T(grad(UInc))))");
    };

    const volVectorField rI1("thmI1", fvc::ddt(rhoI, UIT));
    const volVectorField rI2("thmI2", fvc::div(rhoPhiI, UIT, "div(rhoPhi,Uc)"));
    const volVectorField rI3("thmI3", -(viscT(muOfT(alpha1)) - viscT(muOfT(alphaI))));
    const volVectorField rIp
    (
        "thmIp",
        fvc::reconstruct((-fvc::snGrad(pIrgh) - ghf*fvc::snGrad(rhoI))*mesh.magSf())
    );
    const volVectorField rI("thmrI", rI1 + rI2 + rI3 - rIp);

    // Stats over air cells in 0..0.2 m above the surface
    scalar mxR[2] = {0.0, 0.0}, mxr[2] = {0.0, 0.0}, mxD[2] = {0.0, 0.0};
    scalar ssR[2] = {0.0, 0.0}, ssr[2] = {0.0, 0.0}, ssD[2] = {0.0, 0.0};
    label nT[2] = {0, 0};
    scalar maxUcT = 0.0;
    label cellT = -1;

    const tmp<volScalarField> tUcMag = mag(Uc);
    const scalarField& UcMag = tUcMag().primitiveField();

    forAll(CT, c)
    {
        const scalar zRel = CT[c].z() - incident.etaRow(thmRow[c], tT);

        if (a1T[c] < 0.5 && zRel >= 0.0 && zRel <= 0.2)
        {
            const label zone = (mag(CT[c].x()) < 1.35) ? 0 : 1;
            const vector Rm = Rtot.primitiveField()[c]/rhoT[c];
            const vector rm = rI.primitiveField()[c]/rhoIT[c];
            const vector Dm = Rm - rm;

            mxR[zone] = max(mxR[zone], mag(Rm));
            mxr[zone] = max(mxr[zone], mag(rm));
            mxD[zone] = max(mxD[zone], mag(Dm));
            ssR[zone] += magSqr(Rm);
            ssr[zone] += magSqr(rm);
            ssD[zone] += magSqr(Dm);
            nT[zone]++;

            if (UcMag[c] > maxUcT)
            {
                maxUcT = UcMag[c];
                cellT = c;
            }
        }
    }

    for (label z = 0; z < 2; ++z)
    {
        reduce(nT[z], sumOp<label>());
        reduce(mxR[z], maxOp<scalar>());
        reduce(mxr[z], maxOp<scalar>());
        reduce(mxD[z], maxOp<scalar>());
        reduce(ssR[z], sumOp<scalar>());
        reduce(ssr[z], sumOp<scalar>());
        reduce(ssD[z], sumOp<scalar>());
    }

    const char* zn[2] = {"free", "zone"};

    for (label z = 0; z < 2; ++z)
    {
        if (nT[z] > 0)
        {
            Info<< "THM_STATS t " << tT << " " << zn[z] << " n " << nT[z]
                << " Rtot_max " << mxR[z]
                << " Rtot_rms " << Foam::sqrt(ssR[z]/nT[z])
                << " rI_max " << mxr[z]
                << " rI_rms " << Foam::sqrt(ssr[z]/nT[z])
                << " diff_max " << mxD[z]
                << " diff_rms " << Foam::sqrt(ssD[z]/nT[z]) << endl;
        }
    }

    // Cell of max |Uc| in the set (one processor only in serial runs)
    if (cellT >= 0)
    {
        const vector Rm = Rtot.primitiveField()[cellT]/rhoT[cellT];
        const vector rm = rI.primitiveField()[cellT]/rhoIT[cellT];

        Info<< "THM_CELL t " << tT << " cell " << cellT
            << " x " << CT[cellT].x() << " z " << CT[cellT].z()
            << " alpha1 " << a1T[cellT] << " alphaI " << aIT[cellT]
            << " Uc " << UcMag[cellT]
            << " Rtot " << Rm << " rI " << rm << " diff " << (Rm - rm) << endl;
    }
}
