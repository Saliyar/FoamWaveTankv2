// Logging only (controlDict ucCellLog, default off). No equation, field or
// tolerance is changed. Included at the end of the last outer corrector, so
// incTerm, lam, lamf, faceLim of swenseUcEqn.H are in scope and Uc, pc are
// the final fields of the step (before the zone relaxation Uc *= 1 - w).
//
// For each cell of the box ucCellLogBox (xMin xMax zMin zMax, cell centres
// at the first call) one line per step, per unit volume (N/m3), x and z:
//
//   ddt   fvc::ddt(rho, Uc)
//   conv  fvc::div(rhoPhi, Uc)
//   visc  divDevRhoReff(rho, Uc) & Uc
//   inc   incident mismatch brackets of incTerm (ddt, convection, viscous)
//   clo   cell closure (1 - lam)[ddt(rhoI UI) + div(rhoPhiI UI) - V(alphaI; UI)]
//   face  -reconstruct((1 - lamf) G_f)
//   gpc   -reconstruct(snGrad(pc) |Sf|)
//   grav  -reconstruct(gh_f snGrad(Drho) |Sf|)
//   R     ddt + conv + visc + inc + clo - face - gpc - grav
//
// and for each atmosphere face above the box: pc_b, Uc_b, phiComp/|Sf|,
// snGrad(pc). Output: postProcessing/ucCellLog/<startTime>/ucCellLog.dat
{
    static bool ucInit = true;
    static labelList ucCells;
    static labelList ucFaces;
    static label ucPatch = -1;
    static autoPtr<OFstream> ucOs;

    if (ucInit)
    {
        ucInit = false;

        const scalarList box
        (
            runTime.controlDict().lookupOrDefault<scalarList>
            (
                "ucCellLogBox",
                scalarList({0.9, 1.4, 0.47, 1.0})
            )
        );

        DynamicList<label> cs;
        forAll(mesh.C(), c)
        {
            const vector& C = mesh.C()[c];
            if
            (
                C.x() > box[0] && C.x() < box[1]
             && C.z() > box[2] && C.z() < box[3]
            )
            {
                cs.append(c);
            }
        }
        ucCells.transfer(cs);

        ucPatch = mesh.boundaryMesh().findPatchID("atmosphere");
        DynamicList<label> fs;
        if (ucPatch >= 0)
        {
            const vectorField& Cf = mesh.Cf().boundaryField()[ucPatch];
            forAll(Cf, i)
            {
                if (Cf[i].x() > box[0] && Cf[i].x() < box[1])
                {
                    fs.append(i);
                }
            }
        }
        ucFaces.transfer(fs);

        if (Pstream::master())
        {
            const fileName dir
            (
                runTime.path()/"postProcessing"/"ucCellLog"
               /runTime.timeName()
            );
            mkDir(dir);
            ucOs.reset(new OFstream(dir/"ucCellLog.dat"));
            ucOs()
                << "# C t dt cell x z alpha1 alphaI rho rhoI lam w Ucx Ucz pc"
                << " ddt conv visc inc clo face gpc grav R (x z each)" << nl
                << "# B t dt face x pc_b Ucbx Ucbz phiComp/magSf snGradPc"
                << nl;
        }

        Info<< "ucCellLog: " << ucCells.size() << " cells, "
            << ucFaces.size() << " atmosphere faces" << endl;
    }

    const volVectorField& UIc = incident.UInc();
    const word convC("div(rhoPhi,Uc)");

    const volVectorField tDdt(fvc::ddt(rho, Uc));
    const volVectorField tConv(fvc::div(rhoPhi, Uc, convC));
    const fvVectorMatrix McC(turbulence->divDevRhoReff(rho, Uc));
    const volVectorField tVisc(McC & Uc);

    const volVectorField tClo
    (
        (1.0 - lam)
       *(
            fvc::ddt(rhoI, UIc) + fvc::div(rhoPhiI, UIc, convC)
          - viscUI(muOf(alphaI)())
        )
    );
    const volVectorField tInc(incTerm - tClo);

    const volVectorField tFace(-fvc::reconstruct(faceLim));
    const volVectorField tGpc
    (
        -fvc::reconstruct(fvc::snGrad(pc)*mesh.magSf())
    );
    const volVectorField tGrav
    (
        -fvc::reconstruct(ghf*fvc::snGrad(Drho)*mesh.magSf())
    );

    const volVectorField tR
    (
        tDdt + tConv + tVisc + tInc + tClo - tFace - tGpc - tGrav
    );

    if (Pstream::master())
    {
        OFstream& os = ucOs();
        os.precision(7);
        const scalar tNow = runTime.value();
        const scalar dtNow = runTime.deltaTValue();

        const volVectorField* terms[9] =
        {
            &tDdt, &tConv, &tVisc, &tInc, &tClo, &tFace, &tGpc, &tGrav, &tR
        };

        for (const label c : ucCells)
        {
            os  << "C " << tNow << ' ' << dtNow << ' ' << c << ' '
                << mesh.C()[c].x() << ' ' << mesh.C()[c].z() << ' '
                << alpha1[c] << ' ' << alphaI[c] << ' '
                << rho[c] << ' ' << rhoI[c] << ' '
                << lam[c] << ' ' << w[c] << ' '
                << Uc[c].x() << ' ' << Uc[c].z() << ' ' << pc[c];

            for (const volVectorField* T : terms)
            {
                os << ' ' << (*T)[c].x() << ' ' << (*T)[c].z();
            }
            os << nl;
        }

        if (ucPatch >= 0)
        {
            const fvPatchScalarField& pcb = pc.boundaryField()[ucPatch];
            const fvPatchVectorField& Ucb = Uc.boundaryField()[ucPatch];
            const scalarField snPc(pcb.snGrad());
            const scalarField& phb = phiComp.boundaryField()[ucPatch];
            const scalarField& mSf = mesh.magSf().boundaryField()[ucPatch];
            const vectorField& Cf = mesh.Cf().boundaryField()[ucPatch];

            for (const label i : ucFaces)
            {
                os  << "B " << tNow << ' ' << dtNow << ' ' << i << ' '
                    << Cf[i].x() << ' ' << pcb[i] << ' '
                    << Ucb[i].x() << ' ' << Ucb[i].z() << ' '
                    << phb[i]/mSf[i] << ' ' << snPc[i] << nl;
            }
        }
        os.flush();
    }
}
