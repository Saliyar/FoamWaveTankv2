// Zone source split (waveProperties logZoneSource on, plainBodyZone set).
// Logging only, nothing here enters a solve. Per step:
//   max water |Uc| outside the zone (alpha1 > 0.5, zoneLayer > 0), its cell,
//   and over that cell and its face neighbours the max of
//     cell : |(1-lam) rI_c|, rI_c = ddt(rhoI UI) + div(rhoPhiI UI) - V[muI; UI]
//     face : |R[(1-lam_f) G_f]|  (faceLim, the face part of rI)
//     rest : |incTerm - cell|    (the SWENSE incident term without closure)
//     pres : |R[-ghf snGrad(Drho) - snGrad(pc)]|  (response, for reference)
//   and the same terms integrated (|.| V) over the edge band (zone cells
//   with an outside neighbour, plus layer 1) and over the rest of the zone.
// First step with the max above zoneSourceThreshold: the faces of that cell
// are listed and the fields zsCell, zsFace, zsRest, Uc are written.
{
    const volVectorField zsCell
    (
        "zsCell",
        (1.0 - lam)
       *(
            fvc::ddt(rhoI, UI) + fvc::div(rhoPhiI, UI, convScheme)
          - viscUI(muOf(alphaI)())
        )
    );
    const volVectorField zsFace("zsFace", fvc::reconstruct(faceLim));
    const volVectorField zsRest("zsRest", incTerm - zsCell);
    const volVectorField zsPres
    (
        "zsPres",
        fvc::reconstruct((-ghf*fvc::snGrad(Drho) - fvc::snGrad(pc))*mesh.magSf())
    );

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V();

    scalar uMax = -1;
    label cMax = -1;

    forAll(Uc, c)
    {
        if (zoneLayer[c] > 0 && alpha1[c] > 0.5)
        {
            const scalar u = mag(Uc[c]);
            if (u > uMax)
            {
                uMax = u;
                cMax = c;
            }
        }
    }

    // Edge band and zone interior integrals
    scalar bC = 0, bF = 0, bR = 0, zC = 0, zF = 0, zR = 0;

    forAll(Uc, c)
    {
        bool band = (zoneLayer[c] == 1);

        if (zoneLayer[c] == 0)
        {
            for (const label n : mesh.cellCells()[c])
            {
                if (zoneLayer[n] > 0)
                {
                    band = true;
                }
            }
        }

        if (band)
        {
            bC += mag(zsCell[c])*V[c];
            bF += mag(zsFace[c])*V[c];
            bR += mag(zsRest[c])*V[c];
        }
        else if (zoneLayer[c] == 0)
        {
            zC += mag(zsCell[c])*V[c];
            zF += mag(zsFace[c])*V[c];
            zR += mag(zsRest[c])*V[c];
        }
    }

    if (cMax >= 0)
    {
        scalar sC = mag(zsCell[cMax]), sF = mag(zsFace[cMax]);
        scalar sR = mag(zsRest[cMax]), sP = mag(zsPres[cMax]);

        for (const label n : mesh.cellCells()[cMax])
        {
            sC = max(sC, mag(zsCell[n]));
            sF = max(sF, mag(zsFace[n]));
            sR = max(sR, mag(zsRest[n]));
            sP = max(sP, mag(zsPres[n]));
        }

        Info<< "SWENSE zoneSource: t " << runTime.value()
            << " UcOut " << uMax << " cell " << cMax
            << " x " << C[cMax].x() << " z " << C[cMax].z()
            << " layer " << zoneLayer[cMax]
            << " stencil cell " << sC << " face " << sF
            << " rest " << sR << " pres " << sP
            << " | band cell " << bC << " face " << bF << " rest " << bR
            << " | zone cell " << zC << " face " << zF << " rest " << zR
            << endl;

        if (!zoneSourceCrossed && uMax > zoneSourceThreshold)
        {
            zoneSourceCrossed = true;

            Info<< "SWENSE zoneSource CROSS t " << runTime.value()
                << " cell " << cMax << " Uc " << Uc[cMax]
                << " alpha1 " << alpha1[cMax] << " alphaI " << alphaI[cMax]
                << nl;

            for (const label f : mesh.cells()[cMax])
            {
                if (!mesh.isInternalFace(f))
                {
                    continue;
                }

                const label o = mesh.owner()[f];
                const label nb = (o == cMax ? mesh.neighbour()[f] : o);

                Info<< "    face " << f << " Cf " << mesh.Cf()[f]
                    << " lamf " << lamf[f] << " faceLim " << faceLim[f]
                    << " other cell " << nb << " layer " << zoneLayer[nb]
                    << " alpha1 " << alpha1[nb] << " |cell| "
                    << mag(zsCell[nb]) << " |Uc| " << mag(Uc[nb]) << nl;
            }
            Info<< endl;

            zsCell.write();
            zsFace.write();
            zsRest.write();
            Uc.write();
            alpha1.write();
        }
    }
}
