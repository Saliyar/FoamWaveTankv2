// Logging only (no effect on the solution). Budget of alphaI per 0.1 m bin of
// the free zone (|x| < 1.3 m), accumulated every time step after relaxation.
//   alphaPhiI_f = (rhoPhiI_f - rho2 phiInc_f)/(rho1 - rho2)   (from the advect call)
//   dV_b = sum (alphaI - alphaI.oldTime()) V                  (volume change)
//   F_b  = dt * net alphaPhiI into the bin through its x faces
//   S_b  = dt * sum alphaI div(phiInc) V, split band / bulk
//   R_b  = dV_b - F_b - S_band_b - S_bulk_b                    (closure residual)
// Reports at t = 10, 20, 30, 40 s to budget_bins_<t>.dat (master writes).
{
    const scalar binW = 0.1;
    const label nBins = 26;
    const scalar xFree = 1.3;
    const scalar dtB = runTime.deltaTValue();
    static scalarList cumF(nBins, 0.0), cumSb(nBins, 0.0), cumSk(nBins, 0.0), cumR(nBins, 0.0);
    static scalarList cumV(nBins, 0.0);
    static scalar dyB = -1;

    auto binOf = [&](const scalar x) -> label
    {
        if (mag(x) >= xFree) return -1;
        label b = label(floor((x + xFree)/binW));
        return min(max(b, label(0)), nBins - 1);
    };

    if (dyB < 0)
    {
        // Depth (y-width) of the mesh: total volume over total x-z area of the free zone
        scalar vol = 0, area = 0;
        forAll(mesh.C(), c)
        {
            if (binOf(mesh.C()[c].x()) < 0) continue;
            scalar x0 = GREAT, x1 = -GREAT, z0 = GREAT, z1 = -GREAT;
            for (const label p : mesh.cellPoints(c))
            {
                x0 = min(x0, mesh.points()[p].x()); x1 = max(x1, mesh.points()[p].x());
                z0 = min(z0, mesh.points()[p].z()); z1 = max(z1, mesh.points()[p].z());
            }
            vol += mesh.V()[c];
            area += (x1 - x0)*(z1 - z0);
        }
        dyB = returnReduce(vol, sumOp<scalar>())/max(returnReduce(area, sumOp<scalar>()), SMALL);
    }

    scalarList dVb(nBins, 0.0), Fb(nBins, 0.0), Sbb(nBins, 0.0), Skb(nBins, 0.0), Vnow(nBins, 0.0);

    const volScalarField divPhiIb(fvc::div(phiInc));
    const scalar r12 = rho1.value() - rho2.value();

    forAll(alphaI, c)
    {
        const label b = binOf(mesh.C()[c].x());
        if (b < 0) continue;
        const scalar a = alphaI[c];
        Vnow[b] += a*mesh.V()[c];
        dVb[b] += (a - alphaI.oldTime()[c])*mesh.V()[c];
        const scalar src = a*divPhiIb[c]*mesh.V()[c]*dtB;
        if (a > 1e-6 && a < 1.0 - 1e-6) Sbb[b] += src; else Skb[b] += src;
    }

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    forAll(nei, f)
    {
        const label bo = binOf(mesh.C()[own[f]].x());
        const label bn = binOf(mesh.C()[nei[f]].x());
        if (bo == bn) continue;
        const scalar aphi = (rhoPhiI[f] - rho2.value()*phiInc[f])/r12;
        if (bo >= 0) Fb[bo] -= aphi*dtB;
        if (bn >= 0) Fb[bn] += aphi*dtB;
    }

    // Coupled (processor) faces: the face is not in the internal loop above.
    // Each rank adds only its local side, with the face normal pointing out of
    // the local cell, so the two ranks together give the same flux as an
    // internal face.
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& pp = mesh.boundary()[patchi];

        if (!pp.coupled()) continue;

        const labelUList& fc = pp.faceCells();
        const scalarField aphiP
        (
            (rhoPhiI.boundaryField()[patchi] - rho2.value()*phiInc.boundaryField()[patchi])/r12
        );

        forAll(fc, i)
        {
            const label b = binOf(mesh.C()[fc[i]].x());
            if (b >= 0) Fb[b] -= aphiP[i]*dtB;
        }
    }

    for (label b = 0; b < nBins; b++)
    {
        dVb[b] = returnReduce(dVb[b], sumOp<scalar>());
        Fb[b] = returnReduce(Fb[b], sumOp<scalar>());
        Sbb[b] = returnReduce(Sbb[b], sumOp<scalar>());
        Skb[b] = returnReduce(Skb[b], sumOp<scalar>());
        Vnow[b] = returnReduce(Vnow[b], sumOp<scalar>());
        cumV[b] += dVb[b];
        cumF[b] += Fb[b];
        cumSb[b] += Sbb[b];
        cumSk[b] += Skb[b];
        cumR[b] = cumV[b] - cumF[b] - cumSb[b] - cumSk[b];
    }

    // Per-step trace (master): totals over the free zone, for early checks
    if (Pstream::master())
    {
        static OFstream trace(fileName("budget_trace.dat"));
        static bool first = true;
        if (first)
        {
            trace << "# t dV_total F_total S_band_total S_bulk_total residual_total" << nl;
            first = false;
        }
        trace.precision(10);
        trace << runTime.value() << " " << sum(cumV) << " " << sum(cumF) << " "
              << sum(cumSb) << " " << sum(cumSk) << " " << sum(cumR) << nl;
    }

    const scalar tNow = runTime.value();
    for (const scalar tRep : {10.0, 20.0, 30.0, 40.0})
    {
        if (mag(tNow - tRep) < 0.5*dtB + 1e-9)
        {
            // HOS mean level per bin: eta at the cell-centre x values (all ranks reduce)
            scalarList eS(nBins, 0.0), nE(nBins, 0.0);
            forAll(mesh.C(), c)
            {
                const label b = binOf(mesh.C()[c].x());
                if (b < 0) continue;
                eS[b] += externalWave->eta(point(mesh.C()[c].x(), 0, 0), tNow);
                nE[b] += 1.0;
            }
            for (label b = 0; b < nBins; b++)
            {
                eS[b] = returnReduce(eS[b], sumOp<scalar>());
                nE[b] = returnReduce(nE[b], sumOp<scalar>());
            }
            if (Pstream::master())
            {
                OFstream os(fileName("budget_bins_" + name(label(tRep)) + ".dat"));
                os.precision(10);
                os << "# t " << tNow << " dy " << dyB << " binWidth " << binW << nl;
                os << "# bin xc level_m Emean_m F_m dVcum_m S_band_m S_bulk_m residual_m\n";
                for (label b = 0; b < nBins; b++)
                {
                    const scalar xc = -xFree + (b + 0.5)*binW;
                    const scalar A = binW*dyB;
                    os << xc << " " << (Vnow[b]/A - 5.0) << " " << eS[b]/max(nE[b], 1.0)
                       << " " << cumF[b]/A << " " << cumV[b]/A << " " << cumSb[b]/A
                       << " " << cumSk[b]/A << " " << cumR[b]/A << nl;
                }
            }
        }
    }
}
