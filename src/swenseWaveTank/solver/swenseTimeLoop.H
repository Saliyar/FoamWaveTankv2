// Shared SWENSE time loop. The solver provides alphaStep.H, which advances
// alpha1 and alphaI with the same function call in every outer corrector.

// Per-bin budget logging switch (default off)
const bool budgetLog =
    runTime.controlDict().lookupOrDefault<bool>("budgetLog", false);

while (runTime.run())
{
    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Grid2Grid record update for this time step (as waveFoam does)
    externalWave->step();

    // Pressure-velocity PIMPLE loop. The mesh moves at the first outer
    // corrector, and at every corrector when moveMeshOuterCorrectors is on
    // (as waveFoam). The incident fields follow the mesh (document 16.3.5).
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            mesh.update();

            if (mesh.changing())
            {
                gh = (g & mesh.C()) - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;

                if (mesh.moving())
                {
                    // Swept mesh flux and the Grid2Grid/Air2Air points
                    meshPhi = mesh.phi();
                    incident.movePoints();

                    // Relaxation zones follow the mesh (spec item 8, option a)
                    relax.update();
                }
            }

            #include "swenseIncidentStep.H"

            // Moving mesh: complementary flux from the face velocity, made
            // divergence-free with the incident flux (waveFoam correctPhi)
            if (mesh.changing() && correctPhi)
            {
                #include "swenseCorrectPhi.H"
            }
        }

        // Alpha1 with phi and alphaI with phiIRel, same function (16.3.2)
        #include "alphaStep.H"

        // Relaxation toward alphaH, densities, incident pressure
        #include "swenseAlphaRelax.H"

        // Momentum (Uc) and pressure (pc), SWENSE form
        #include "swenseUcEqn.H"

        while (pimple.correct())
        {
            #include "swensePcEqn.H"
        }

        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }

    // Consistency check of the discrete invariance (document Prop. 16.1)
    Info<< "SWENSE check: max|alpha1 - alphaI| = "
        << max(mag(alpha1 - alphaI)).value()
        << "  max|Uc| = " << max(mag(Uc)).value()
        << "  max|div(phiI)| = "
        << max(mag(fvc::div(phiInc))).value()
        << "  alpha1 [" << min(alpha1).value() << ", "
        << max(alpha1).value() << "]"
        << "  alphaI [" << min(alphaI).value() << ", "
        << max(alphaI).value() << "]" << endl;

    // Logging only: per-bin alphaI budget (no effect on the solution).
    // Off unless budgetLog true is set in controlDict.
    if (budgetLog)
    {
        #include "swenseBudgetLog.H"
    }

    // Relaxation of the complementary velocity in the zone
    Uc *= (1.0 - w);
    Uc.correctBoundaryConditions();

    U == Uc + incident.UInc();
    p == pc + pStar;
    p_rgh == p - rho*gh;

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}
