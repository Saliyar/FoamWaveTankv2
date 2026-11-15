// SWENSE time loop, momentum blend (momentum_swense_Rev1, Part II).
// Order per time step:
//   wave models at the new time, relaxation target alphaH
//   1a  alphaI with the incident flux and incident velocity of level n
//   2   alpha1 with the total flux and total velocity of level n
//       relaxation of both fractions, densities, Drho
//   1b-1d  incident blend, raw face flux, projection (new level)
//   3   momentum (Uc) and pressure (pc)
//   4   total fields
// The solver provides alphaStep.H (steps 1a and 2).

// Per-bin budget logging switch (default off)
const bool budgetLog =
    runTime.controlDict().lookupOrDefault<bool>("budgetLog", false);

while (runTime.run())
{
    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Grid2Grid record update for this time step (as waveFoam does), then
    // the water and air wave models and the relaxation target
    externalWave->step();

    // Static mesh: the incident level is set once per time step, before the
    // PIMPLE loop (unchanged order)
    if (!mesh.dynamic())
    {
        incident.evaluate(runTime.value());

        alphaH.primitiveFieldRef() =
            incident.blendAlpha(runTime.value()).primitiveField();

        // Steps 1a and 2: both fractions by the same routine, with the fluxes
        // and velocities of the previous level
        #include "alphaStep.H"

        // Relaxation toward alphaH, densities, Drho
        #include "swenseAlphaRelax.H"

        // Steps 1b to 1d: incident field of the new level
        #include "swenseIncidentStep.H"
    }

    // Step 3: momentum (Uc) and pressure (pc)
    while (pimple.loop())
    {
        // Moving mesh (dyn port): mesh update at the first outer corrector,
        // and at each one when moveMeshOuterCorrectors is on (as waveIsoFoam).
        // The incident level follows the mesh: geometry, then flux, then the
        // incident steps of the new level.
        if (mesh.dynamic() && (pimple.firstIter() || moveMeshOuterCorrectors))
        {
            mesh.update();

            if (mesh.changing())
            {
                gh = (g & mesh.C()) - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;

                if (mesh.moving())
                {
                    // Swept mesh flux, incident points and relaxation weight
                    meshPhi = mesh.phi();
                    incident.movePoints();
                    relax.update();
                }
            }

            incident.evaluate(runTime.value());

            alphaH.primitiveFieldRef() =
                incident.blendAlpha(runTime.value()).primitiveField();

            #include "alphaStep.H"
            #include "swenseAlphaRelax.H"
            #include "swenseIncidentStep.H"

            // Complementary flux from the face velocity on the new mesh
            if (mesh.changing() && correctPhi)
            {
                #include "swenseCorrectPhi.H"
            }
        }

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

    // Logging only: theorem residuals (switch theoremLog)
    if (theoremOn)
    {
        #include "swenseTheorem.H"
    }

    // Logging only: term budget of the Uc equation (switch budgetTerms)
    if (budgetOn)
    {
        #include "swenseBudget.H"
    }

    // Discrete invariance (Theorem 6) and incident divergence (G1, G2)
    Info<< "SWENSE check: max|alpha1 - alphaI| = "
        << max(mag(alpha1 - alphaI)).value()
        << "  max|Uc| = " << max(mag(Uc)).value()
        << "  max|pc| = " << max(mag(pc)).value()
        << "  max|phiComp| = " << max(mag(phiComp)).value()
        << "  max|div(phiI)| = "
        << max(mag(fvc::div(phiInc))).value()
        << "  alpha1 [" << min(alpha1).value() << ", "
        << max(alpha1).value() << "]"
        << "  alphaI [" << min(alphaI).value() << ", "
        << max(alphaI).value() << "]" << endl;

    // Level of alphaI against the wave-model elevation (G3), log only
    #include "swenseLevelLog.H"

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
    p == incident.pInc() + pc + Drho*gh;
    p_rgh == p - rho*gh;

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}
