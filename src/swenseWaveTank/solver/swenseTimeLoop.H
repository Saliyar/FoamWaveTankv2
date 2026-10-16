// Shared SWENSE time loop. The solver provides alphaStep.H, which advances
// alpha1 and alphaI with the same function call in every outer corrector.

while (runTime.run())
{
    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Grid2Grid record update for this time step (as waveFoam does)
    externalWave->step();

    // Incident fields are blended with the advected alphaI, which is the
    // value at the start of this step (document Rev1, 16.3.5). The HOS target
    // alphaH is the relaxation target of both volume fractions (16.3.3).
    incident.update(runTime.value(), alphaI);
    alphaH.primitiveFieldRef() =
        incident.blendAlpha(runTime.value()).primitiveField();

    // Incident flux of alphaI and total flux (SWENSE)
    phiInc = incident.UfInc() & mesh.Sf();
    phi = phiComp + phiInc;

    while (pimple.loop())
    {
        // Alpha1 with phi and alphaI with phiInc, same function (16.3.2)
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

    // Logging only: per-bin alphaI budget (no effect on the solution)
    #include "swenseBudgetLog.H"

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
