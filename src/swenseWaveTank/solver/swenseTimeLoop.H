// Shared SWENSE time loop. The solver provides alphaStep.H (its alpha
// transport) in its own directory.

while (runTime.run())
{
    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Grid2Grid record update for this time step (as waveFoam does)
    externalWave->step();

    // Incident fields, incident flux and total flux (SWENSE)
    incident.update(runTime.value());
    phiInc = incident.UfInc() & mesh.Sf();
    phi = phiComp + phiInc;

    while (pimple.loop())
    {
        // Alpha transport: solver-specific (VoF or isoAdvector)
        #include "alphaStep.H"

        // Relaxation toward the incident water fraction, density update
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
