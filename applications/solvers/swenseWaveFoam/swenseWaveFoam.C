/*---------------------------------------------------------------------------*\
Application
    swenseWaveFoam

Description
    Two-phase (VoF) incompressible solver with the SWENSE decomposition of the
    flow into incident (wave) and complementary parts:

        U = Uc + UInc,   p = pc + pStar,   pStar = (rho/rho1) pInc

    The incident fields come from swenseIncident (Grid2Grid water, Air2Air air),
    the complementary fields are solved in a PIMPLE loop, and Uc is relaxed to
    zero in the relaxation zone with the farfield weight of swenseRelaxation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "subCycle.H"
#include "adjustPhi.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "swenseIncident.H"
#include "swenseRelaxation.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "alphaControls.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Grid2Grid record update for this time step (as waveFoam does)
        externalWave->step();

        // Incident fields, incident flux and total flux
        incident.update(runTime.value());
        phiInc = incident.UfInc() & mesh.Sf();
        phi = phiComp + phiInc;

        while (pimple.loop())
        {
            // Alpha transport (OpenFOAM VoF MULES), then relaxation toward the
            // incident water fraction, then density update
            #include "alphaEqnSubCycle.H"

            alpha1 = (1.0 - w)*alpha1 + w*incident.alphaInc();

            {
                const volScalarField limitedAlpha1
                (
                    "limitedAlpha1",
                    min(max(alpha1, scalar(0)), scalar(1))
                );
                rho == limitedAlpha1*rho1 + (1.0 - limitedAlpha1)*rho2;
                mixture.correct();
            }

            // Incident pressure scaled to the mixture density
            pStar == rho/rho1*incident.pInc();

            #include "UcEqn.H"

            while (pimple.correct())
            {
                #include "pcEqn.H"
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

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
