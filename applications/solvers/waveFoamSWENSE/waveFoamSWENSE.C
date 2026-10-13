/*---------------------------------------------------------------------------*\
Application
    waveFoamSWENSE

Description
    Two-phase (VoF) incompressible wave solver with the SWENSE decomposition:

        U = Uc + UInc,   p = pc + pStar,   pStar = (rho/rho1) pInc

    The incident fields come from swenseIncident (Grid2Grid water, Air2Air
    air). The alpha transport is the OpenFOAM VoF MULES step. The rest of the
    solver is shared with waveIsoFoamSWENSE (src/swenseWaveTank/solver).

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
#include "externalWaveForcing.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "swenseSetup.H"
    #include "initContinuityErrs.H"
    #include "createSolverFields.H"

    Info<< "\nStarting time loop\n" << endl;

    #include "swenseTimeLoop.H"

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
