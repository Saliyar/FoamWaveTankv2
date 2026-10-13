/*---------------------------------------------------------------------------*\
Application
    waveIsoFoamSWENSE

Description
    Two-phase (VoF) incompressible wave solver with the SWENSE decomposition,
    using the isoAdvector interface advection for alpha (as waveIsoFoam).
    Everything except the alpha step is shared with waveFoamSWENSE
    (src/swenseWaveTank/solver). The alpha advection uses the total SWENSE flux
    phi = phiComp + phiInc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
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
