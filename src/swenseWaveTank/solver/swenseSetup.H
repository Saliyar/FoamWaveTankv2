// Shared SWENSE set-up: fields, mixture, incident and relaxation objects.
// Included by both waveFoamSWENSE and waveIsoFoamSWENSE. The alpha-specific
// objects are created in the solver's createSolverFields.H.

// Gravity (OpenFOAM v2412 standard), waveProperties and the Grid2Grid
// external wave forcing (waves2Foam), as waveFoam does
#include "readGravitationalAcceleration.H"
#include "readWaveProperties.H"
#include "createExternalWaveForcing.H"

// Total velocity, complementary velocity and complementary pressure
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject("U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading field Uc\n" << endl;
volVectorField Uc
(
    IOobject("Uc", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject("pc", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// Total flux, initialised from the total velocity
surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::flux(U)
);

// Two-phase mixture (VoF) and densities
Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

surfaceScalarField rhoPhi
(
    IOobject("rhoPhi", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(rho)*phi
);

// Transport model (laminar in the reference case), as waveFoam does
typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType(rho, U, phi, rhoPhi, mixture)
);

// SWENSE incident fields and relaxation weight
Info<< "Creating SWENSE incident fields\n" << endl;
swenseIncident incident(mesh, externalWave(), rho1);

Info<< "Creating SWENSE relaxation\n" << endl;
const dictionary& farfieldDict =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .subDict("farfieldCoeffs");
swenseRelaxation relax(mesh, farfieldDict.subDict("relaxationZone"));
const volScalarField& w = relax.weight();

// Incident pressure scaled to the mixture density, and the total pressure
volScalarField pStar
(
    IOobject("pStar", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    rho/rho1*incident.pInc()
);

// Reference height (constant/hRef, default 0) and gravity potential
uniformDimensionedScalarField hRef
(
    IOobject("hRef", runTime.constant(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::NO_WRITE),
    dimensionedScalar("hRef", dimLength, 0.0)
);

dimensionedScalar ghRef
(
    mag(g.value()) > SMALL
  ? g & (cmptMag(g.value())/mag(g.value()))*hRef
  : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0.0)
);

volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

volScalarField p
(
    IOobject("p", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE),
    pc + pStar
);

volScalarField p_rgh
(
    IOobject("p_rgh", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE),
    p - rho*gh
);

// Reference cell for pc from the pimple dictionary
label pcRefCell = 0;
scalar pcRefValue = 0.0;
setRefCell(pc, pimple.dict(), pcRefCell, pcRefValue);

if (pc.needReference())
{
    pc += dimensionedScalar("pc", pc.dimensions(), pcRefValue - getRefCellValue(pc, pcRefCell));
    p = pc + pStar;
    p_rgh = p - rho*gh;
}

// Complementary flux and incident flux
surfaceScalarField phiComp
(
    IOobject("phiComp", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::flux(Uc)
);

surfaceScalarField phiInc
(
    IOobject("phiInc", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    incident.UfInc() & mesh.Sf()
);

phi = phiComp + phiInc;
rhoPhi = fvc::interpolate(rho)*phi;

// Interface-pressure and momentum helper fields
surfaceVectorField Ucf(IOobject("Ucf", runTime.timeName(), mesh,
    IOobject::NO_READ, IOobject::NO_WRITE), fvc::interpolate(Uc));

bool includeDivUInc
(
    pimple.dict().lookupOrDefault<Switch>("includeDivUInc", false)
);

