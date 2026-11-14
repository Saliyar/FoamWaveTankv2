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

// Start-up check (document Section 16, Part 1 rows 12-13): the complementary
// fields pc and Uc must not read the total U or the total flux phi. Such
// boundary conditions give a non-zero complementary value when the flow is
// exactly incident. Allowed: fixedValue, zeroGradient, fixedFluxPressure, etc.
{
    const wordList forbidden
    (
        {
            "totalPressure",
            "pressureInletOutletVelocity",
            "pressureInletVelocity",
            "pressureInletOutletParSlipVelocity",
            "flowRateInletVelocity",
            "inletOutlet",
            "outletInlet"
        }
    );

    forAll(pc.boundaryField(), patchi)
    {
        if (forbidden.found(pc.boundaryField()[patchi].type()))
        {
            FatalErrorInFunction
                << "pc on patch " << pc.boundaryField()[patchi].patch().name()
                << " uses " << pc.boundaryField()[patchi].type()
                << ", which reads the total U or phi. Use fixedValue or"
                << " zeroGradient for the complementary pressure."
                << exit(FatalError);
        }
    }

    forAll(Uc.boundaryField(), patchi)
    {
        if (forbidden.found(Uc.boundaryField()[patchi].type()))
        {
            FatalErrorInFunction
                << "Uc on patch " << Uc.boundaryField()[patchi].patch().name()
                << " uses " << Uc.boundaryField()[patchi].type()
                << ", which reads the total U or phi. Use fixedValue or"
                << " zeroGradient for the complementary velocity."
                << exit(FatalError);
        }
    }
}

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

// Blended incident fields (water and air, weight alphaI)
Info<< "Creating SWENSE incident fields\n" << endl;
swenseIncident incident(mesh, externalWave(), rho1, rho2);

// This revision (momentum_swense_Rev1, Part II) is for a fixed mesh: the
// spec does not cover moving meshes
if (mesh.dynamic())
{
    FatalErrorInFunction
        << "The momentum-blend SWENSE (Rev1) is implemented for a static mesh"
        << " only" << exit(FatalError);
}

// Wave models at the start time: eta for the initial alphaI and alphaH, and
// the water and air fields for the initial incident step at the end of this
// file
externalWave->step();
incident.evaluate(runTime.value());

// Incident flux (projected, spec step 1d): the flux of the alphaI transport
// and of the alphaI boundary inletOutlet. Set by the initial incident step
surfaceScalarField phiInc
(
    IOobject("phiInc", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    incident.UfInc() & mesh.Sf()
);

// Mesh flux (swept volumes per time step) of the moving mesh. Zero on a
// static mesh, so every relative flux below is the absolute flux there.
// The total flux phi and the incident flux phiIRel are relative to the mesh
// (as in the plain ALE solvers), while phiInc and phiComp stay absolute.
surfaceScalarField meshPhi
(
    IOobject("meshPhi", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("meshPhi", dimVolume/dimTime, 0.0)
);

// Relative incident flux phiInc - meshPhi, the flux that advects alphaI
surfaceScalarField phiIRel
(
    IOobject("phiIRel", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    phiInc - meshPhi
);

// Point positions at t = 0: a cell or a column that has not moved since then
// keeps its level correction and its relaxation weight (spec items 7 and 8)
const pointField points0(mesh.points());

// Moving-mesh switch: the outer correctors update the mesh (as waveFoam)
const bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
);

// Advected incident volume fraction alphaI (document Section 14.2). Its
// start value is the blend weight, and alpha starts equal to it (Theorem 1).
// Copy of alpha1 so that alphaI has the same boundary types (swenseWaveAlpha
// at the farfield, which takes the blend weight on the boundary faces)
// Boundary types come from 0/alphaI (the farfield swenseWaveAlpha, the
// atmosphere inletOutlet with phi phiInc, the walls); values are set below
volScalarField alphaI
(
    IOobject("alphaI", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// On a restart (start time > 0) alphaI and alpha1 keep the values
// read from the time directory; the blend initialisation is only for a start from t = 0.
if (runTime.value() < SMALL)
{
    alphaI.primitiveFieldRef() =
        incident.blendAlpha(runTime.value()).primitiveField();
    alphaI.correctBoundaryConditions();

    alpha1.primitiveFieldRef() = alphaI.primitiveField();
    alpha1.correctBoundaryConditions();
}
alpha2 = 1.0 - alpha1;
mixture.correct();

// Densities from the same function of the volume fraction (document 16.1)
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE),
    swenseRho(alpha1, rho1, rho2)
);
rho.oldTime();

volScalarField rhoI
(
    IOobject("rhoI", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    swenseRho(alphaI, rho1, rho2)
);
rhoI.oldTime();

// Density difference Drho = rho - rhoI (spec Section 12.2)
volScalarField Drho
(
    IOobject("Drho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    rho - rhoI
);

// HOS target for the relaxation zones, same initial value as alphaI
volScalarField alphaH
(
    IOobject("alphaH", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    alphaI
);

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

// SWENSE relaxation weight
Info<< "Creating SWENSE relaxation\n" << endl;
const dictionary& farfieldDict =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .subDict("farfieldCoeffs");
swenseRelaxation relax(mesh, farfieldDict.subDict("relaxationZone"));
const volScalarField& w = relax.weight();

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

// Total pressure p = pInc + pc + Drho gh (spec Section 12.4); pInc is set by
// the initial incident step below
volScalarField p
(
    IOobject("p", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE),
    incident.pInc() + pc + Drho*gh
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
    p = incident.pInc() + pc + Drho*gh;
    p_rgh = p - rho*gh;
}

// Complementary flux and incident flux
// phiComp is read from the time directory when present (written by
// AUTO_WRITE); otherwise it starts as fvc::flux(Uc), as before.
surfaceScalarField phiComp
(
    IOobject("phiComp", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    fvc::flux(Uc)
);


// Incident mass flux of alphaI (returned by its transport, document 16.3)
surfaceScalarField rhoPhiI
(
    IOobject("rhoPhiI", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(rhoI)*phiIRel
);

phi = phiComp + phiInc - meshPhi;
rhoPhi = fvc::interpolate(rho)*phi;

// Interface-pressure and momentum helper fields
surfaceVectorField Ucf(IOobject("Ucf", runTime.timeName(), mesh,
    IOobject::NO_READ, IOobject::NO_WRITE), fvc::interpolate(Uc));

// Ucf rebuilt from phiComp as in swensePcEqn.H (interpolate Uc, normal part from phiComp)
{
    Ucf = fvc::interpolate(Uc);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Ucf += n*(phiComp/mesh.magSf() - (n & Ucf));
}

// Potential of the incident flux projection (spec eq. 31): zero on the
// atmosphere patch, zero gradient on every other boundary
wordList psiProjTypes
(
    mesh.boundaryMesh().size(),
    zeroGradientFvPatchScalarField::typeName
);
forAll(mesh.boundaryMesh(), patchi)
{
    if (mesh.boundaryMesh()[patchi].name() == "atmosphere")
    {
        psiProjTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

volScalarField psiProj
(
    IOobject("psiProj", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("psiProj", dimDensity*dimArea/dimTime, 0.0),
    psiProjTypes
);

// The projection uses the flux of its matrix (as p_rgh in the plain solvers)
mesh.setFluxRequired(psiProj.name());

// Logging only: term budget switch and the raw-blend gradient of UInc (set in
// swenseIncidentStep.H before the projection)
const bool budgetOn = runTime.controlDict().lookupOrDefault<bool>("budgetTerms", false);
const bool theoremOn = runTime.controlDict().lookupOrDefault<bool>("theoremLog", false);

// Limited closure of the incident residual in the Uc and pc equations
// (waveProperties limitedClosure, default off)
const bool limitedClosure =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("limitedClosure", false);
scalarField budgetRawGradMag(mesh.nCells(), 0.0);

// Initial incident step (spec steps 1b to 1d) with the initial alphaI: the
// incident flux and the corrected incident velocity of level 0
#include "swenseIncidentStep.H"

incident.UIncRW().oldTime();

phi = phiComp + phiInc - meshPhi;
rhoPhi = fvc::interpolate(rho)*phi;
rhoPhiI = fvc::interpolate(rhoI)*phiIRel;

U == Uc + incident.UInc();
p = incident.pInc() + pc + Drho*gh;
p_rgh = p - rho*gh;
