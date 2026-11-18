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

// Moving meshes (dyn port, dyn_spec.md): a dynamic mesh runs the incident
// step inside the PIMPLE loop, after mesh.update(). A static mesh keeps the
// original order (see swenseTimeLoop.H).

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

// Flux correction after a mesh motion (swenseCorrectPhi.H), as waveIsoFoam
const bool correctPhi
(
    pimple.dict().lookupOrDefault<Switch>("correctPhi", false)
);

// Diagnostic switches for the moving-mesh terms (waveProperties keys, default
// off: the normal solver path is unchanged). Used to find the term behind a
// gauge difference against a static run.
const dictionary& swenseSwitchDict =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties");

// a) alphaI and alpha1 transport with phiI and phi (absolute), not minus meshPhi
const bool swenseNoRelFlux
(
    swenseSwitchDict.lookupOrDefault<Switch>("swenseNoRelFlux", false)
);

// b) Uc equation: incident time derivative without the moving-mesh volume factor
const bool swenseNoMeshDdtUc
(
    swenseSwitchDict.lookupOrDefault<Switch>("swenseNoMeshDdtUc", false)
);

// c) incident projection skipped on a moved mesh (phiInc = phiRaw)
const bool swenseNoProjMesh
(
    swenseSwitchDict.lookupOrDefault<Switch>("swenseNoProjMesh", false)
);

// d) correctPhi flux correction skipped after motion
const bool swenseNoCorrectPhi
(
    swenseSwitchDict.lookupOrDefault<Switch>("swenseNoCorrectPhi", false)
);

// e) relaxation-weight update frozen at its step-0 value
const bool swenseNoWeightUpdate
(
    swenseSwitchDict.lookupOrDefault<Switch>("swenseNoWeightUpdate", false)
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

// Flux correction after a mesh motion (swenseCorrectPhi.H): the target is
// div(phiComp + phiInc) = 0. The face coefficient is the rAUf of the last
// pc solve (set at the end of swensePcEqn.H), as interFoam's correctPhi.H
// uses fvc::interpolate(rAU()); before the first pc solve it is 1 with the
// dimensions of 1/A, as interFoam's initCorrectPhi.H.
const bool includeDivUInc = true;

surfaceScalarField rAUfPrev
(
    IOobject("rAUfPrev", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("rAUfPrev", dimTime/rho.dimensions(), 1.0)
);

volScalarField pcCorr
(
    IOobject("pcCorr", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("pcCorr", pc.dimensions(), 0.0),
    pc.boundaryField().types()
);

// Absolute total flux phi + meshPhi, used by the alpha1 transport when
// swenseNoRelFlux is on (updated at the start of alphaStep.H)
surfaceScalarField phiAbsTot
(
    IOobject("phiAbsTot", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    phi + meshPhi
);

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

// Density-ratio cut of the limited closure: lam = 0 where
// max(rho,rhoI)/min(rho,rhoI) > limitedClosureRatio (default 2)
const scalar limitedClosureRatio =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<scalar>("limitedClosureRatio", 2.0);

// Continuous limited closure (waveProperties limitedClosureContinuous,
// default off): lam = 1 - |alpha1 - alphaI| instead of the 0/1 cut, i.e.
// the share of the incident residual kept is the mismatched mass
// |rho - rhoI|/(rho1 - rho2). lam = 1 where alpha1 = alphaI (residual
// dropped, invariance kept), lam = 0 where the phases are swapped (residual
// kept, total equation exact there). No free parameter, and linear in the
// fields, so no step switches the O(rho g dz/d) interface residual on or off.
// (A density-ratio form min(rho,rhoI)/max(rho,rhoI) was tried first and
// rejected: 1% of water in an incident-air cell already gives lam = 0.1.)
// limitedClosureRatio is not used in this mode.
const bool limitedClosureContinuous =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("limitedClosureContinuous", false);

// Plain closure in a body zone (waveProperties plainBodyZone, default none):
// the name of a cellZone (made in the case, e.g. with topoSet, from the cells
// that the body motion can move). lam = 0 in every zone cell, so lam_f = 0 on
// every face with the owner or the neighbour in the zone; outside the zone
// lam is the limited closure as before. Needs limitedClosure on.
const word plainBodyZoneName =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<word>("plainBodyZone", "none");

label plainBodyZoneID = -1;

if (plainBodyZoneName != "none")
{
    plainBodyZoneID = mesh.cellZones().findZoneID(plainBodyZoneName);

    if (plainBodyZoneID < 0)
    {
        FatalErrorInFunction
            << "plainBodyZone " << plainBodyZoneName
            << " is not a cellZone of the mesh; cellZones: "
            << mesh.cellZones().names() << exit(FatalError);
    }

    if (!limitedClosure)
    {
        FatalErrorInFunction
            << "plainBodyZone " << plainBodyZoneName
            << " needs limitedClosure on" << exit(FatalError);
    }

    Info<< "SWENSE plainBodyZone " << plainBodyZoneName << ": "
        << returnReduce(mesh.cellZones()[plainBodyZoneID].size(), sumOp<label>())
        << " cells with lam = 0" << endl;
}

// Write the incident residual rI at write times (waveProperties
// logIncidentResidual, default off). Logging only.
const bool logIncidentResidual =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("logIncidentResidual", false);

// Incident pressure balance (waveProperties incidentPressureBalance, default
// off; momentum_swense Part V, items V1 to V5). After the last pressure
// corrector of a step the flux of the incident residual r_I is balanced by a
// pressure correction deltaHeld:
//   laplacian(rAUf, deltaHeld) = div(Phi[r_I]),  pIncHeld = pInc - deltaHeld
// (swenseBalanceFlux.H, swenseBalanceSolve.H). p is rebuilt with pIncHeld
// (swensePcEqn.H), the closure face term of the next step uses pIncHeld
// (swenseUcEqn.H), and swenseBalanceLog.H logs the result. Off: nothing is
// built, solved or changed.
const bool incidentPressureBalance =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("incidentPressureBalance", false);

// deltaHeld boundary types: the constraint type on constraint patches
// (empty), fixedValue 0 where pc fixes its value, fixed gradient elsewhere
// (gradient set before each solve so that the boundary flux of
// Phi[r_I] - rAUf snGrad(deltaHeld) |Sf| is zero, as fixedFluxPressure)
wordList deltaHeldTypes
(
    pc.boundaryField().size(),
    fixedGradientFvPatchScalarField::typeName
);

forAll(pc.boundaryField(), patchi)
{
    const word& ptype = mesh.boundaryMesh()[patchi].type();

    if (polyPatch::constraintType(ptype))
    {
        deltaHeldTypes[patchi] = ptype;
    }
    else if (pc.boundaryField()[patchi].fixesValue())
    {
        deltaHeldTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

volScalarField deltaHeld
(
    IOobject
    (
        "deltaHeld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        incidentPressureBalance ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("deltaHeld", pc.dimensions(), 0.0),
    deltaHeldTypes
);

// Body patches (Uc of type swenseBodyVelocity), for the logged integral of
// deltaHeld n dS
DynamicList<label> bodyPatchList;
forAll(Uc.boundaryField(), patchi)
{
    if (Uc.boundaryField()[patchi].type() == "swenseBodyVelocity")
    {
        bodyPatchList.append(patchi);
    }
}
const labelList bodyPatches(bodyPatchList);

if (incidentPressureBalance)
{
    // V2 takes the face flux of the deltaHeld equation
    mesh.setFluxRequired(deltaHeld.name());

    Info<< "SWENSE incidentPressureBalance on; deltaHeld boundary types "
        << deltaHeldTypes << "; body patches " << bodyPatches << endl;
}

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
