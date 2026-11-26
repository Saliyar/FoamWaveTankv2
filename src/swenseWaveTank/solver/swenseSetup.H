// Shared SWENSE set-up: fields, mixture, incident and relaxation objects.
// Included by both waveFoamSWENSE and waveIsoFoamSWENSE. The alpha-specific
// objects are created in the solver's createSolverFields.H.

// Gravity (OpenFOAM v2412 standard), waveProperties and the Grid2Grid
// external wave forcing (waves2Foam), as waveFoam does
#include "readGravitationalAcceleration.H"
#include "readWaveProperties.H"

// External wave forcing. rampMode factor (default): the forcing of
// externalForcing, as createExternalWaveForcing.H. rampMode height: the
// height-ramped analytic theory (swenseHeightRampForcing), the one place that
// supplies the ramped incident. Then the start-up consistency check of an
// analytic incident at the full height: kinematic and mass-flux residual over
// one wavelength / (omega H/2), stop above incidentResidualMax (default 0.005)
autoPtr<Foam::waveTheories::externalWaveForcing> externalWave
(
    swenseHeightRampForcing::newForcing(runTime, mesh)
);
swenseHeightRampForcing::residualCheck(externalWave(), mesh);

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
// Atmosphere brake: pc totalPressure with U Uc and phi phiComp, and Uc
// pressureInletOutletVelocity with phi phiComp, read only complementary fields
// (zero when the flow is exactly incident, so pc_b = p0 and pure waves are
// unchanged); they are allowed. With the total U or phi they stay forbidden.
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
        const fvPatchScalarField& pf = pc.boundaryField()[patchi];
        if (isA<totalPressureFvPatchScalarField>(pf))
        {
            const totalPressureFvPatchScalarField& tp =
                refCast<const totalPressureFvPatchScalarField>(pf);
            if (tp.UName() == "Uc" && tp.phiName() == "phiComp")
            {
                Info<< "SWENSE: pc on patch " << pf.patch().name()
                    << " totalPressure with U " << tp.UName() << " phi "
                    << tp.phiName() << " allowed (atmosphere brake)" << endl;
                continue;
            }
        }
        if (forbidden.found(pf.type()))
        {
            FatalErrorInFunction
                << "pc on patch " << pf.patch().name()
                << " uses " << pf.type()
                << ", which reads the total U or phi: pc would then be"
                << " non-zero when the flow is exactly incident. Use"
                << " fixedValue or zeroGradient, or totalPressure with"
                << " U Uc and phi phiComp (complementary fields only)."
                << exit(FatalError);
        }
    }

    forAll(Uc.boundaryField(), patchi)
    {
        const fvPatchVectorField& pf = Uc.boundaryField()[patchi];
        if (isA<pressureInletOutletVelocityFvPatchVectorField>(pf))
        {
            const pressureInletOutletVelocityFvPatchVectorField& pv =
                refCast<const pressureInletOutletVelocityFvPatchVectorField>(pf);
            if (pv.phiName() == "phiComp")
            {
                Info<< "SWENSE: Uc on patch " << pf.patch().name()
                    << " pressureInletOutletVelocity with phi "
                    << pv.phiName() << " allowed (atmosphere brake)" << endl;
                continue;
            }
        }
        if (forbidden.found(pf.type()))
        {
            FatalErrorInFunction
                << "Uc on patch " << pf.patch().name()
                << " uses " << pf.type()
                << ", which reads the total U or phi: Uc would then be"
                << " non-zero when the flow is exactly incident. Use"
                << " fixedValue or zeroGradient, or"
                << " pressureInletOutletVelocity with phi phiComp"
                << " (complementary flux only)."
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
    // alphaHExact on: exact wetted fraction (still water: equals the blend
    // weight when blendThickness = dz of the still-water row)
    alphaI.primitiveFieldRef() =
        incident.targetAlpha(runTime.value())().primitiveField();

    if (swenseIncident::alphaHExact(mesh.thisDb()))
    {
        Info<< "SWENSE alphaHExact: initial max|exact - blend| = "
            << gMax(mag(alphaI.primitiveField()
                 - incident.blendAlpha(runTime.value()).primitiveField()))
            << endl;
    }
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

// Edge of the body zone (waveProperties plainBodyZoneEdge, default min):
//   min  : lam_f = min of owner and neighbour (zone-edge faces get lam_f = 0)
//   max  : lam_f = max on the zone-edge faces (edge faces keep lam_f = 1)
//   ramp : lam rises 0 -> 1 over plainBodyZoneRampCells cell layers outside
//          the zone (layer k: lam = k/(N+1)); lam_f = min as before
const word plainBodyZoneEdge =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<word>("plainBodyZoneEdge", "min");

const label plainBodyZoneRampCells =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<label>("plainBodyZoneRampCells", 3);

if
(
    plainBodyZoneEdge != "min" && plainBodyZoneEdge != "max"
 && plainBodyZoneEdge != "ramp"
)
{
    FatalErrorInFunction
        << "plainBodyZoneEdge must be min, max or ramp, not "
        << plainBodyZoneEdge << exit(FatalError);
}

// Face residual in the body zone (waveProperties plainBodyZoneFace):
//   full       : (1-lam_f) G_f, as before (default)
//   solenoidal : (1-lam_f) (G_f - snGrad(phiG)), lap(phiG) = div(G) per step
const word plainBodyZoneFace =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<word>("plainBodyZoneFace", "full");

if (plainBodyZoneFace != "full" && plainBodyZoneFace != "solenoidal")
{
    FatalErrorInFunction
        << "plainBodyZoneFace must be full or solenoidal, not "
        << plainBodyZoneFace << exit(FatalError);
}
Info<< "SWENSE plainBodyZoneFace " << plainBodyZoneFace << endl;

// Topological layer of each cell from the zone: 0 in the zone, k for the
// k-th face-neighbour layer outside, labelMax beyond (serial: cellCells)
labelList zoneLayer(mesh.nCells(), labelMax);

if (plainBodyZoneID >= 0)
{
    labelList front(mesh.cellZones()[plainBodyZoneID]);

    for (const label c : front)
    {
        zoneLayer[c] = 0;
    }

    for (label k = 1; k <= max(plainBodyZoneRampCells, label(4)); ++k)
    {
        DynamicList<label> next;

        for (const label c : front)
        {
            for (const label n : mesh.cellCells()[c])
            {
                if (zoneLayer[n] == labelMax)
                {
                    zoneLayer[n] = k;
                    next.append(n);
                }
            }
        }

        front.transfer(next);
    }

    Info<< "SWENSE plainBodyZoneEdge " << plainBodyZoneEdge;
    if (plainBodyZoneEdge == "ramp")
    {
        Info<< " over " << plainBodyZoneRampCells << " cells";
    }
    Info<< endl;
}

// Per-step log of the zone source split (waveProperties logZoneSource,
// default off): max water |Uc| outside the zone and, around that cell, the
// cell closure (1-lam) rI_c, the face term R[(1-lam_f) G_f] and the rest of
// the incident term. Logging only.
const bool logZoneSource =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("logZoneSource", false);
scalar zoneSourceThreshold =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<scalar>("zoneSourceThreshold", 0.05);
bool zoneSourceCrossed = false;

// Discrete incident hydrostatics (waveProperties incidentHydrostatic
// discrete, default off; static or moving mesh). pInc = pHydI + pDyn with
//   pDyn       = p_side - rho_side g (eta - z)   (swenseIncident::sidePressure)
//   pHydI_rgh  = pHydI - rhoI gh, marched down each column from the top:
//                pHydI_rgh(lower) = pHydI_rgh(upper) - gh_f (rhoI(lower) - rhoI(upper))
//                (the gh_f of G_f on that face), top cell = analytic
//                rho_side g (eta - z) - rhoI gh.
// The vertical face residual G_f then has no hydrostatic part, bit for bit.
// swenseHydroDiscrete.H rebuilds the fields; the closure face term and p use
// pIncD. Off: nothing is built.
const word incidentHydrostatic =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<word>("incidentHydrostatic", "analytic");

if (incidentHydrostatic != "analytic" && incidentHydrostatic != "discrete")
{
    FatalErrorInFunction
        << "incidentHydrostatic must be analytic or discrete, not "
        << incidentHydrostatic << exit(FatalError);
}

const bool hydroDiscrete = (incidentHydrostatic == "discrete");

// Atmosphere value of pc consistent with p = 0 there (waveProperties
// pcAtmosphereConsistent, default off). p = pInc + pc + Drho gh, so with the
// incident pressure zero on the atmosphere patch pc_b = -Drho_b gh_b. Set in
// swensePcEqn.H before constrainPressure; pc stays fixedValue. Off: pc_b is
// the value of 0/pc (0), unchanged.
const bool pcAtmosphereConsistent =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("pcAtmosphereConsistent", false);

const label pcAtmPatch = mesh.boundaryMesh().findPatchID("atmosphere");

if (pcAtmosphereConsistent && pcAtmPatch < 0)
{
    FatalErrorInFunction
        << "pcAtmosphereConsistent needs a patch named atmosphere"
        << exit(FatalError);
}

Info<< "SWENSE pcAtmosphereConsistent " << pcAtmosphereConsistent << endl;

// Column tables, built once on the mesh at set-up (topology: the cell
// through the top face of each cell, found by the face normal at set-up;
// -1 where the top face is a boundary: atmosphere, or a body wall). Cells in
// order of falling set-up z, which keeps every cell after the one above it
// whatever the later motion
labelList colAbove(mesh.nCells(), -1);
labelList colFace(mesh.nCells(), -1);
labelList colOrder;

volScalarField pHydA
(
    IOobject("pHydA", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(p.dimensions(), Zero)
);
volScalarField pDynD
(
    IOobject("pDynD", runTime.timeName(), mesh,
        IOobject::NO_READ, hydroDiscrete ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(p.dimensions(), Zero)
);
volScalarField pHrgh
(
    IOobject("pHydIrgh", runTime.timeName(), mesh,
        IOobject::NO_READ, hydroDiscrete ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(p.dimensions(), Zero)
);
volScalarField pIncD
(
    IOobject("pIncD", runTime.timeName(), mesh,
        IOobject::NO_READ, hydroDiscrete ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    incident.pInc()
);

if (hydroDiscrete)
{
    if (mag(g.value().x()) > SMALL || mag(g.value().y()) > SMALL)
    {
        FatalErrorInFunction
            << "incidentHydrostatic discrete needs g along z" << exit(FatalError);
    }

    const vectorField& Cc = mesh.C().primitiveField();
    const scalar tolN = 0.999;

    for (label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        const vector nf = mesh.Sf()[f]/mesh.magSf()[f];

        if (mag(nf.z()) < tolN)
        {
            continue;
        }

        label lo = mesh.owner()[f];
        label up = mesh.neighbour()[f];

        if (Cc[up].z() < Cc[lo].z())
        {
            Swap(lo, up);
        }

        if (colAbove[lo] >= 0)
        {
            FatalErrorInFunction
                << "incidentHydrostatic discrete: cell " << lo
                << " has two cells above" << exit(FatalError);
        }

        colAbove[lo] = up;
        colFace[lo] = f;
    }

    // A column cut by a processor boundary is not allowed (serial or
    // whole-column decomposition)
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if (!pp.coupled())
        {
            continue;
        }

        forAll(pp, i)
        {
            const vector nf = pp.faceAreas()[i]/mag(pp.faceAreas()[i]);

            if (mag(nf.z()) >= tolN)
            {
                FatalErrorInFunction
                    << "incidentHydrostatic discrete: patch " << pp.name()
                    << " cuts a column; use serial or whole-column"
                    << " decomposition" << exit(FatalError);
            }
        }
    }

    scalarField negZ(mesh.nCells());
    forAll(negZ, c)
    {
        negZ[c] = -Cc[c].z();
    }
    colOrder = sortedOrder(negZ);

    label nTop = 0;
    forAll(colAbove, c)
    {
        if (colAbove[c] < 0)
        {
            ++nTop;
        }
    }

    Info<< "SWENSE incidentHydrostatic discrete: column starts " << nTop
        << (mesh.dynamic() ? " (moving mesh: topology fixed at set-up)" : "")
        << endl;
}

// Column level correction of the incident surface (controlDict
// levelCorrection, default off; swenseColumnLevel.H, called after the
// relaxation in swenseAlphaRelax.H). With plainBodyZone or a moving mesh it
// is allowed only in columns wholly outside the zone (checked at set-up).
const bool levelCorrection =
    runTime.controlDict().lookupOrDefault<bool>("levelCorrection", false);
Info<< "SWENSE levelCorrection " << levelCorrection << endl;

// Exact incident reset (controlDict alphaIExactReset, default off;
// swenseExactReset.H, called in swenseAlphaRelax.H after the relaxation, in
// place of the level correction). alphaI is set to swenseIncident::exactAlpha
// outside plainBodyZone, and alpha1 takes the same increment. The relaxation
// target must be the same function (waveProperties alphaHExact on), otherwise
// the zones and the reset pull alphaI toward different values every step.
const bool alphaIExactReset =
    runTime.controlDict().lookupOrDefault<bool>("alphaIExactReset", false);
Info<< "SWENSE alphaIExactReset " << alphaIExactReset << endl;
if (alphaIExactReset && !swenseIncident::alphaHExact(mesh.thisDb()))
{
    FatalErrorInFunction
        << "controlDict alphaIExactReset needs waveProperties alphaHExact on:"
        << " the reset and the relaxation target must both be"
        << " swenseIncident::exactAlpha" << exit(FatalError);
}
if (alphaIExactReset && levelCorrection)
{
    Info<< "SWENSE alphaIExactReset on: levelCorrection is skipped" << endl;
}

// Cells left unchanged by the reset (plainBodyZone) and running totals
boolList exactResetSkip(alphaIExactReset ? mesh.nCells() : 0, false);
if (alphaIExactReset && plainBodyZoneID >= 0)
{
    for (const label c : mesh.cellZones()[plainBodyZoneID])
    {
        exactResetSkip[c] = true;
    }
    Info<< "SWENSE alphaIExactReset: cells in plainBodyZone left unchanged: "
        << returnReduce(mesh.cellZones()[plainBodyZoneID].size(), sumOp<label>())
        << endl;
}
scalar exactResetTotal = 0;
scalar exactResetClipTotal = 0;

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

// Incident alphaI on body walls (waveProperties incidentAlphaWall, default
// off): the alphaI patch of every body patch must be swenseIncidentAlpha.
// alphaStep.H evaluates these patches after alphaI = alphaI.oldTime(), so
// inflow faces (phiIRel into the fluid) carry the wetted face fraction at
// t^{n+1/2} into the upwind boundary flux of isoAdvector. Off: the patch
// condition is exactly zeroGradient and nothing extra is evaluated.
const bool incidentAlphaWall =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("incidentAlphaWall", false);

if (incidentAlphaWall)
{
    forAll(bodyPatches, i)
    {
        const fvPatchScalarField& ab = alphaI.boundaryField()[bodyPatches[i]];
        if (ab.type() != "swenseIncidentAlpha")
        {
            FatalErrorInFunction
                << "incidentAlphaWall on: alphaI on body patch "
                << mesh.boundary()[bodyPatches[i]].name()
                << " must be swenseIncidentAlpha, not " << ab.type()
                << exit(FatalError);
        }
    }
    Info<< "SWENSE incidentAlphaWall on, body patches " << bodyPatches
        << ", advecting flux "
        << (swenseNoRelFlux ? "phiInc (set phi phiInc in 0/alphaI)" : "phiIRel")
        << endl;
}

// Incident reset in the body zone (waveProperties incidentResetInZone,
// default off): after the alphaI step, alphaI = exact wetted cell fraction in
// the plainBodyZone cells, with rhoI, uI, pInc rebuilt there and Uc, pc
// shifted so that the totals are unchanged (swenseResetInZone.H).
const bool incidentResetInZone =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("incidentResetInZone", false);

if (incidentResetInZone && plainBodyZoneID < 0)
{
    FatalErrorInFunction
        << "incidentResetInZone needs plainBodyZone" << exit(FatalError);
}

// Per-step log of the smallest body-wall to zone-edge distance in cells
// (waveProperties logZoneEdge, default off; swenseZoneEdgeLog.H)
const bool logZoneEdge =
    mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        .lookupOrDefault<Switch>("logZoneEdge", false);
bool zoneEdgeBelow1 = false;

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
