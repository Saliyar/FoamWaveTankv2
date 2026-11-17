/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseBodyVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    exactWallFlux_(false)
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    exactWallFlux_(dict.getOrDefault<Switch>("exactWallFlux", false))
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const swenseBodyVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    exactWallFlux_(ptf.exactWallFlux_)
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const swenseBodyVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    exactWallFlux_(ptf.exactWallFlux_)
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const swenseBodyVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    exactWallFlux_(ptf.exactWallFlux_)
{}


// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::swenseBodyVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Registered incident velocity (swenseIncident::UInc_), on this patch
    const volVectorField& UInc =
        db().lookupObject<volVectorField>("UInc");

    // Registered mesh flux (swept volumes per step, zero on a static mesh).
    const surfaceScalarField& meshPhi =
        db().lookupObject<surfaceScalarField>("meshPhi");

    const label patchi = patch().index();

    // Full wall velocity as movingWallVelocity: the tangential part from the
    // face centre motion, (Cf - Cf_old)/dt, and the normal part replaced by
    // meshPhi/magSf. The total U on the body is this wall velocity, so
    // Uc = Uwall - UInc.
    const vectorField n(patch().nf());
    const scalarField& magSf = patch().magSf();
    const vectorField& UIb = UInc.boundaryField()[patchi];

    // Face-centre velocity of the wall, zero on a static mesh
    vectorField Up(patch().size(), Zero);

    if (patch().boundaryMesh().mesh().moving())
    {
        const fvMesh& mesh = patch().boundaryMesh().mesh();
        const polyPatch& pp = patch().patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());
        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        const scalar deltaT = mesh.time().deltaTValue();
        Up = (pp.faceCentres() - oldFc)/deltaT;
    }

    vectorField Ucb(patch().size(), Zero);

    if (exactWallFlux_ && db().foundObject<surfaceScalarField>("phiInc"))
    {
        // Option exactWallFlux: normal part (meshPhi - phiInc)/|Sf|, so that
        // phiComp + phiInc = meshPhi on the patch (phiComp_b = Uc_b . Sf
        // through constrainPressure); tangential part of Uwall - UInc. With
        // phiInc = UInc = 0 this is Uwall term by term.
        const surfaceScalarField& phiInc =
            db().lookupObject<surfaceScalarField>("phiInc");

        const scalarField Unc
        (
            (meshPhi.boundaryField()[patchi] - phiInc.boundaryField()[patchi])
           /(magSf + VSMALL)
        );

        Ucb = Up + n*(Unc - (n & Up));
        Ucb -= (UIb - n*(n & UIb));
    }
    else
    {
        // Off (default): Uc = Uwall - UInc. The relative wall flux is then
        // (UfInc - UInc).Sf, the normal part of the incident projection
        // correction of the wall cell (zero without an incident field).
        if (patch().boundaryMesh().mesh().moving())
        {
            const scalarField Un
            (
                meshPhi.boundaryField()[patchi]/(magSf + VSMALL)
            );

            Ucb = Up + n*(Un - (n & Up));
        }

        Ucb = Ucb - UIb;
    }

    operator==(Ucb);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::swenseBodyVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    os.writeEntryIfDifferent<Switch>("exactWallFlux", false, exactWallFlux_);
}


// * * * * * * * * * * * * * * * Build Run-Time Selection * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        swenseBodyVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
