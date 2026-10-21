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
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
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
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const swenseBodyVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf)
{}


Foam::swenseBodyVelocityFvPatchVectorField::
swenseBodyVelocityFvPatchVectorField
(
    const swenseBodyVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
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
    vectorField Uwall(patch().size(), Zero);

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
        const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

        const vectorField n(patch().nf());
        const scalarField& magSf = patch().magSf();
        const scalarField Un
        (
            meshPhi.boundaryField()[patchi]/(magSf + VSMALL)
        );

        Uwall = Up + n*(Un - (n & Up));
    }

    operator==(Uwall - UInc.boundaryField()[patchi]);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::swenseBodyVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
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
