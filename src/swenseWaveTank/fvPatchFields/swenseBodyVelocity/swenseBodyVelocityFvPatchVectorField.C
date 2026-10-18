/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseBodyVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

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

    operator==(-UInc.boundaryField()[patch().index()]);

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
