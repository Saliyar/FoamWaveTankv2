/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseWaveAlphaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "swenseSurface.H"

// * * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::swenseWaveAlphaFvPatchScalarField::swenseWaveAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::swenseWaveAlphaFvPatchScalarField::swenseWaveAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::swenseWaveAlphaFvPatchScalarField::swenseWaveAlphaFvPatchScalarField
(
    const swenseWaveAlphaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::swenseWaveAlphaFvPatchScalarField::swenseWaveAlphaFvPatchScalarField
(
    const swenseWaveAlphaFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::swenseWaveAlphaFvPatchScalarField::swenseWaveAlphaFvPatchScalarField
(
    const swenseWaveAlphaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::swenseWaveAlphaFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();
    const scalar t = db().time().value();

    const waveTheories::externalWaveForcing& ext =
        mesh.lookupObject<waveTheories::externalWaveForcing>
        (
            "externalWaveForcing"
        );

    swenseSurface surf(mesh, ext, t);

    const scalarField& magSf = this->patch().magSf();
    const label start = this->patch().patch().start();

    forAll(magSf, facei)
    {
        refValue()[facei] =
            min(1.0, max(0.0, surf.wetArea(facei + start)/magSf[facei]));
    }

    refGrad() = 0.0;
    valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::swenseWaveAlphaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * Runtime selection * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        swenseWaveAlphaFvPatchScalarField
    );
}

// ************************************************************************* //
