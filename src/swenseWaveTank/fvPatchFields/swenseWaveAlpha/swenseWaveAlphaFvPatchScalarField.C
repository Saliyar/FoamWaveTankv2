/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseWaveAlphaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "swenseIncident.H"
#include "swenseRamp.H"
#include "swenseSurface.H"
#include "IOdictionary.H"

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

    // Blend weight of the water incident field at the face centres, the
    // same alphaI as the cell field (swenseIncident::blendWeight)
    const scalar delta = swenseIncident::blendThickness(mesh.thisDb());
    const scalar Tsoft = swenseRamp::Tsoft(mesh.thisDb());
    const scalar f = swenseRamp::factor(Tsoft, t);
    const vectorField& Cf = this->patch().Cf();

    // alphaHExact on: exact wetted face fraction (the cut of swenseSurface)
    if (swenseIncident::alphaHExact(mesh.thisDb()))
    {
        swenseSurface surf(mesh, ext, t);
        const pointField& P = mesh.points();
        const faceList& F = mesh.faces();
        const label start = this->patch().start();

        forAll(Cf, facei)
        {
            const scalar A = F[start + facei].mag(P);
            refValue()[facei] = A > VSMALL
              ? min(max(surf.wetArea(start + facei)/A, scalar(0)), scalar(1))
              : 0.0;
        }

        refGrad() = 0.0;
        valueFraction() = 1.0;
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    forAll(Cf, facei)
    {
        const point& x = Cf[facei];

        const scalar eta =
            Tsoft > 0 ? f*ext.eta(x, t) : ext.eta(x, t);

        refValue()[facei] =
            swenseIncident::blendWeight(x.z() - eta, delta);
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
