/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseWaveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "swenseSurface.H"
#include "swenseRamp.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::swenseWaveVelocityFvPatchVectorField::
swenseWaveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


Foam::swenseWaveVelocityFvPatchVectorField::
swenseWaveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1.0;
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::swenseWaveVelocityFvPatchVectorField::
swenseWaveVelocityFvPatchVectorField
(
    const swenseWaveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::swenseWaveVelocityFvPatchVectorField::
swenseWaveVelocityFvPatchVectorField
(
    const swenseWaveVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchVectorField(ptf)
{}


Foam::swenseWaveVelocityFvPatchVectorField::
swenseWaveVelocityFvPatchVectorField
(
    const swenseWaveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::swenseWaveVelocityFvPatchVectorField::updateCoeffs()
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

    const dictionary& pDict =
        mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .subDict(this->patch().name() + "Coeffs");

    if (pDict.lookupOrDefault<scalar>("rampTime", 0.0) > SMALL)
    {
        FatalErrorInFunction
            << "swenseWaveVelocity does not support rampTime on patch "
            << this->patch().name() << exit(FatalError);
    }

    const vector wind = pDict.lookupOrDefault<vector>("wind", vector::zero);

    swenseSurface surf(mesh, ext, t);
    const scalar Tsoft = swenseRamp::Tsoft(mesh.thisDb());
    const scalar f = swenseRamp::factor(Tsoft, t);

    const scalarField& magSf = this->patch().magSf();
    const label start = this->patch().patch().start();

    vectorField& patchU = refValue();

    forAll(magSf, facei)
    {
        const label faceID = facei + start;
        const scalar wetArea = surf.wetArea(faceID);

        if (wetArea < SMALL*magSf[facei])
        {
            patchU[facei] = wind;
        }
        else
        {
            patchU[facei] = ext.U(surf.wetCentre(faceID), t);

            if (Tsoft > 0)
            {
                patchU[facei] *= f;
            }
        }
    }

    valueFraction() = 1.0;

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::swenseWaveVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * Runtime selection * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        swenseWaveVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
