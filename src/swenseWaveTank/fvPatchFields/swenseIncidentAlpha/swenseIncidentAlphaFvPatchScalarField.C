/*---------------------------------------------------------------------------*\
    swenseIncidentAlpha: alphaI on body patches, incident wetted fraction on
    inflow faces, zeroGradient on outflow faces (see the header).
\*---------------------------------------------------------------------------*/

#include "swenseIncidentAlphaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swenseIncidentAlphaFvPatchScalarField::swenseIncidentAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phiIRel"),
    surf_(),
    nIn_(0),
    wetIn_(0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::swenseIncidentAlphaFvPatchScalarField::swenseIncidentAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phiIRel")),
    surf_(),
    nIn_(0),
    wetIn_(0)
{
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }
    refValue() = *this;
}


Foam::swenseIncidentAlphaFvPatchScalarField::swenseIncidentAlphaFvPatchScalarField
(
    const swenseIncidentAlphaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    surf_(),
    nIn_(0),
    wetIn_(0)
{}


Foam::swenseIncidentAlphaFvPatchScalarField::swenseIncidentAlphaFvPatchScalarField
(
    const swenseIncidentAlphaFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    surf_(),
    nIn_(ptf.nIn_),
    wetIn_(ptf.wetIn_)
{}


Foam::swenseIncidentAlphaFvPatchScalarField::swenseIncidentAlphaFvPatchScalarField
(
    const swenseIncidentAlphaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    surf_(),
    nIn_(ptf.nIn_),
    wetIn_(ptf.wetIn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swenseIncidentAlphaFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const bool on =
        mesh.thisDb().foundObject<IOdictionary>("waveProperties")
     && mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
            .getOrDefault<Switch>("incidentAlphaWall", false);

    refGrad() = 0.0;
    nIn_ = 0;
    wetIn_ = 0;

    if
    (
        !on
     || !db().foundObject<surfaceScalarField>(phiName_)
     || !mesh.thisDb().foundObject<waveTheories::externalWaveForcing>
        (
            "externalWaveForcing"
        )
    )
    {
        // zeroGradient
        valueFraction() = 0.0;
        refValue() = patchInternalField();
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    // Mid-step time of the step that ends at the current time
    const scalar tMid =
        db().time().value() - 0.5*db().time().deltaTValue();

    if (!surf_)
    {
        surf_.reset
        (
            new swenseSurface
            (
                mesh,
                mesh.thisDb().lookupObject<waveTheories::externalWaveForcing>
                (
                    "externalWaveForcing"
                ),
                tMid
            )
        );
    }
    surf_->setTime(tMid);

    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();
    const label start = patch().start();

    forAll(phip, i)
    {
        if (phip[i] < 0)
        {
            const label facei = start + i;
            const scalar A = faces[facei].mag(pts);
            const scalar Aw = surf_->wetArea(facei);

            refValue()[i] = A > VSMALL ? min(max(Aw/A, 0.0), 1.0) : 0.0;
            valueFraction()[i] = 1.0;

            ++nIn_;
            wetIn_ += refValue()[i];
        }
        else
        {
            refValue()[i] = 0.0;
            valueFraction()[i] = 0.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::swenseIncidentAlphaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phiIRel", phiName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        swenseIncidentAlphaFvPatchScalarField
    );
}

// ************************************************************************* //
