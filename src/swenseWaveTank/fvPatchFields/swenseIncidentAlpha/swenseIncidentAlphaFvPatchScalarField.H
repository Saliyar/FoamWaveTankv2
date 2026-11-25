/*---------------------------------------------------------------------------*\
Class
    Foam::swenseIncidentAlphaFvPatchScalarField

Description
    alphaI boundary condition for body patches (waveProperties
    incidentAlphaWall, default off).

    The incident flux phiIRel crosses a body wall (the incident wave does not
    see the body). With zeroGradient, a face where phiIRel flows from the wall
    into the fluid imports the wall cell's own alphaI, so the wall is a source
    and sink of alphaI. This condition gives those faces the incident value
    instead:
        phiIRel_f <  0 (into the fluid):  alphaI_f = A_wet(eta, t^{n+1/2}) / |S_f|
        phiIRel_f >= 0 (out of the fluid): alphaI_f = alphaI_c (zeroGradient)
    A_wet is the exact wetted area of the face below the incident surface
    (ramped eta of the externalWaveForcing, as in swenseIncident::etaAt), cut
    on the current face points (moving mesh), at t^{n+1/2} = t - dt/2 of the
    step that is advected. isoAdvector takes the upwind boundary value on
    inflow faces (upwind flux of the patch value) and corrects only outflow
    boundary faces, so this value is what is advected.

    Mixed type: valueFraction 1 on inflow faces, 0 (refGrad 0) on outflow
    faces. With the switch off all faces have valueFraction 0, i.e. exactly
    zeroGradient.

Usage
    box
    {
        type    swenseIncidentAlpha;
        phi     phiIRel;    // optional, the flux that advects alphaI
    }

SourceFiles
    swenseIncidentAlphaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef swenseIncidentAlphaFvPatchScalarField_H
#define swenseIncidentAlphaFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "swenseSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class swenseIncidentAlphaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class swenseIncidentAlphaFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the flux that advects alphaI
        word phiName_;

        //- Cut-face geometry against the incident surface (built on first use)
        autoPtr<swenseSurface> surf_;

        //- Inflow faces and wetted inflow fraction of the last update (log)
        label nIn_;
        scalar wetIn_;


public:

    //- Runtime type information
    TypeName("swenseIncidentAlpha");


    // Constructors

        swenseIncidentAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        swenseIncidentAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        swenseIncidentAlphaFvPatchScalarField
        (
            const swenseIncidentAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        swenseIncidentAlphaFvPatchScalarField
        (
            const swenseIncidentAlphaFvPatchScalarField&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new swenseIncidentAlphaFvPatchScalarField(*this)
            );
        }

        swenseIncidentAlphaFvPatchScalarField
        (
            const swenseIncidentAlphaFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new swenseIncidentAlphaFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Inflow faces and sum of their wetted area of the last update
        label nInflow() const { return nIn_; }
        scalar wetInflow() const { return wetIn_; }

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
