/*---------------------------------------------------------------------------*\
Class
    Foam::swenseWaveAlphaFvPatchScalarField

Description
    Water fraction boundary condition for the SWENSE solver: the value on each
    boundary face is the wetted fraction of the face (wetted area / face area),
    clipped to [0, 1]. Mixed type with valueFraction 1.

    The wetted fraction uses the Grid2Grid external wave forcing registered as
    "externalWaveForcing".

SourceFiles
    swenseWaveAlphaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef swenseWaveAlphaFvPatchScalarField_H
#define swenseWaveAlphaFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "externalWaveForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class swenseWaveAlphaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class swenseWaveAlphaFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data


public:

    //- Runtime type information
    TypeName("swenseWaveAlpha");


    // Constructors

        //- Construct from patch and internal field
        swenseWaveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swenseWaveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given patch field onto a new patch
        swenseWaveAlphaFvPatchScalarField
        (
            const swenseWaveAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swenseWaveAlphaFvPatchScalarField
        (
            const swenseWaveAlphaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new swenseWaveAlphaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swenseWaveAlphaFvPatchScalarField
        (
            const swenseWaveAlphaFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new swenseWaveAlphaFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
