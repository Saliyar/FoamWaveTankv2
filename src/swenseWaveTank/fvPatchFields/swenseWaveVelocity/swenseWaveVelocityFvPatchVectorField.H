/*---------------------------------------------------------------------------*\
Class
    Foam::swenseWaveVelocityFvPatchVectorField

Description
    Velocity boundary condition for the SWENSE solver.

    For each boundary face: if the face is dry (no wetted area), the value is
    the wind velocity; otherwise the value is the wave velocity at the centre
    of the wetted part of the face. Mixed type with valueFraction 1.

    The wave velocity is the Grid2Grid external wave forcing registered as
    "externalWaveForcing". The wind and rampTime are read from the
    waveProperties sub-dictionary named <patch>Coeffs. The ramp is not
    supported: rampTime must be 0 (or absent).

SourceFiles
    swenseWaveVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef swenseWaveVelocityFvPatchVectorField_H
#define swenseWaveVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"
#include "externalWaveForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class swenseWaveVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class swenseWaveVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data


public:

    //- Runtime type information
    TypeName("swenseWaveVelocity");


    // Constructors

        //- Construct from patch and internal field
        swenseWaveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swenseWaveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given patch field onto a new patch
        swenseWaveVelocityFvPatchVectorField
        (
            const swenseWaveVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swenseWaveVelocityFvPatchVectorField
        (
            const swenseWaveVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new swenseWaveVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swenseWaveVelocityFvPatchVectorField
        (
            const swenseWaveVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new swenseWaveVelocityFvPatchVectorField(*this, iF)
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
