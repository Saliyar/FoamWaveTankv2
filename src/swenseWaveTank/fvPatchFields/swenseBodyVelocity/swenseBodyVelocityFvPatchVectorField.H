/*---------------------------------------------------------------------------*\
Class
    Foam::swenseBodyVelocityFvPatchVectorField
Description
    Velocity boundary condition for a body in the SWENSE solver (fixed value).
    The value is -UInc on the patch, so that the total velocity U = Uc + UInc
    is zero on the body (no-slip for the total field). UInc is the registered
    incident velocity field of the SWENSE solver, evaluated on the patch.
SourceFiles
    swenseBodyVelocityFvPatchVectorField.C
\*---------------------------------------------------------------------------*/

#ifndef swenseBodyVelocityFvPatchVectorField_H
#define swenseBodyVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class swenseBodyVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class swenseBodyVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("swenseBodyVelocity");


    // Constructors

        //- Construct from patch and internal field
        swenseBodyVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swenseBodyVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given patch field onto a new patch
        swenseBodyVelocityFvPatchVectorField
        (
            const swenseBodyVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swenseBodyVelocityFvPatchVectorField
        (
            const swenseBodyVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new swenseBodyVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swenseBodyVelocityFvPatchVectorField
        (
            const swenseBodyVelocityFvPatchVectorField&,
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
                new swenseBodyVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients: value = -UInc on the patch
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
