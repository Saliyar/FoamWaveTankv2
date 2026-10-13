/*---------------------------------------------------------------------------*\
Class
    Foam::swenseRelaxation

Description
    Relaxation weight for the SWENSE solver: a farfield weight built from the
    distance of each cell to the farfield patches.

    For each farfield patch p with distance d (wave-patch distance of the cell),
    farfield distance D_p and blending distance b_p:
        main_p    = 1 if d < D_p - b_p, else 0
        blend_p   = 1 if d < D_p, times (1 - main_p), else 0
        scaled_p  = blend_p (1 - (blend_p d - (D_p - b_p)) / b_p)
    s = max over p of (main_p + scaled_p)
    w = (exp(s^k) - 1) / (exp(1) - 1),   k = spatialExponent (default 3.5)

    The weight is 1 in the farfield and 0 away from it, with a blend between.
    It is static (the mesh does not move).

SourceFiles
    swenseRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef swenseRelaxation_H
#define swenseRelaxation_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swenseRelaxation Declaration
\*---------------------------------------------------------------------------*/

class swenseRelaxation
{
    // Private data

        const fvMesh& mesh_;

        //- Relaxation weight, 0 in the interior, 1 in the farfield
        volScalarField w_;

    // Private member functions

        //- Disallow default bitwise copy
        swenseRelaxation(const swenseRelaxation&);
        void operator=(const swenseRelaxation&);


public:

    // Constructors

        //- Construct from mesh and the relaxation dictionary, which holds
        //  farfieldPatchNames, farfieldDistances, blendingDistances and
        //  optionally spatialExponent (default 3.5)
        swenseRelaxation(const fvMesh& mesh, const dictionary& dict);


    // Member Functions

        //- Relaxation weight field
        const volScalarField& weight() const { return w_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
