/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of FoamWaveTank, built on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::floaterMotionRestraints::mooringStiffnessMatrix

Description
    Linear 6-DoF mooring restraint for FloatStepper (floaterMotion):

        F = F0 - C (xi - xi0) - B dxi/dt

    with xi = (x y z rx ry rz), the translation of the reference point and the
    rotation of the body (rotation vector of Q Q0^T, which equals the
    roll/pitch/yaw angles for small rotations).  F = (Fx Fy Fz Mx My Mz) is
    the force and the moment about the current position of the reference
    point, in global axes.

    It is an alternative to the MoorDyn restraints (moorDynR1, moorDynR2),
    meant for code-to-code comparison with potential-flow and BEM tools that
    model the mooring as a 6x6 stiffness matrix.  The matrix can be derived
    from a MoorDyn input file with the moorDynStiffnessMatrix utility, which
    writes an include file with C and F0.

    Units: C rows/columns are SI, i.e. N/m, N/rad, N m/m and N m/rad.

Usage
    \verbatim
    restraints
    {
        mooring
        {
            floaterMotionRestraint  mooringStiffnessMatrix;

            // Optional. Default: the body's initial centre of rotation
            referencePoint  (0 0 -0.016);

            // 6x6, row i = force/moment component i
            stiffness
            (
                (  C11 C12 C13 C14 C15 C16 )
                ...
            );

            // Optional. Force/moment at xi = xi0. Default: zero
            preload  (Fx Fy Fz Mx My Mz);

            // Optional 6x6 linear damping. Default: zero
            // damping  ( (...) ... );

            // Optional. Print xi and F every call. Default: true
            report   true;
        }
    }
    \endverbatim

    The output of moorDynStiffnessMatrix can be included directly:
    \verbatim
        mooring
        {
            floaterMotionRestraint  mooringStiffnessMatrix;
            #include "$FOAM_CASE/constant/mooringStiffness"
        }
    \endverbatim

SourceFiles
    mooringStiffnessMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef mooringStiffnessMatrix_H
#define mooringStiffnessMatrix_H

#include "floaterMotionRestraint.H"
#include "SquareMatrix.H"
#include "point.H"
#include "tensor.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace floaterMotionRestraints
{

/*---------------------------------------------------------------------------*\
                   Class mooringStiffnessMatrix Declaration
\*---------------------------------------------------------------------------*/

class mooringStiffnessMatrix
:
    public floaterMotionRestraint
{
    // Private Data

        //- Reference point in the initial configuration. Unset: use the
        //  body's initial centre of rotation
        point refPoint0_;

        //- True if referencePoint was given
        bool userRefPoint_;

        //- Stiffness matrix C (6x6)
        SquareMatrix<scalar> C_;

        //- Damping matrix B (6x6)
        SquareMatrix<scalar> B_;

        //- Force and moment at the reference configuration
        FixedList<scalar, 6> F0_;

        //- Print the state every call
        bool report_;


    // Private Member Functions

        //- Read a 6x6 matrix given as a list of 6 lists of 6 scalars
        static SquareMatrix<scalar> readMatrix
        (
            const dictionary& dict,
            const word& key,
            const bool mandatory
        );

        //- Rotation vector (axis * angle) of a rotation tensor
        static vector rotationVector(const tensor& R);


public:

    //- Runtime type information
    TypeName("mooringStiffnessMatrix");


    // Constructors

        //- Construct from components
        mooringStiffnessMatrix
        (
            const word& name,
            const dictionary& rBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<floaterMotionRestraint> clone() const
        {
            return autoPtr<floaterMotionRestraint>
            (
                new mooringStiffnessMatrix(*this)
            );
        }


    //- Destructor
    virtual ~mooringStiffnessMatrix() = default;


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const floaterMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& rBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace floaterMotionRestraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
