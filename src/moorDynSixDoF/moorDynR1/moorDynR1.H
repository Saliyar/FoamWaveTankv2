/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::sixDoFRigidBodyMotionRestraints::moorDynR1

Description
    sixDoFRigidBodyMotionRestraints model.  Mooring dynamics resolved.

    Port of the FloatStepper floaterMotionRestraints::moorDynR1 to
    sixDoFRigidBodyMotion.  Physics unchanged: X = CoM - initialCoM,
    Euler XYZ angles, v, omega.

    MoorDyn v1 cannot save/restore its state, so LinesCalc is called at most
    once per time step.  Repeated calls with the same time (outer correctors
    with moveMeshOuterCorrectors) return the force of the first call.

    Calling MoorDyn: a lumped-mass based mooring dynamics model by
    M. Hall and A. Goupee, “Validation of a lumped-mass mooring line model with
    DeepCwind semisubmersible model test data,” Ocean Engineering, vol. 104,
    pp. 590–603, Aug. 2015.

    Chen, H., & Hall, M. (2022). CFD simulation of floating body motion with
    mooring dynamics: Coupling MoorDyn with OpenFOAM. Applied Ocean Research,
    124, 103210. https://doi.org/10.1016/j.apor.2022.103210

SourceFiles
    moorDynR1.C

\*---------------------------------------------------------------------------*/

#ifndef moorDynR1_H
#define moorDynR1_H

#include "sixDoFRigidBodyMotionRestraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class moorDynR1 Declaration
\*---------------------------------------------------------------------------*/

class moorDynR1
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- Initialize MoorDyn at the first time step ()
        mutable bool initialized_;

        //- Time of the last LinesCalc call
        mutable scalar calcTime_;

        //- Number of cached returns within the current time step
        mutable label nCacheUse_;

        //- Force and moment of the last LinesCalc call
        mutable vector cachedForce_;
        mutable vector cachedMoment_;

        //- Allow a restart (MoorDyn v1 cannot save its state, so the lines
        //  restart from static equilibrium); stop otherwise (default)
        bool allowStaticRestart_;


public:

    //- Runtime type information
    TypeName("moorDynR1");


    // Constructors

        //- Construct from components
        moorDynR1
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new moorDynR1(*this)
            );
        }


    //- Destructor
    virtual ~moorDynR1();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sixDoFRigidBodyMotionRestraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
