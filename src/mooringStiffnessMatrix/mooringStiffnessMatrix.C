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

\*---------------------------------------------------------------------------*/

#include "mooringStiffnessMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "floaterMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace floaterMotionRestraints
{
    defineTypeNameAndDebug(mooringStiffnessMatrix, 0);

    addToRunTimeSelectionTable
    (
        floaterMotionRestraint,
        mooringStiffnessMatrix,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::SquareMatrix<Foam::scalar>
Foam::floaterMotionRestraints::mooringStiffnessMatrix::readMatrix
(
    const dictionary& dict,
    const word& key,
    const bool mandatory
)
{
    SquareMatrix<scalar> M(6, Zero);

    if (!mandatory && !dict.found(key))
    {
        return M;
    }

    const List<scalarList> rows(dict.get<List<scalarList>>(key));

    if (rows.size() != 6)
    {
        FatalIOErrorInFunction(dict)
            << "Entry '" << key << "' must have 6 rows, found "
            << rows.size() << exit(FatalIOError);
    }

    forAll(rows, i)
    {
        if (rows[i].size() != 6)
        {
            FatalIOErrorInFunction(dict)
                << "Row " << i << " of '" << key << "' must have 6 entries, "
                << "found " << rows[i].size() << exit(FatalIOError);
        }
        for (label j = 0; j < 6; ++j)
        {
            M(i, j) = rows[i][j];
        }
    }

    return M;
}


Foam::vector
Foam::floaterMotionRestraints::mooringStiffnessMatrix::rotationVector
(
    const tensor& R
)
{
    // Axial vector of the skew part: sin(theta)*axis
    const vector s
    (
        0.5*(R.zy() - R.yz()),
        0.5*(R.xz() - R.zx()),
        0.5*(R.yx() - R.xy())
    );

    const scalar sinTheta = mag(s);
    const scalar cosTheta = 0.5*(tr(R) - 1.0);
    const scalar theta = ::atan2(sinTheta, cosTheta);

    if (sinTheta < 1e-12)
    {
        // theta ~ 0 (theta ~ pi is not a mooring use case)
        return s;
    }

    return (theta/sinTheta)*s;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::mooringStiffnessMatrix::mooringStiffnessMatrix
(
    const word& name,
    const dictionary& rBMRDict
)
:
    floaterMotionRestraint(name, rBMRDict),
    refPoint0_(Zero),
    userRefPoint_(false),
    C_(6, Zero),
    B_(6, Zero),
    F0_(Zero),
    report_(true)
{
    read(rBMRDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::floaterMotionRestraints::mooringStiffnessMatrix::restrain
(
    const floaterMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    const point p0 =
        userRefPoint_ ? refPoint0_ : motion.initialCentreOfRotation();

    // Current position of the material point that was at p0
    const point p = motion.transform(p0);

    // Generalised displacement and velocity of the reference point
    const vector dx = p - p0;
    const vector dtheta =
        rotationVector(motion.orientation() & motion.initialQ().T());
    const vector u = motion.velocity(p);
    const vector w = motion.omega();

    FixedList<scalar, 6> xi, xidot, F;
    for (direction k = 0; k < 3; ++k)
    {
        xi[k] = dx[k];
        xi[k+3] = dtheta[k];
        xidot[k] = u[k];
        xidot[k+3] = w[k];
    }

    for (label i = 0; i < 6; ++i)
    {
        scalar Fi = F0_[i];
        for (label j = 0; j < 6; ++j)
        {
            Fi -= C_(i, j)*xi[j] + B_(i, j)*xidot[j];
        }
        F[i] = Fi;
    }

    restraintForce = vector(F[0], F[1], F[2]);
    restraintMoment = vector(F[3], F[4], F[5]);

    // The moment is about p; floaterMotion adds (p - CoR) ^ F
    restraintPosition = p;

    if (report_)
    {
        Info<< "mooringStiffnessMatrix " << name_
            << " t = " << motion.time().value()
            << " xi " << xi
            << " force " << restraintForce
            << " moment " << restraintMoment
            << endl;
    }
}


bool Foam::floaterMotionRestraints::mooringStiffnessMatrix::read
(
    const dictionary& rBMRDict
)
{
    floaterMotionRestraint::read(rBMRDict);

    userRefPoint_ = rBMRCoeffs_.readIfPresent("referencePoint", refPoint0_);

    C_ = readMatrix(rBMRCoeffs_, "stiffness", true);
    B_ = readMatrix(rBMRCoeffs_, "damping", false);

    F0_ = Zero;
    if (rBMRCoeffs_.found("preload"))
    {
        const scalarList f(rBMRCoeffs_.get<scalarList>("preload"));
        if (f.size() != 6)
        {
            FatalIOErrorInFunction(rBMRCoeffs_)
                << "Entry 'preload' must have 6 entries (Fx Fy Fz Mx My Mz), "
                << "found " << f.size() << exit(FatalIOError);
        }
        forAll(f, i)
        {
            F0_[i] = f[i];
        }
    }

    report_ = rBMRCoeffs_.getOrDefault<bool>("report", true);

    return true;
}


void Foam::floaterMotionRestraints::mooringStiffnessMatrix::write
(
    Ostream& os
) const
{
    if (userRefPoint_)
    {
        os.writeEntry("referencePoint", refPoint0_);
    }

    for (const word key : {"stiffness", "damping"})
    {
        const SquareMatrix<scalar>& M = (key == "stiffness") ? C_ : B_;

        os.writeKeyword(key) << nl << token::BEGIN_LIST << nl;
        for (label i = 0; i < 6; ++i)
        {
            os  << "    " << token::BEGIN_LIST;
            for (label j = 0; j < 6; ++j)
            {
                os  << ' ' << M(i, j);
            }
            os  << ' ' << token::END_LIST << nl;
        }
        os  << token::END_LIST << token::END_STATEMENT << nl;
    }

    os.writeKeyword("preload") << F0_ << token::END_STATEMENT << nl;
    os.writeEntry("report", report_);
}


// ************************************************************************* //
