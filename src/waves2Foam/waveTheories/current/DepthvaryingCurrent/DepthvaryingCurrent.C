/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/



#include "DepthvaryingCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace waveTheories
{

defineTypeNameAndDebug(DepthvaryingCurrent, 0);
addToRunTimeSelectionTable(waveTheory, DepthvaryingCurrent, dictionary);


// Constructors
DepthvaryingCurrent::DepthvaryingCurrent
(
    const word& subDictName,
    const fvMesh& mesh
)
:
    waveTheory(subDictName, mesh),
    U_(vector(coeffDict_.lookup("U"))),
    UF_(vector(coeffDict_.lookup("UF"))),
    ZC_(readScalar(coeffDict_.lookup("ZC"))),
    Tsoft_(readScalar(coeffDict_.lookup("Tsoft"))),
    localSeaLevel_
    (
        coeffDict_.lookupOrDefault<scalar>("localSeaLevel", seaLevel_)
    ),
    profileType_(coeffDict_.lookup("profileType")),
    H_(coeffDict_.lookupOrDefault<scalar>("H", mag(ZC_))),
    z0_(coeffDict_.lookupOrDefault<scalar>
        (
            "z0",
            0.5*(localSeaLevel_ + ZC_)
        )
    )
{
    // ---- Validate profileType immediately -----------------

    if
    (
        profileType_ != "constant"
     && profileType_ != "depthVarying"
     && profileType_ != "parabolic"
    )
    {
        FatalIOErrorInFunction(coeffDict_)
            << "Invalid profileType \"" << profileType_ << "\"" << nl
            << "Valid profileType options are:" << nl
            << "  - constant" << nl
            << "  - depthVarying" << nl
            << "  - parabolic" << nl
            << exit(FatalIOError);
    }
}


// Print coefficients
void DepthvaryingCurrent::printCoeffs()
{
    Info<< "Loading current model: " << typeName << nl
        << "  Profile type        : " << profileType_ << nl
        << "  Deep current U      : " << U_ << nl
        << "  Surface current UF : " << UF_ << nl
        << "  Transition depth ZC: " << ZC_ << nl
        << "  Soft-start Tsoft   : " << Tsoft_ << nl << endl;
}


// Soft start
scalar DepthvaryingCurrent::factor(const scalar& time) const
{
    if (Tsoft_ <= 0.0) return 1.0;

    return sin(constant::mathematical::pi/2.0
         * min(time, Tsoft_) / Tsoft_);
}


// Free surface
scalar DepthvaryingCurrent::eta
(
    const point&,
    const scalar&
) const
{
    return localSeaLevel_;
}


// Excess pressure
scalar DepthvaryingCurrent::pExcess
(
    const point&,
    const scalar&
) const
{
    return referencePressure(localSeaLevel_);
}


// Velocity profile (CRITICAL FUNCTION)
vector DepthvaryingCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    const scalar z = x.z();
    vector UC(Zero);

    // ---------------- CONSTANT CURRENT ----------------
    if (profileType_ == "constant")
    {
        UC = U_;
    }

    // ---------------- DEPTH-VARYING (CORRECTED) --------
    else if (profileType_ == "depthVarying")
    {
        if (z <= ZC_)
        {
            UC = U_;
        }
        else
        {
            scalar w = (z - ZC_)/(localSeaLevel_ - ZC_);
            w = max(min(w, 1.0), 0.0);

            UC = (1.0 - w)*U_ + w*UF_;
        }
    }

    // ---------------- PARABOLIC PROFILE ----------------
    else if (profileType_ == "parabolic")
    {
        scalar eta = (z - z0_)/H_;
        scalar Ux  = mag(UF_)*(1.0 - eta*eta);

        UC = vector(max(Ux, 0.0), 0, 0);
    }

    else
    {
        FatalErrorInFunction
            << "Unknown profileType: " << profileType_ << nl
            << "Valid types are: constant | depthVarying | parabolic"
            << exit(FatalError);
    }

    return UC * factor(time);
}


Switch DepthvaryingCurrent::conflictTSoftInitialise() const
{
    return false;
}

} // End namespace waveTheories
} // End namespace Foam


