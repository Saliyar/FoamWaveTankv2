/*---------------------------------------------------------------------------*\
Description
    Volume-fraction transport and density for the SWENSE solvers (document
    Rev1, Section 16.1 and 16.3).

    swenseRho(s)
        rho(s) = rho2 + (rho1 - rho2) clamp(s, 0, 1)        (eq. 16.1)
    The same function builds rho from alpha and rhoI from alphaI.

    swenseAdvanceVoF
        Advances a volume fraction from its old-time value with a given volume
        flux, with the MULES settings of its solver dictionary (the same
        settings for every field), and returns its mass flux
        rhoPhi = alphaPhi (rho1 - rho2) + phi rho2. It is called for alpha
        with the total flux phi and for alphaI with the incident flux phiI.

SourceFiles
    swenseVoFAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef swenseVoFAdvection_H
#define swenseVoFAdvection_H

#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Density of a fluid with volume fraction s of the first phase
tmp<volScalarField> swenseRho
(
    const volScalarField& s,
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2
);

//- Advance the volume fraction alpha from its old-time value with the flux
//  phi (VoF scheme of OpenFOAM alphaEqn.H, Euler time). Returns the mass flux
void swenseAdvanceVoF
(
    volScalarField& alpha,
    const surfaceScalarField& phi,
    const volVectorField& U,
    surfaceScalarField& rhoPhi,
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
