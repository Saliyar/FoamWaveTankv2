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
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "StokesIIIncident.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(StokesIIIncident, 0);
addToRunTimeSelectionTable(waveTheory, StokesIIIncident, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

StokesIIIncident::StokesIIIncident
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H_(0),
    T_(0),
    angle_(0),
    phase_(1.5*M_PI),
    Rlg_(tensor::I),
    Rgl_(tensor::I),
    zBottom_(0),
    h_(0),
    L_(0),
    K_(0),
    omega_(0)
{
    const word waveDictName(coeffDict_.lookup("waveDict"));
    const dictionary& wd = subDict(waveDictName);

    const word esiType(wd.lookup("waveModel"));
    if (esiType != "StokesII")
    {
        FatalIOErrorInFunction(wd)
            << "StokesIIIncident reads the ESI wave dictionary " << waveDictName
            << ", whose waveModel must be StokesII (found " << esiType << ")"
            << exit(FatalIOError);
    }

    H_ = readScalar(wd.lookup("waveHeight"));
    T_ = readScalar(wd.lookup("wavePeriod"));
    angle_ = degToRad(readScalar(wd.lookup("waveAngle")));
    phase_ = wd.lookupOrDefault<scalar>("wavePhase", 1.5*M_PI);

    // Patch frame as ESI waveModel::initialiseGeometry (waveModel.C:57-120):
    // x = -mean patch normal, z = -g, y = z ^ x
    const label patchi = mesh_.boundaryMesh().findPatchID(waveDictName);
    if (patchi < 0)
    {
        FatalErrorInFunction
            << "StokesIIIncident: no patch named " << waveDictName
            << " (the ESI wave dictionary is named after its patch)"
            << exit(FatalError);
    }

    const polyPatch& pp = mesh_.boundaryMesh()[patchi];
    const vector xl = normalised(-gAverage(pp.faceAreas()));
    const vector zl = -g_/mag(g_);
    const vector yl = zl ^ xl;
    Rlg_ = tensor(xl, yl, zl);
    Rgl_ = Rlg_.T();

    const vectorField CpLocal(Rgl_ & pp.localPoints());
    boundBox bb(CpLocal, true);
    zBottom_ = bb.min().z();

    // Depth: seaLevel above the lowest patch point (Rev1 needs seaLevel 0
    // and hydrostatics from z = 0, so the still level is seaLevel)
    h_ = seaLevel_ - zBottom_;

    scalar hRef = 0;
    if (wd.readIfPresent("waterDepthRef", hRef) && mag(hRef - h_) > 1e-9*h_)
    {
        FatalIOErrorInFunction(wd)
            << "waterDepthRef " << hRef << " differs from seaLevel minus the"
            << " lowest point of patch " << waveDictName << " (" << h_ << ")"
            << exit(FatalIOError);
    }

    // Wave length as ESI StokesI::waveLength (StokesIWaveModel.C:74-89)
    const scalar L0 = mag(g_)*T_*T_/(2.0*M_PI);
    L_ = L0;
    for (int i = 1; i <= 100; i++)
    {
        L_ = L0*Foam::tanh(2.0*M_PI*h_/L_);
    }

    K_ = 2.0*M_PI/L_;
    omega_ = 2.0*M_PI/T_;

    Info<< "StokesIIIncident from " << waveDictName << ": H " << H_
        << " m, T " << T_ << " s, angle " << radToDeg(angle_)
        << " deg, phase " << phase_ << ", depth " << h_ << " m (bottom z "
        << zBottom_ << "), L " << L_ << " m, K " << K_ << " 1/m, c "
        << omega_/K_ << " m/s; no ramp" << endl;
}


void StokesIIIncident::printCoeffs()
{
    Info<< "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar StokesIIIncident::phaseTot(const point& xl, const scalar time) const
{
    // StokesIIWaveModel.C:65 with Kx = K cos(angle), Ky = K sin(angle)
    return
        K_*Foam::cos(angle_)*xl.x() + K_*Foam::sin(angle_)*xl.y()
      - omega_*time + phase_;
}


scalar StokesIIIncident::eta
(
    const point& x,
    const scalar& time
) const
{
    const point xl(Rgl_ & x);
    const scalar th = phaseTot(xl, time);
    const scalar sigma = Foam::tanh(K_*h_);

    // StokesIIWaveModel.C:71-73
    const scalar e =
        H_*0.5*Foam::cos(th)
      + K_*H_*H_/4.0*(3.0 - sigma*sigma)/(4.0*pow3(sigma))*Foam::cos(2.0*th);

    // Local z is the vertical (-g); still level = bottom + h = seaLevel
    return zBottom_ + h_ + e;
}


scalar StokesIIIncident::pExcess
(
    const point& x,
    const scalar& time
) const
{
    const point xl(Rgl_ & x);
    const scalar th = phaseTot(xl, time);

    // Height above the bottom (ESI z, z + h of the textbook)
    const scalar Z = xl.z() - zBottom_;
    const scalar rg = rhoWater_*mag(g_);

    // First order
    scalar res =
        rg*H_/2.0*Foam::cosh(K_*Z)/Foam::cosh(K_*h_)*Foam::cos(th);

    // Second order (as stokesSecond.C pExcess, Z + h there = Z here)
    res +=
        1.0/8.0*rg*K_*sqr(H_)/Foam::sinh(2.0*K_*h_)
       *(
            (3.0*Foam::cosh(2.0*K_*Z)/sqr(Foam::sinh(K_*h_)) - 1.0)
           *Foam::cos(2.0*th)
          - Foam::cosh(2.0*K_*Z) + 1.0
        );

    return res + referencePressure();
}


vector StokesIIIncident::U
(
    const point& x,
    const scalar& time
) const
{
    const point xl(Rgl_ & x);
    const scalar th = phaseTot(xl, time);
    const scalar z = xl.z() - zBottom_;
    const scalar k = K_;
    const scalar h = h_;
    const scalar H = H_;
    const scalar omega = omega_;

    // StokesIIWaveModel.C:96-103, z measured from the bottom
    scalar u =
        H*0.5*omega*Foam::cos(th)*Foam::cosh(k*z)/Foam::sinh(k*h)
      + 3.0/4.0*H*H/4.0*omega*k*Foam::cosh(2.0*k*z)/pow4(Foam::sinh(k*h))
       *Foam::cos(2.0*th);
    scalar w =
        H*0.5*omega*Foam::sin(th)*Foam::sinh(k*z)/Foam::sinh(k*h)
      + 3.0/4.0*H*H/4.0*omega*k*Foam::sinh(2.0*k*z)/pow4(Foam::sinh(k*h))
       *Foam::sin(2.0*th);
    const scalar v = u*Foam::sin(angle_);
    u *= Foam::cos(angle_);

    // Local to global, as waveModel.C:402
    return Rlg_ & vector(u, v, w);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
