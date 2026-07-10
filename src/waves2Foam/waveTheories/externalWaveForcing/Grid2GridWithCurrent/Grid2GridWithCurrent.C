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

#include "Grid2GridWithCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "Grid2Grid.H"

namespace Foam
{
namespace waveTheories
{

defineTypeNameAndDebug(Grid2GridWithCurrent, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    Grid2GridWithCurrent,
    externalWaveForcing
);


// * * * * * * * CONSTRUCTOR * * * * * * * //

Grid2GridWithCurrent::Grid2GridWithCurrent
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh),

    waveProps_(io.db().lookupObject<IOdictionary>("waveProperties")),
    coeffDict_(waveProps_.subDict("Grid2GridCoeffs")),

    seaLevel_(readScalar(waveProps_.lookup("seaLevel"))),
    g_(vector::zero),
    rhoWater_(1000),

    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft", 0.0)),

    referenceTime_
    (
        coeffDict_.lookupOrDefault<scalar>("referenceTime", 0.0)
    ),
    referencePosition_
    (
        coeffDict_.lookupOrDefault<vector>("referencePosition", vector::zero)
    ),
    referenceThetaRad_
    (
        coeffDict_.lookupOrDefault<scalar>("referenceTheta", 0.0)*M_PI/180.0
    ),

    Grid2GridInputFileName_(coeffDict_.lookup("Grid2GridInputFileName")),
    transMatrixToOF_(tensor::zero),
    transMatrixFromOF_(tensor::zero),

    // ✅ CURRENT
    Uc_(coeffDict_.lookupOrDefault<vector>("Uc", vector::zero)),
    UF_(coeffDict_.lookupOrDefault<vector>("UF", vector::zero)),
    ZC_(coeffDict_.lookupOrDefault<scalar>("ZC", 0.0)),
    currentTsoft_(coeffDict_.lookupOrDefault<scalar>("currentTsoft", 0.0)),
    profileType_(coeffDict_.lookupOrDefault<word>("profileType", "constant")),
    H_(coeffDict_.lookupOrDefault<scalar>("H", mag(ZC_))),
    z0_(coeffDict_.lookupOrDefault<scalar>("z0", 0.5*(seaLevel_ + ZC_)))
{
    mappingTensors(referenceThetaRad_);

    // ✅ Read density
    IOdictionary transProp
    (
        IOobject
        (
            "transportProperties",
            "constant",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary sD(transProp.subDict(Foam::waves2Foam::waterPhase()));
    rhoWater_ = dimensionedScalar(sD.lookup("rho")).value();

    // ✅ Init Grid2Grid
    string path = mesh_.time().constant()/string(Grid2GridInputFileName_);
    const char* filePath = path.c_str();

    indexGrid2Grid_ = -1;

    Foam::Grid2Grid::__modgrid2grid_MOD_initializegrid2griddict
    (
        filePath,
        (int*)&indexGrid2Grid_
    );

    Info << "Grid2GridWithCurrent initialized (wave + current)" << endl;
}


// * * * * * * PRIVATE FUNCTIONS * * * * * //

scalar Grid2GridWithCurrent::factor(const scalar& time) const
{
    scalar f(1.0);

    if (Tsoft_ > 0.0)
    {
        f = sin(2*M_PI/(4.0*Tsoft_) * min(Tsoft_, time));
    }

    return f;
}


// ✅ Required (copied from Grid2Grid)
double Grid2GridWithCurrent::getGrid2GridTIme(const scalar& time) const
{
    return time + referenceTime_;
}


// ✅ Required (copied from Grid2Grid)
void Grid2GridWithCurrent::getGrid2GridPosition
(
    const point& x,
    double& X,
    double& Y,
    double& Z
) const
{
    vector xx = transMatrixFromOF_ & (x - referencePosition_);

    X = xx.x();
    Y = xx.y();
    Z = xx.z() - seaLevel_;
}


// ✅ CURRENT MODEL
vector Grid2GridWithCurrent::currentU
(
    const point& x,
    const scalar& time
) const
{
    scalar z = x.z();
    vector UC(vector::zero);

    if (profileType_ == "constant")
    {
        UC = Uc_;
    }
    else if (profileType_ == "depthVarying")
    {
        if (z <= ZC_)
        {
            UC = Uc_;
        }
        else
        {
            scalar w = (z - ZC_)/(seaLevel_ - ZC_);
            w = max(min(w, 1.0), 0.0);

            UC = (1.0 - w)*Uc_ + w*UF_;
        }
    }
    else if (profileType_ == "parabolic")
    {
        scalar eta = (z - z0_)/H_;
        scalar Ux  = mag(UF_)*(1.0 - eta*eta);

        UC = vector(max(Ux, 0.0), 0, 0);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown profileType: " << profileType_
            << exit(FatalError);
    }

    scalar fac =
        (currentTsoft_ > 0.0)
      ? sin(M_PI/2.0 * min(time, currentTsoft_) / currentTsoft_)
      : 1.0;

    return UC * fac;
}


// ✅ Mapping (original style, NO new includes)
void Grid2GridWithCurrent::mappingTensors(const scalar& theta)
{
    IOobject gHeader
    (
        "g",
        mesh_.time().constant(),   // correct location
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    uniformDimensionedVectorField gField(gHeader);
    vector g = gField.value();

    g_ = g;
    g /= mag(g);

    if (mag(g.z() + 1.0) > SMALL)
    {
        FatalErrorInFunction
            << "Gravity must be in -Z"
            << exit(FatalError);
    }

    transMatrixFromOF_.xx() = cos(theta);
    transMatrixFromOF_.xy() = -sin(theta);
    transMatrixFromOF_.yx() = sin(theta);
    transMatrixFromOF_.yy() = cos(theta);
    transMatrixFromOF_.zz() = 1;

    transMatrixToOF_ = inv(transMatrixFromOF_);
}


// * * * * * * MAIN FUNCTIONS * * * * * //

void Grid2GridWithCurrent::step()
{
    double time = getGrid2GridTIme( mesh_.time().value() );

    #pragma omp critical(grid2grid_step)
    {
        Foam::Grid2Grid::__modgrid2grid_MOD_correctgrid2grid
        (
            (int*)&indexGrid2Grid_,
            &time
        );
    }
}


scalar Grid2GridWithCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
    double X, Y, Z;
    getGrid2GridPosition(x, X, Y, Z);

    double t = getGrid2GridTIme(time);

    double etaG2G;

    #pragma omp critical(grid2grid_eta)
    {
        Foam::Grid2Grid::__modgrid2grid_MOD_gethoseta
        (
            (int*)&indexGrid2Grid_,
            &X, &Y,
            &t,
            &etaG2G
        );
    }

    scalar eta = etaG2G;

    eta *= factor(time);
    eta += seaLevel_;

    return eta;
}


scalar Grid2GridWithCurrent::pExcess
(
    const point& x,
    const scalar& time
) const
{
    double X, Y, Z;
    getGrid2GridPosition(x, X, Y, Z);

    double t = getGrid2GridTIme(time);

    double pd;

    #pragma omp critical(grid2grid_pressure)
    {
        Foam::Grid2Grid::__modgrid2grid_MOD_gethospd
        (
            (int*)&indexGrid2Grid_,
            &X, &Y, &Z,
            &t,
            &pd
        );
    }

    scalar p = pd;

    p *= factor(time);
    p += rhoWater_ * mag(g_) * seaLevel_;

    return p;
}


// ✅ FINAL VELOCITY
vector Grid2GridWithCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    double X, Y, Z;
    getGrid2GridPosition(x, X, Y, Z);

    double t = getGrid2GridTIme(time);

    double Ux, Uy, Uz;

    #pragma omp critical(grid2grid_velocity)
    {
        Foam::Grid2Grid::__modgrid2grid_MOD_gethosu
        (
            (int*)&indexGrid2Grid_,
            &X, &Y, &Z,
            &t,
            &Ux, &Uy, &Uz
        );
    }

    vector Uw(Ux, Uy, Uz);

    Uw *= factor(time);
    Uw = transMatrixToOF_ & Uw;

    // ✅ CURRENT ADDITION
    vector Ucurr = currentU(x, time);

    return Uw + Ucurr;
}


// * * * * * * INTERFACE * * * * * //

bool Grid2GridWithCurrent::writeData(Ostream&) const
{
    return true;
}

bool Grid2GridWithCurrent::implementPressure() const
{
    return true;
}

void Grid2GridWithCurrent::close()
{
}

} // namespace waveTheories
} // namespace Foam