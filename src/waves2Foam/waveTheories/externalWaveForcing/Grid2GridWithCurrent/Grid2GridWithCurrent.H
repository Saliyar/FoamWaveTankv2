#ifndef Grid2GridWithCurrent_H
#define Grid2GridWithCurrent_H

#include "externalWaveForcing.H"

namespace Foam
{
namespace waveTheories
{

class Grid2GridWithCurrent
:
    public externalWaveForcing
{
    // ==============================
    // ORIGINAL GRID2GRID VARIABLES
    // ==============================

    const IOdictionary& waveProps_;
    dictionary coeffDict_;

    scalar seaLevel_;
    vector g_;
    scalar rhoWater_;

    scalar Tsoft_;

    scalar referenceTime_;
    vector referencePosition_;
    scalar referenceThetaRad_;

    word Grid2GridInputFileName_;

    tensor transMatrixToOF_;
    tensor transMatrixFromOF_;

    label indexGrid2Grid_;

    // ==============================
    // ✅ CURRENT VARIABLES
    // ==============================
    vector Uc_;     // deep current
    vector UF_;     // surface current
    scalar ZC_;
    scalar currentTsoft_;
    word profileType_;
    scalar H_;
    scalar z0_;

    // ==============================
    // PRIVATE FUNCTIONS
    // ==============================

    scalar factor(const scalar& time) const;

    void mappingTensors(const scalar& theta);

    // ✅ MUST IMPLEMENT (YOU MISSED THESE)
    double getGrid2GridTIme(const scalar& time) const;

    void getGrid2GridPosition
    (
        const point& x,
        double& X,
        double& Y,
        double& Z
    ) const;

    // ✅ CURRENT MODEL
    vector currentU(const point& x, const scalar& time) const;

public:

    TypeName("Grid2GridWithCurrent");

    Grid2GridWithCurrent
    (
        IOobject io,
        Time& rT,
        const fvMesh& mesh
    );

    void step();
    void close();

    scalar eta(const point& x, const scalar& time) const;

    scalar pExcess(const point& x, const scalar& time) const;

    vector U(const point& x, const scalar& time) const;

    virtual bool writeData(Ostream& os) const;

    virtual bool implementPressure() const;
};

} 
}

#endif
