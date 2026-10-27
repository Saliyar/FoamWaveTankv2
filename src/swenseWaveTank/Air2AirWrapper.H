/*---------------------------------------------------------------------------*\
  Air2AirWrapper: interface between waveGenerator and Air2AirEnhanced
\*---------------------------------------------------------------------------*/

#ifndef Air2AirWrapper_H
#define Air2AirWrapper_H

#include "Air2AirEnhanced.H"
#include "pointField.H"
#include "vectorField.H"
#include "scalarField.H"
#include "dictionary.H"
#include "autoPtr.H"

namespace Foam
{

class Air2AirWrapper
{
    bool initialized_;
    autoPtr<Air2AirEnhanced> solver_;

    Air2AirWrapper(const Air2AirWrapper&);
    void operator=(const Air2AirWrapper&);

public:

    //- Construct from the Air2Air sub-dictionary of waveProperties
    //  (order, cutoff, cutoffC, cutoffKp, cutoffLength, pressureModel, rhoAir, verbose)
    Air2AirWrapper(const string& modesFilePath, const dictionary& dict);

    ~Air2AirWrapper();

    //- source theory only: see Air2AirEnhanced::setTheory
    void setTheory
    (
        scalar h,
        scalar g,
        scalar c,
        scalar Lp,
        const std::vector<scalar>& eta,
        int M,
        const std::string& cutoff,
        scalar cutoffC
    )
    {
        solver_().setTheory(h, g, c, Lp, eta, M, cutoff, cutoffC);
    }

    //- Set the air coefficients for CFD time t [s]
    void update(scalar t);

    //- z is measured from the mean water level (z - seaLevel)
    vector velocity(scalar x, scalar y, scalar z) const;
    void velocity(const pointField& p, vectorField& result) const;

    scalar pressure(scalar x, scalar y, scalar z) const;
    void pressure(const pointField& p, scalarField& result) const;

    const Air2AirEnhanced& solver() const { return solver_(); }
};

} // End namespace Foam

#endif
