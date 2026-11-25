/*---------------------------------------------------------------------------*\
  Air2AirWrapper implementation
\*---------------------------------------------------------------------------*/

#include "Air2AirWrapper.H"
#include "error.H"
#include "messageStream.H"

namespace Foam
{

Air2AirWrapper::Air2AirWrapper(const string& modesFilePath, const dictionary& dict)
:
    initialized_(false),
    solver_(new Air2AirEnhanced(modesFilePath, dict))
{}

Air2AirWrapper::~Air2AirWrapper()
{}

void Air2AirWrapper::update(scalar t)
{
    solver_->update(t);
    initialized_ = true;
}

vector Air2AirWrapper::velocity(scalar x, scalar y, scalar z) const
{
    if (!initialized_)
    {
        FatalErrorInFunction << "Call update(t) before querying" << exit(FatalError);
    }
    return solver_->velocity(x, y, z);
}

void Air2AirWrapper::velocity(const pointField& p, vectorField& result) const
{
    if (result.size() != p.size()) result.setSize(p.size());
    forAll(p, i)
    {
        result[i] = velocity(p[i].x(), p[i].y(), p[i].z());
    }
}

scalar Air2AirWrapper::pressure(scalar x, scalar y, scalar z) const
{
    if (!initialized_)
    {
        FatalErrorInFunction << "Call update(t) before querying" << exit(FatalError);
    }
    return solver_->pressure(x, y, z);
}

void Air2AirWrapper::pressure(const pointField& p, scalarField& result) const
{
    if (result.size() != p.size()) result.setSize(p.size());
    forAll(p, i)
    {
        result[i] = pressure(p[i].x(), p[i].y(), p[i].z());
    }
}

} // End namespace Foam
