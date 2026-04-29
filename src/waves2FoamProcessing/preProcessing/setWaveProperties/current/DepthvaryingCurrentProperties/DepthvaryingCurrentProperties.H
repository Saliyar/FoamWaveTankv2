/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::DepthvaryingCurrentProperties

Description
    Property wrapper for the DepthvaryingCurrent wave theory.

    This class does not implement any physics. It provides the dictionary
    interface required by waves2Foam to construct the DepthvaryingCurrent
    wave theory at run-time.

\*---------------------------------------------------------------------------*/

#ifndef DepthvaryingCurrentProperties_H
#define DepthvaryingCurrentProperties_H

#include "setWaveProperties.H"

namespace Foam
{

class DepthvaryingCurrentProperties
:
    public setWaveProperties
{
public:

    //- Runtime type information
    TypeName("DepthvaryingCurrentProperties");

    // Constructors
    DepthvaryingCurrentProperties
    (
        const Time&,
        dictionary&,
        vector,
        bool
    );

    //- Write wave properties
    virtual void set(Ostream&);
};

} // End namespace Foam

#endif