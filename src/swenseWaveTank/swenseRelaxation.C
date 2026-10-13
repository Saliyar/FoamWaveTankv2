/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseRelaxation.H"
#include "patchWave.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * Constructors * * * * * * * * //

Foam::swenseRelaxation::swenseRelaxation
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    w_
    (
        IOobject("swenseWeight", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("swenseWeight", dimless, 0.0)
    )
{
    const wordList patchNames(dict.lookup("farfieldPatchNames"));
    const scalarList farDist(dict.lookup("farfieldDistances"));
    const scalarList blendDist(dict.lookup("blendingDistances"));
    const scalar spatialExp
    (
        dict.lookupOrDefault<scalar>("spatialExponent", 3.5)
    );

    if (patchNames.size() != farDist.size() || farDist.size() != blendDist.size())
    {
        FatalErrorInFunction
            << "farfieldPatchNames, farfieldDistances and blendingDistances "
            << "must have the same size" << exit(FatalError);
    }

    scalarField s(mesh.nCells(), 0.0);

    forAll(patchNames, patchI)
    {
        const label patchID = mesh.boundaryMesh().findPatchID(patchNames[patchI]);

        if (patchID < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patchNames[patchI] << nl
                << "Available patch names: " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        labelHashSet patchIDs;
        patchIDs.insert(patchID);

        // Distance of each cell from this patch (not corrected at walls)
        patchWave waveDistance(mesh, patchIDs, false);
        const scalarField& d = waveDistance.distance();

        const scalar Dff = farDist[patchI];
        const scalar b = blendDist[patchI];

        forAll(s, celli)
        {
            // main: 1 in the farfield away from the blend band
            const scalar main = pos((Dff - b) - d[celli]);

            // blend: 1 in the band, then linear
            const scalar blend = pos(Dff - d[celli])*(1.0 - main);

            const scalar scaled =
                blend*(1.0 - (blend*d[celli] - (Dff - b))/b);

            s[celli] = max(s[celli], main + scaled);
        }
    }

    forAll(s, celli)
    {
        w_[celli] =
            (Foam::exp(Foam::pow(s[celli], spatialExp)) - 1.0)
           /(Foam::exp(1.0) - 1.0);
    }
}

// ************************************************************************* //
