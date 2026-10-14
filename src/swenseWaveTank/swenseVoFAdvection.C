/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "swenseVoFAdvection.H"
#include "fvCFD.H"
#include "MULES.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "gaussConvectionScheme.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::swenseRho
(
    const volScalarField& s,
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2
)
{
    return rho2 + (rho1 - rho2)*min(max(s, scalar(0)), scalar(1));
}


void Foam::swenseAdvanceVoF
(
    volScalarField& alpha,
    const surfaceScalarField& phi,
    const volVectorField& U,
    surfaceScalarField& rhoPhi,
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2
)
{
    const fvMesh& mesh = alpha.mesh();
    const dictionary& dict = mesh.solverDict(alpha.name());

    const label nAlphaCorr(readLabel(dict.lookup("nAlphaCorr")));
    const bool MULESCorr(dict.lookupOrDefault<Switch>("MULESCorr", false));
    const scalar cAlpha(dict.lookupOrDefault<scalar>("cAlpha", 0.0));
    const scalar icAlpha(dict.lookupOrDefault<scalar>("icAlpha", 0.0));
    const scalar scAlpha(dict.lookupOrDefault<scalar>("scAlpha", 0.0));

    if (dict.lookupOrDefault<label>("nAlphaSubCycles", 1) != 1)
    {
        FatalErrorInFunction
            << "nAlphaSubCycles must be 1 for the SWENSE volume fractions"
            << exit(FatalError);
    }

    if (dict.lookupOrDefault<Switch>("alphaApplyPrevCorr", false))
    {
        FatalErrorInFunction
            << "alphaApplyPrevCorr is not supported for the SWENSE volume fractions"
            << exit(FatalError);
    }

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    // Advective form (document Rev2, 16.7): d(alpha)/dt + div(alpha phi)
    // = alpha div(phi), with the discrete divergence of the same face flux.
    // Source form for MULES: Sp = div(phi) (implicit in the predictor and in
    // the explicit update). In the limited correction the source is already
    // in psi, so Su = -Sp psi, as OpenFOAM alphaEqn.H does; otherwise the
    // source would be applied twice.
    const volScalarField divPhi("divPhi", fvc::div(phi));
    zeroField Su;

    // Restart every call from the old-time value (document Section 16.3.2)
    alpha = alpha.oldTime();

    // Compression coefficient, with the optional terms of alphaEqn.H
    surfaceScalarField phic(cAlpha*mag(phi/mesh.magSf()));

    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (cAlpha*icAlpha)*fvc::interpolate(mag(U));
    }

    if (scAlpha > 0)
    {
        phic +=
            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }

    // No compression at non-coupled boundary faces (inlets, outlets etc.)
    surfaceScalarField::Boundary& phicBf = phic.boundaryFieldRef();

    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    );

    // MULES flux from the upwind predictor (MULESCorr)
    surfaceScalarField alphaPhi10
    (
        "alphaPhi10",
        phi*fvc::interpolate(alpha)
    );

    if (MULESCorr)
    {
        fvScalarMatrix alphaEqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha)
          == fvm::Sp(divPhi, alpha)
        );

        alphaEqn.solve(dict);

        alphaPhi10 = phi*upwind<scalar>(mesh, phi).interpolate(alpha);
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        const volScalarField alpha2("alpha2", 1.0 - alpha);

        // Interface normal from the gradient of this alpha
        const surfaceVectorField gradAlphaf(fvc::interpolate(fvc::grad(alpha)));
        const surfaceScalarField nHatf
        (
            (gradAlphaf/(mag(gradAlphaf) + deltaN)) & mesh.Sf()
        );

        const surfaceScalarField phir(phic*nHatf);

        surfaceScalarField alphaPhiUn
        (
            fvc::flux(phi, alpha, alphaScheme)
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(alphaPhiUn - alphaPhi10);
            volScalarField alpha0("alpha0", alpha);

            MULES::correct
            (
                geometricOneField(),
                alpha,
                alphaPhiUn,
                talphaPhiCorr.ref(),
                divPhi,
                (-divPhi*alpha)(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhiCorr();
            }
            else
            {
                alpha = 0.5*alpha + 0.5*alpha0;
                alphaPhi10 += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            alphaPhi10 = alphaPhiUn;

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha,
                phi,
                alphaPhi10,
                divPhi,
                Su,
                oneField(),
                zeroField()
            );
        }
    }

    rhoPhi = alphaPhi10*(rho1 - rho2) + phi*rho2;
}

// ************************************************************************* //
