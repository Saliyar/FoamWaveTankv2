// SWENSE relaxation of alpha toward the incident water fraction, then the
// limited density and mixture update, and the scaled incident pressure.
// Applied after the alpha step of either solver.

alpha1 = (1.0 - w)*alpha1 + w*incident.alphaInc();

{
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1, scalar(0)), scalar(1))
    );
    rho == limitedAlpha1*rho1 + (1.0 - limitedAlpha1)*rho2;
    mixture.correct();
}

// Incident pressure scaled to the mixture density
pStar == rho/rho1*incident.pInc();
