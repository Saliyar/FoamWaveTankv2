// Incident pressure balance, V5 (momentum_swense Part V): log of every step.
//   max |deltaHeld| over the cells
//   max over the cells of |sum_f PhiRIStar_f| (net flux of r_I* per cell)
//   on every body patch: integral of deltaHeld n dS (n = Sf/|Sf|, out of the
//   fluid, the convention of the forces function object, so the pressure
//   force on the body changes by minus this when p uses pIncHeld) and its
//   moment about the origin
{
    const scalarField netStar
    (
        fvc::div(PhiRIStar)().primitiveField()*mesh.V()
    );

    Info<< "SWENSE balance: max|deltaHeld| = "
        << gMax(mag(deltaHeld.primitiveField())())
        << "  max|sum_f PhiRI*| = " << gMax(mag(netStar)())
        << "  max|sum_f PhiRI| = "
        << gMax(mag(fvc::div(PhiRI)().primitiveField()*mesh.V())())
        << endl;

    forAll(bodyPatches, i)
    {
        const label patchi = bodyPatches[i];

        const vectorField dSf
        (
            deltaHeld.boundaryField()[patchi]
           *mesh.Sf().boundaryField()[patchi]
        );

        const vector Fd = gSum(dSf);
        const vector Md = gSum(mesh.Cf().boundaryField()[patchi] ^ dSf);

        Info<< "SWENSE balance body " << mesh.boundaryMesh()[patchi].name()
            << ": int deltaHeld n dS = " << Fd
            << "  moment (origin) = " << Md << endl;
    }
}
