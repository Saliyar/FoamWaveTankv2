/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of FoamWaveTank, built on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    moorDynStiffnessMatrix

Group
    grpPreProcessingUtilities

Description
    Linearise a MoorDyn v2 mooring system with one coupled body about a body
    position and write the 6x6 stiffness matrix C and the preload F0 for the
    mooringStiffnessMatrix restraint.

    The body state is xi = (x y z roll pitch yaw), as used by the MoorDyn
    BODY coupling (moorDynR2): position of the body reference point and
    rotation angles in rad.  For each DoF j the body is placed at xi0 +/- d e_j,
    the MoorDyn initial-condition solver relaxes the lines, and the static
    force F (Fx Fy Fz Mx My Mz, moment about the reference point) is read
    back.  Then

        F0       = F(xi0)
        C(:, j)  = -(F(xi0 + d e_j) - F(xi0 - d e_j)) / (2 d)

    MoorDyn sets the mass, volume and inertia of a Coupled body to zero, so
    F is the line force only; the CFD body carries its own mass.  The
    linearisation is quasi-static: line inertia and drag are not in C.

    With -sweep, the MoorDyn force and the linear prediction are also written
    for a range of offsets of one DoF, to show where the linear model holds.

Usage
    moorDynStiffnessMatrix -position '(0 0 -0.016 0 0 0)' [OPTIONS]

    -input <file>       MoorDyn input file (default Mooring/lines.txt)
    -position <list>    xi0, 6 values (mandatory)
    -delta <scalar>     translation step in m (default 1e-3)
    -angleDelta <s>     rotation step in rad (default 1e-3)
    -output <file>      include file to write (default constant/mooringStiffness)
    -sweep '(dof amp n)'  write sweep_<dof>.dat with 2n+1 offsets in [-amp, amp]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "FixedList.H"
#include "SquareMatrix.H"
#include "IStringStream.H"

#include "MoorDyn2.h"

using namespace Foam;

typedef FixedList<scalar, 6> vector6;

static const char* dofNames[6] = {"x", "y", "z", "rx", "ry", "rz"};
static const char* forceNames[6] = {"Fx", "Fy", "Fz", "Mx", "My", "Mz"};


// Static force of the mooring system on the coupled body at body state x
vector6 mooringForce(const fileName& input, const vector6& x)
{
    MoorDyn system = MoorDyn_Create(input.c_str());
    if (!system)
    {
        FatalErrorInFunction
            << "MoorDyn could not read " << input << exit(FatalError);
    }
    MoorDyn_SetVerbosity(system, MOORDYN_ERR_LEVEL);
    MoorDyn_SetLogLevel(system, MOORDYN_ERR_LEVEL);

    unsigned int nDof = 0;
    MoorDyn_NCoupledDOF(system, &nDof);
    if (nDof != 6)
    {
        FatalErrorInFunction
            << "Expected one coupled body (6 coupled DoF), MoorDyn reports "
            << nDof << exit(FatalError);
    }

    double X[6], XD[6] = {0, 0, 0, 0, 0, 0}, F[6];
    for (label i = 0; i < 6; ++i)
    {
        X[i] = x[i];
    }

    if (MoorDyn_Init(system, X, XD) != MOORDYN_SUCCESS)
    {
        FatalErrorInFunction
            << "MoorDyn initialisation failed at " << x << exit(FatalError);
    }

    // One mooring time step with the body held still, to get the forces
    double dt = 0;
    MoorDyn_GetDt(system, &dt);
    double t = 0;
    if (MoorDyn_Step(system, X, XD, F, &t, &dt) != MOORDYN_SUCCESS)
    {
        FatalErrorInFunction
            << "MoorDyn step failed at " << x << exit(FatalError);
    }

    MoorDyn_Close(system);

    vector6 f;
    for (label i = 0; i < 6; ++i)
    {
        f[i] = F[i];
    }
    return f;
}


vector6 readVector6(const argList& args, const word& opt)
{
    const scalarList l(args.get<scalarList>(opt));
    if (l.size() != 6)
    {
        FatalErrorInFunction
            << "-" << opt << " needs 6 values, got " << l.size()
            << exit(FatalError);
    }
    vector6 v;
    forAll(l, i)
    {
        v[i] = l[i];
    }
    return v;
}


void writeMatrix(Ostream& os, const word& key, const SquareMatrix<scalar>& M)
{
    os  << key << nl << '(' << nl;
    for (label i = 0; i < 6; ++i)
    {
        os  << "    (";
        for (label j = 0; j < 6; ++j)
        {
            os  << ' ' << M(i, j);
        }
        os  << " )    // " << forceNames[i] << nl;
    }
    os  << ");" << nl;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Derive the 6x6 mooring stiffness matrix and preload from a MoorDyn"
        " input file (one coupled body) by central differences"
    );
    argList::noParallel();
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    argList::addOption("input", "file", "MoorDyn input file");
    argList::addOption("position", "list", "Body state (x y z rx ry rz)");
    argList::addOption("delta", "scalar", "Translation step [m]");
    argList::addOption("angleDelta", "scalar", "Rotation step [rad]");
    argList::addOption("output", "file", "Output include file");
    argList::addOption("sweep", "(dof amp n)", "Linearity sweep of one DoF");

    argList args(argc, argv, false, false, false);

    if (!args.found("position"))
    {
        FatalErrorInFunction
            << "-position '(x y z rx ry rz)' is mandatory" << exit(FatalError);
    }

    const fileName input =
        args.getOrDefault<fileName>("input", "Mooring/lines.txt");
    const fileName output =
        args.getOrDefault<fileName>("output", "constant/mooringStiffness");
    const scalar dLin = args.getOrDefault<scalar>("delta", 1e-3);
    const scalar dAng = args.getOrDefault<scalar>("angleDelta", 1e-3);
    const vector6 x0 = readVector6(args, "position");

    Info<< "MoorDyn input : " << input << nl
        << "xi0           : " << x0 << nl
        << "steps         : " << dLin << " m, " << dAng << " rad" << nl
        << endl;

    const vector6 F0 = mooringForce(input, x0);
    Info<< "F0 = " << F0 << endl;

    SquareMatrix<scalar> C(6, Zero);

    for (label j = 0; j < 6; ++j)
    {
        const scalar d = (j < 3 ? dLin : dAng);
        vector6 xp(x0), xm(x0);
        xp[j] += d;
        xm[j] -= d;

        const vector6 Fp = mooringForce(input, xp);
        const vector6 Fm = mooringForce(input, xm);

        for (label i = 0; i < 6; ++i)
        {
            C(i, j) = -(Fp[i] - Fm[i])/(2*d);
        }
        Info<< "  column " << dofNames[j] << " done" << endl;
    }

    // Remove round-off: entries below 1e-9 of the largest become exact zeros
    {
        scalar Cmax = 0;
        for (label i = 0; i < 6; ++i)
        {
            for (label j = 0; j < 6; ++j)
            {
                Cmax = max(Cmax, mag(C(i, j)));
            }
        }
        for (label i = 0; i < 6; ++i)
        {
            for (label j = 0; j < 6; ++j)
            {
                if (mag(C(i, j)) <= 1e-9*Cmax)
                {
                    C(i, j) = 0;
                }
            }
        }
    }

    Info<< nl << "C =" << nl;
    for (label i = 0; i < 6; ++i)
    {
        Info<< "  ";
        for (label j = 0; j < 6; ++j)
        {
            Info<< ' ' << C(i, j);
        }
        Info<< nl;
    }

    // Asymmetry relative to the largest entry
    scalar maxC = SMALL, maxAsym = 0;
    for (label i = 0; i < 6; ++i)
    {
        for (label j = 0; j < 6; ++j)
        {
            maxC = max(maxC, mag(C(i, j)));
            maxAsym = max(maxAsym, mag(C(i, j) - C(j, i)));
        }
    }
    Info<< nl << "max |C - C^T| / max |C| = " << maxAsym/maxC
        << "  (a static mooring matrix about a point with preload need not be"
        << " symmetric)" << nl << endl;

    {
        mkDir(output.path());
        OFstream os(output);
        os  << "// Mooring stiffness matrix derived by moorDynStiffnessMatrix"
            << nl
            << "// MoorDyn input : " << input << nl
            << "// xi0           : " << x0 << nl
            << "// steps         : " << dLin << " m, " << dAng << " rad"
            << nl
            << "// F = preload - stiffness*(xi - xi0),"
            << " xi = (x y z rx ry rz), F = (Fx Fy Fz Mx My Mz)" << nl << nl;
        os.precision(10);
        os  << "referencePoint  (" << x0[0] << ' ' << x0[1] << ' ' << x0[2]
            << ");" << nl << nl;
        writeMatrix(os, "stiffness", C);
        os  << nl << "preload  " << F0 << ";" << nl;
        Info<< "Wrote " << output << endl;
    }

    if (args.found("sweep"))
    {
        const scalarList s(args.get<scalarList>("sweep"));
        if (s.size() != 3)
        {
            FatalErrorInFunction
                << "-sweep needs (dof amplitude n)" << exit(FatalError);
        }
        const label dof = label(s[0]);
        const scalar amp = s[1];
        const label n = max(label(s[2]), label(1));

        const fileName sweepFile =
            output.path()/("sweep_" + word(dofNames[dof]) + ".dat");
        OFstream os(sweepFile);
        os.precision(10);
        os  << "# offset";
        for (label i = 0; i < 6; ++i)
        {
            os  << " MD_" << forceNames[i];
        }
        for (label i = 0; i < 6; ++i)
        {
            os  << " LIN_" << forceNames[i];
        }
        os  << nl;

        for (label k = -n; k <= n; ++k)
        {
            const scalar off = amp*k/n;
            vector6 x(x0);
            x[dof] += off;
            const vector6 F = mooringForce(input, x);

            os  << off;
            for (label i = 0; i < 6; ++i)
            {
                os  << ' ' << F[i];
            }
            for (label i = 0; i < 6; ++i)
            {
                os  << ' ' << F0[i] - C(i, dof)*off;
            }
            os  << nl;
        }
        Info<< "Wrote " << sweepFile << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
