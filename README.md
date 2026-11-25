# FoamWaveTank

A numerical wave tank for **OpenFOAM v2506** (openfoam.com), built on
[waves2Foam](https://doi.org/10.1002/fld.2726). It imposes nonlinear, phase-resolved wave fields from
the Higher-Order Spectral (HOS) solvers HOS-NWT and HOS-Ocean in the CFD domain, adds ocean currents,
and simulates floating bodies with [FloatStepper](https://doi.org/10.1098/rsos.231453) and moorings with
[MoorDyn](https://github.com/FloatingArrayDesign/MoorDyn).

Developed and maintained by Sithik Aliyar, DTU Wind and Energy Systems.
Please cite the papers listed under [Citation](#citation) when you use the code.

**Documentation:** [FoamWaveTank Manual](doc/manual/FoamWaveTank_Manual.pdf) (what FoamWaveTank adds);
the waves2Foam Manual (N. G. Jacobsen, 2017) for the relaxation zones, wave theories and utilities it inherits.

![FoamWaveTank methodology and validation](Sithik_foamU2025.jpg)

---

## The solvers

| Solver | Interface capturing | Wave zones | Body motion | Manual |
|---|---|---|---|---|
| `waveFoam` | MULES | yes | `sixDoFRigidBodyMotion` | Sec. 4.1 |
| `waveIsoFoam` | isoAdvector | yes | `sixDoFRigidBodyMotion` | Sec. 4.2 |
| `floatStepperWaveIsoFoam` | isoAdvector | yes | FloatStepper | Sec. 4.4 |
| `floatStepper` | isoAdvector | no | FloatStepper | Sec. 4.4 |

All solvers accept every waves2Foam wave theory and the HOS forcing, and have a `ddtCorr on|off;`
switch in `system/fvSchemes` (default `on`, manual Sec. 4.3).

What FoamWaveTank adds to waves2Foam:

| Feature | Keyword / library | Manual |
|---|---|---|
| HOS wave field (HOS-NWT, HOS-Ocean via Grid2Grid) | `externalForcing Grid2Grid;` | Ch. 3 |
| HOS waves plus a current (constant, depth-varying, parabolic) | `externalForcing Grid2GridWithCurrent;` | Sec. 3.5 |
| Current only, in a relaxation zone | `waveType DepthvaryingCurrent;` | Sec. 3.6 |
| MoorDyn v1/v2 moorings for FloatStepper | `floaterMotionRestraint moorDynR1/moorDynR2` | Ch. 5 |
| MoorDyn v1/v2 moorings for `sixDoFRigidBodyMotion` | `sixDoFRigidBodyMotionRestraint moorDynR1/moorDynR2` | Sec. 5.3 |
| Linear 6x6 mooring stiffness matrix for FloatStepper, derived from a MoorDyn file, for comparison with potential-flow/BEM codes | `floaterMotionRestraint mooringStiffnessMatrix`, utility `moorDynStiffnessMatrix` | -- |
| Installation with prerequisite and result checks | `./Allwmake` | Ch. 2 |

---

## Contents

| Directory | Contents |
|---|---|
| `src/waves2Foam/` | waves2Foam library, including the HOS (`externalWaveForcing/Grid2Grid*`) and current (`current/DepthvaryingCurrent`) theories |
| `src/moorDynSixDoF/` | MoorDyn restraints for `sixDoFRigidBodyMotion` (`libmoorDynR1SixDoF`, `libmoorDynR2SixDoF`) |
| `src/mooringStiffnessMatrix/` | linear 6x6 mooring restraint for FloatStepper (`libmooringStiffnessMatrix`) |
| `src/waves2Foam*/` | other waves2Foam libraries (processing, sampling, GABC) |
| `applications/solvers/` | `waveFoam`, `waveIsoFoam`, and `floatStepperWaveIsoFoam/` (FloatStepper, git submodule, with MoorDyn) |
| `applications/utilities/` | waves2Foam pre- and post-processing utilities (`setWaveField`, `setWaveParameters`, `relaxationZoneLayout`, `waveGaugesNProbes`, `moorDynStiffnessMatrix`, ...) |
| `ThirdParty/` | Grid2Grid (HOS wrapper, downloaded and built by `Allwmake`) and fenton4Foam |
| `tutorials/` | example cases (below) |
| `doc/manual/` | the manual, LaTeX source and PDF |
| `bin/` | environment file (`bashrc`, created from `bashrc.org`) |

---

## Installation

```bash
## Step 1: Load OpenFOAM v2506 in your terminal
source /usr/lib/openfoam/openfoam2506/etc/bashrc     # Ubuntu/WSL package install
# module load openfoam.com/2506                       # HPC server (module names differ per cluster)

## Step 2: Clone the repository
git clone https://github.com/Saliyar/FoamWaveTankv2.git
cd FoamWaveTankv2

## Step 3: Compile (10-20 minutes)
./Allwmake
```

`./Allwmake` first detects where it is running (HPC server, Ubuntu/Debian or RHEL-type Linux) and checks
the prerequisites: compilers, make, CMake, git, FFTW, GSL and OpenFOAM v2506. If something is missing it
compiles nothing, lists what is missing and prints the command that installs it, for example on Ubuntu:

```bash
sudo apt update && sudo apt install -y build-essential gfortran cmake git libfftw3-dev libgsl-dev
```

(On Ubuntu it offers to run this for you.) If OpenFOAM v2506 is not installed, it prints the install
commands as well:

```bash
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt install -y openfoam2506-default
```

The build ends with `COMPILATION DONE` and the folders of the solvers and libraries
(`$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN`), or with `COMPILATION FAILED` and the parts that did not build.
The openfoam.org versions of OpenFOAM are not supported.

## Updating

```bash
cd FoamWaveTankv2
git pull
./Allwmake
```

`Allwmake` also updates FloatStepper. Clones made before 23 November 2026 have the FloatStepper submodule
in the old folder `applications/solvers/floatStepperWaveFoam`; after updating, delete it with
`rm -rf applications/solvers/floatStepperWaveFoam`.

---

## Tutorials

Copy a tutorial to a working folder and run `./Allrun` there (`./Allclean` resets it).

| Tutorial | Solver | Shows | Needs |
|---|---|---|---|
| `tutorials/waveFoam` | `waveFoam` | HOS-NWT waves through inlet and outlet relaxation zones | `constant/modes_HOS_SWENSE.dat` |
| `tutorials/waveIsoFoam` | `waveIsoFoam` | the same case with isoAdvector | `constant/modes_HOS_SWENSE.dat` |
| `tutorials/HOS_with_Current` | `waveFoam` | HOS waves with a depth-varying current | `constant/modes_HOS_SWENSE.dat` |
| `tutorials/CurrentWithProfile` | `waveIsoFoam` | a current without waves (constant, depth-varying or parabolic profile) | -- |
| `tutorials/floatStepperWaveIsoFoam` | `floatStepperWaveIsoFoam` | a box moored with MoorDyn v2 in regular waves (about 5 min on one core) | -- |
| `tutorials/floatStepperWaveIsoFoam_stiffnessMatrix` | `floatStepperWaveIsoFoam` | the same box with a linear 6x6 mooring stiffness matrix derived from the MoorDyn file; `Allrun` derives C, runs the C-matrix and MoorDyn cases side by side and plots the comparison (about 20 min on two cores); see its `README.md` | -- |
| `tutorials/postProcessingScript` | MATLAB | wave comparison, spectra, peak statistics | -- |

The HOS file is produced by HOS-NWT ([repository](https://gitlab.com/lheea/HOS-NWT)) and is not part of
this repository.

---

## Citation

If you use FoamWaveTankv2, please cite (GitHub also shows this under **"Cite this repository"**, from
[CITATION.cff](CITATION.cff)):

**HOS–OpenFOAM coupling (FoamWaveTank):**

> S. Aliyar, G. Ducrozet, B. Bouscasse, F. Bonnefoy, V. Sriram, P. Ferrant (2022). Numerical coupling
> strategy using HOS-OpenFOAM-MoorDyn for OC3 Hywind SPAR type platform. *Ocean Engineering*, 263, 112206.
> https://doi.org/10.1016/j.oceaneng.2022.112206

**Floating bodies with FloatStepper** (`floatStepper`, `floatStepperWaveIsoFoam`):

> J. Roenby, S. Aliyar, H. Bredmose (2024). A robust algorithm for computational floating body dynamics.
> *Royal Society Open Science*, 11, 231453. https://doi.org/10.1098/rsos.231453

**Moorings with MoorDyn:**

> S. Aliyar, G. Ducrozet, B. Bouscasse, F. Bonnefoy, V. Sriram, P. Ferrant (2022). Numerical coupling
> strategy using HOS-OpenFOAM-MoorDyn for OC3 Hywind SPAR type platform. *Ocean Engineering*, 263, 112206.
> https://doi.org/10.1016/j.oceaneng.2022.112206

> H. Chen, M. Hall (2022). CFD simulation of floating body motion with mooring dynamics: Coupling MoorDyn
> with OpenFOAM. *Applied Ocean Research*, 124, 103210. https://doi.org/10.1016/j.apor.2022.103210

> M. Hall, A. Goupee (2015). Validation of a lumped-mass mooring line model with DeepCwind semisubmersible
> model test data. *Ocean Engineering*, 104, 590–603.

**Wave generation and absorption (waves2Foam):**

> N. G. Jacobsen, D. R. Fuhrman, J. Fredsøe (2012). A wave generation toolbox for the open-source CFD
> library: OpenFoam®. *International Journal for Numerical Methods in Fluids*, 70(9), 1073–1088.
> https://doi.org/10.1002/fld.2726

When using HOS wave fields, please also cite [HOS-NWT](https://gitlab.com/lheea/HOS-NWT) and
[Grid2Grid](https://gitlab.com/lheea/Grid2Grid) as described in their repositories.
A dedicated FoamWaveTankv2 paper will be added here when published.

---

## Use of generative AI

Claude Code, with Anthropic's Claude Opus models, was used in the development of this repository,
including the installation script (`Allwmake`), the MoorDyn restraints for `sixDoFRigidBodyMotion`,
bug fixes in the solvers and utilities, the tutorial fixes, the build tests, the manual and this README.
All code, cases and documentation produced with the tool have been reviewed by the author, who takes
full responsibility for them.

---

## License

FoamWaveTankv2 is free software, released under the **GNU General Public License v3.0 or later**
(see [LICENSE](LICENSE)), the same license family as OpenFOAM and waves2Foam on which it is built.
Third-party components (waves2Foam, Grid2Grid, MoorDyn, FloatStepper) keep their own licenses.

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM
software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.

---

## Contact and bug reports

Please report bugs and questions on the issue tracker of this repository.

Sithik Aliyar, DTU Wind and Energy Systems, Technical University of Denmark.
