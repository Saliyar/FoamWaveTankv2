# FoamWaveTankv2

## Overview

**FoamWaveTankv2** is a coupled numerical framework integrating a **Higher-Order Spectral (HOS)** wave model with **OpenFOAM** for high-fidelity wave–structure interaction simulations. The coupling is implemented in a format compatible with **waves2Foam**, enabling efficient transfer of nonlinear wave fields into CFD simulations.

The framework leverages the **Grid2Grid** methodology for consistent interpolation and mapping between the HOS solver and the OpenFOAM computational domain. This allows accurate representation of nonlinear wave kinematics and dynamics within a CFD-based wave tank environment.


---

## Key Features

* Coupling of **HOS wave solver** with **OpenFOAM (.com version)**
* Compatibility with **waves2Foam-style workflows**
* Grid-to-grid interpolation using **Grid2Grid**
* Designed for **nonlinear wave propagation and wave–structure interaction**
* Modular structure for easy extension and testing
* Applicable to offshore and coastal engineering problems

---

## Requirements

* **OpenFOAM v2506 (openfoam.com / ESI version)**. The openfoam.org versions are not supported.
* A C/C++/Fortran compiler, make, CMake, git, FFTW3 and GSL.
* Linux: Ubuntu/Debian (also WSL on Windows), RHEL/Rocky/Fedora, or an HPC server with environment modules.

waves2Foam, Grid2Grid, MoorDyn and FloatStepper are included and compiled by `./Allwmake`.
HOS-NWT ([repository](https://gitlab.com/lheea/HOS-NWT)) is only needed to generate your own HOS wave fields.

---

## Installation

```bash
## Step 1: Load OpenFOAM v2506 in your terminal
source /usr/lib/openfoam/openfoam2506/etc/bashrc     # Ubuntu/WSL package install
# module load openfoam.com/2506                       # HPC server (module names differ per cluster)

## Step 2: Clone the repository
git clone https://github.com/Saliyar/FoamWaveTankv2.git
cd FoamWaveTankv2

## Step 3: Compile (about 10 minutes)
./Allwmake
```

`./Allwmake` first detects where it is running (HPC server, Ubuntu/Debian or RHEL-type Linux) and checks
the prerequisites. If something is missing it compiles nothing, lists what is missing and prints the
command that installs it, for example on Ubuntu:

```bash
sudo apt update && sudo apt install -y build-essential gfortran cmake git libfftw3-dev libgsl-dev
```

(On Ubuntu it offers to run this for you.) If OpenFOAM v2506 is not installed yet, it prints the install
commands as well:

```bash
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt install -y openfoam2506-default
```

Install what it asks for, then run `./Allwmake` again. The build ends with `COMPILATION DONE` and lists
where the solvers and libraries were installed (`$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN`), or with
`COMPILATION FAILED` and the names of the parts that did not build.

---

## Usage

1. Prepare the HOS input and simulation parameters.
2. Set up the OpenFOAM case following the **waves2Foam structure**.
3. Configure coupling parameters (Grid2Grid mapping, time stepping, etc.).
4. Run the simulation:

   ```bash
   ./Allrun
   ```

---

## Applications

This framework is intended for:

* Nonlinear wave propagation studies
* Wave–structure interaction problems
* Floating offshore wind turbine simulations
* Coastal and offshore hydrodynamics research

---

## Current Status

⚠️ **Work in Progress**

This repository is actively under development. Features, structure, and implementation details may change. Validation studies and documentation are ongoing.

![FoamWaveTank methodology and Validation](Sithik_foamU2025.jpg)
---

## Contribution & Feedback

Contributions, suggestions, and feedback are highly welcome.
Please feel free to open an issue or submit a pull request.

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

> M. Hall, A. Goupee (2015). Validation of a lumped-mass mooring line model with DeepCwind semisubmersible
> model test data. *Ocean Engineering*, 104, 590–603.

> H. Chen, M. Hall (2022). CFD simulation of floating body motion with mooring dynamics: Coupling MoorDyn
> with OpenFOAM. *Applied Ocean Research*, 124, 103210. https://doi.org/10.1016/j.apor.2022.103210

**Wave generation and absorption (waves2Foam):**

> N. G. Jacobsen, D. R. Fuhrman, J. Fredsøe (2012). A wave generation toolbox for the open-source CFD
> library: OpenFoam®. *International Journal for Numerical Methods in Fluids*, 70(9), 1073–1088.
> https://doi.org/10.1002/fld.2726

When using HOS wave fields, please also cite [HOS-NWT](https://gitlab.com/lheea/HOS-NWT) and
[Grid2Grid](https://gitlab.com/lheea/Grid2Grid) as described in their repositories.
A dedicated FoamWaveTankv2 paper will be added here when published.

---

## License

FoamWaveTankv2 is free software, released under the **GNU General Public License v3.0 or later**
(see [LICENSE](LICENSE)), the same license family as OpenFOAM and waves2Foam on which it is built.
Third-party components (waves2Foam, Grid2Grid, MoorDyn, FloatStepper) keep their own licenses.

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM
software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.

---

## Author

Sithik Aliyar
DTU Wind and Energy Systems

