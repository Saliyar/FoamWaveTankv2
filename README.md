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

If you use this work in your research, please consider citing the associated publication (to be added).

---

## License

(To be specified)

---

## Author

Sithik Aliyar
DTU Wind and Energy Systems

