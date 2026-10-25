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

* **OpenFOAM (.com version)** (validated for versions ≥ v2024)
* **waves2Foam framework**
* **Grid2Grid library/toolbox**
* Linux-based environment (recommended)

---

## Installation

You need OpenFOAM (.com), a C/C++ and Fortran compiler, CMake and GSL. Grid2Grid comes with the repository as a git submodule, and FFTW (which Grid2Grid needs) is found or built automatically, so there are no paths to edit.

```bash
## Step 0: Prerequisites
# Ubuntu/Debian:  sudo apt install build-essential gfortran cmake libgsl-dev
# RHEL/Rocky:     sudo dnf install gcc gcc-c++ gcc-gfortran cmake gsl-devel
# Cluster:        module load the OpenFOAM, gcc, cmake and gsl modules
#                 (an FFTW module is optional, see below)

## Step 1: Load OpenFOAM, e.g.
source /usr/lib/openfoam/openfoam2412/etc/bashrc

## Step 2: Clone with submodules
git clone --recurse-submodules https://github.com/Saliyar/FoamWaveTankv2.git
cd FoamWaveTankv2

## Step 3: Compile (uses all cores; ./Allwmake -j 4 to limit)
./Allwmake
```

The build stops early with a clear message if something is missing. The full output is saved in `log.build`.

**How FFTW is handled.** `ThirdParty/Allwmake` checks, in order, `FFTW3_LIB_PATH` (if you set it), a previous local build in `ThirdParty/fftw`, loaded modules and the standard system folders. Every `libfftw3.a` it finds is test-linked into a shared library first, because a static FFTW built without `-fPIC` cannot be linked into `libGrid2Grid.so`. If none works, FFTW 3.3.10 is downloaded and built once into `ThirdParty/fftw` (a few minutes).

**No internet on the build node?** Download https://www.fftw.org/fftw-3.3.10.tar.gz elsewhere and copy it into `ThirdParty/` before running `./Allwmake`.

**Rebuilding Grid2Grid.** It is skipped when `$FOAM_USER_LIBBIN/libGrid2Grid.so` already exists. Run `ThirdParty/Allclean` to force a rebuild (`ThirdParty/Allclean -all` also removes the local FFTW).

Related tools: [HOS-NWT](https://gitlab.com/lheea/HOS-NWT) generates the wave fields, and [Grid2Grid](https://gitlab.com/lheea/Grid2Grid) maps them onto the OpenFOAM mesh.

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

