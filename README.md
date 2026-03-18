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
> **Note:** Before installing FoamWaveTankv2, make sure you have the following dependencies installed:
> 
> - **HOS-NWT (Higher-Order Spectral Nonlinear Wave Toolkit)**  
>   GitHub / documentation: [HOS-NWT Repository](https://gitlab.com/lheea/HOS-NWT)
> 
> - **Grid2Grid (Grid-to-Grid interpolation toolbox)**  
>   GitHub / documentation: [Grid2Grid Repository](https://gitlab.com/lheea/Grid2Grid)

```bash
## Step 0: Load required modules (If you are in Server/install these libraries)
module load OpenFoam/v2206/gcc-11.4.0-binutils-2.40-openmpi-4.1.6
module load fftw/3.310

## Step 1: Clone the repository
git clone https://github.com/Saliyar/FoamWaveTankv2.git
cd FoamWaveTankv2

## Step 2: Edit bin/bashrc to set WAVES_DIR
# This will open nano; update WAVES_DIR to your FoamWaveTankv2 path, save and exit
nano bin/bashrc

## Step 3: Update Grid2Grid FFTW path (Update with your path)
sed -i 's|set(FFTW3_LIB_PATH  .*|set(FFTW3_LIB_PATH  /appl9/fftw/3.3.10/XeonE5-2660v3/gcc-13.3.0/lib)|' Grid2Grid/Grid2GridOptions.cmake

## Step 4: Compile
./Allwmake
   ```

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

