# Phase 0: facts (read from 8848341, no runs)

Base commit: 8848341 (head of branch SWENSE). Work branch: SWENSE_analytic, created
as a git worktree at FoamWaveTankv2_SWENSE_analytic so the MoorDyn_sixDoF working tree is untouched.

## 1. swenseIncident takes any externalWaveForcing; Grid2Grid fast path is a dynamic_cast
Confirmed. `swenseIncident.C:196`: `g2g_ = dynamic_cast<const waveTheories::Grid2Grid*>(&ext_);`
Non-Grid2Grid forcings take the per-point path (`swenseIncident.C:558-562`: `ext_.eta(etaPts_[r], t)`,
and the per-point `ext_.U`, `ext_.pExcess` calls at lines 340, 354).

## 2. Grid2GridCoeffs is read unconditionally
Constructor (`swenseIncident.C:89-106`) reads `refPos_` (referencePosition), `theta_` (referenceTheta, degrees),
`refTime_` (referenceTime) from `waveProperties/Grid2GridCoeffs`. Tsoft is checked at `swenseIncident.C:167`.
Every use:
- `refPos_`: `toWave` (x - refPos_, line 303), `waterP` (line 352), `airUP` (line 370), `buildTables`
  (line 303/367), `pExcess` path (line 485), `pExcess` at 585/595, 676.
- `theta_`: `toWave` (line 317-318), `toMesh` (line 327-328) only.
- `refTime_`: `tAir = t + refTime_` (line 543) only.
- Also the Air2Air block reads `farfieldCoeffs/Air2Air/modesFile` unconditionally (line 178-180); this is Phase 2.

## 3. Units and convention of pExcess
- swenseIncident: `waterP = rhoWater_*(ext_.pExcess(x,t) - magG_*zw)` (line 354). Dimensionally, `pExcess`
  is per unit density, [m2/s2]; the hydrostatic term `g z` is subtracted after it, so the consumed quantity is
  `pExcess = p/rho + g z` (Bernoulli-form excess; hydrostatic included, then removed by the `-g z` term).
- streamFunction::pExcess (streamFunction.C:167-192): `res = R_ - Z|g| - 0.5|v|^2 - (x-ref).g`, then
  `res *= rhoWater_`. With returnZ = z - seaLevel and g pointing down, the Z and (x-ref).g terms
  combine to leave `R_ - 0.5|v|^2 + |g| seaLevel` with R_ = 0.5|v_c|^2 + |g|(E - seaLevel). So
  `theory.pExcess = rho*(p/rho + g z)` in Pa.
  **Adapter conversion: `pExcess_swense(x,t) = theory.pExcess(x,t) / rhoWater`.** Then
  `rho*(pExcess_swense - g z) = p`, as required.
- Grid2Grid: `Grid2Grid.C:437`: `pExcess = pr*f + rhoWater_*|g|*seaLevel`. The second term is Pa
  (it contains rho), so Grid2Grid returns Pa in its own code. The Fortran `getHOSPd`
  (ThirdParty/Grid2Grid/src/libGrid2Grid/modGrid2Grid.f90:388) is documented as
  `p = -rho phi_t - 0.5 rho|U|^2` (doc/interface.tex:239), which is Pa.
- **Open inconsistency (not resolved here, not changed):** swenseIncident subtracts `g z` (m2/s2) from
  Grid2Grid's value, which dimensionally needs per-density input, while Grid2Grid.C returns Pa. The HOS path
  must stay bit-identical, so this is flagged, not fixed. A measurement (HOS `gethospd` value vs rho g a at a
  point) would settle which one is right. Not done in this task.

## 4. Who reads the forcing (targets of farfield relaxation and BCs)
All go through the single `externalWaveForcing` object looked up by name "externalWaveForcing":
- `swenseWaveAlphaFvPatchScalarField.C:86-102`: `ext.eta(x,t)` (blendWeight).
- `swenseWaveVelocityFvPatchVectorField.C:91-128`: `ext.U(...)`.
- `swenseSurface.C:47`: `ext_.eta(pp, t_)`.
- `swenseIncident` (the water side and the eta cache).
Only `Air2AirWrapper`/`Air2AirEnhanced` read the modes file (`readModes`, Air2AirEnhanced.C:53, 388).
Nothing else reads Grid2Grid or the modes file directly. Relaxation zones use swenseIncident.

## 5. Air2AirEnhanced mode source and basis
- Constructor `Air2AirEnhanced(modesFile, dict)` opens the file (Air2AirEnhanced.C:209-212).
- `readModes(it, e, et, Mx, My)` (Air2AirEnhanced.C:53) reads eta and eta_t modes per record.
- `update(t)` (header line 103): the record window is a Grid2Grid-style cache (`cache_`, `slot`, nInterp=4,
  Air2AirEnhanced.C:271-424); coefficients are solved once per HOS record in a two-slot cache, interpolated
  linearly in time, and dA/dt is the record difference over dt_HOS (Air2AirEnhanced.H:19-20).
- Basis in analyticAirCoeffs.H: cosine tank basis on [0, Lx] x [0, Ly]:
  `eta = sum e_mn cos(kx_m x) cos(ky_n y)`, `kx_m = m pi/Lx` (lines 11-14). Phase 2 needs a periodic
  variant (cos and sin, k_m = 2 pi m / L) beside it. Not started.

## 6. setWaveParameters streamFunction options
From `src/waves2FoamProcessing/preProcessing/setWaveProperties/regular/streamFunctionProperties/streamFunctionProperties.C`:
- `height` (H), `depth` (h), `specifyPeriod` (true: `period`; false: `waveLength`).
- `specifyEuler` true: `eulerVelocity`; false: `stokesVelocity`.
- `N` (number of Fourier modes), `Niter` (Newton iterations).
- `direction` (read at line 207).
- Writes `streamFunctionCoeffs` with `depth`, `period`, `uBar`, `nModes`, `omega`.
The theory reads `<name>Coeffs/waveType` through `waveTheory::New(name, mesh)`
(newWaveTheory.C:52-58), so the dictionary `waveTheoryForcingCoeffs` with `waveType streamFunction`
works if the coefficients are copied in from the setWaveParameters output.
streamFunction's Tsoft defaults to the period (streamFunction.C:63); with Tsoft 0, `factor()` returns 1 once
time >= Tstart (streamFunction.C:90-102).

## Phase 0 also needs (not in the spec list, found while reading)
- The worktree's ThirdParty/Grid2Grid is an empty submodule; the built libGrid2Grid.so and the module
  environment are on /zhome (see doc/swense_port/04_build_and_run_commands.md). Build is not yet reproduced
  from this worktree.
