# SWENSE analytic incident field: result (R0 accepted relative; R1 pass; R2 pass under revised rule)

Reference wave (HOS regular case, modes_HOS_SWENSE_clean.dat, x = 0 gauge):
T = 1.82721 s, H = 0.175838 m, h = 5 m, zero mean current.
Script: doc/swense_analytic/water_check/hos_gauge.py.

| Check | Value | Rule | Result |
| --- | --- | --- | --- |
| Air: periodic vs cosine basis | max abs dA 1.4e-18 | 1e-12 | pass |
| Air: A_t rotation vs central difference | relative 4.7e-6 | 1e-5 | pass |
| Air: kinematic mismatch, orders 1 / 2 / 3 | 1.18e-2 / 1.56e-3 / 6.22e-4 | decreasing, order 3 < 1e-2 (revised) | pass |
| Water streamFunction (N 20): div U, 50 points | 1.3e-9 of k a omega | 1e-6 | pass |
| Water streamFunction (N 20): kinematic FSBC | 1.12e-6 of a omega | 1e-5 (revised) | pass |
| Water streamFunction (N 20): dynamic FSBC | 5.7e-6 of g a | 1e-5 | pass |
| Water streamFunction (N 30, Niter 100): kinematic | 1.116e-6 (unchanged) | info | info |
| Water stokesFifth: div / kinematic / dynamic | 1.3e-9 / 8.3e-5 / 1.9e-5 | info only | info |
| Water stokesSecond: div / kinematic / dynamic | 1.3e-9 / 1.7e-2 / 6.8e-3 | info only | info |
| R0: HOS case, before (8848341) vs after, 10 steps, fields at t = 0.059 s | max abs diff 1.5e-11 (p), 2.5e-12 (pc), 7.3e-12 (p_rgh), 6.7e-16 (alpha), 0 (U); relative to field max: p 3.1e-16, p_rgh 7.4e-13 | absolute 1e-14 as written; accepted as relative round-off (user decision) | pass (relative) |
| R1: streamFunction case, 200 steps (N 30 coefficients, source theory, cutoff auto, cutoffC 1, order 3) | max abs(alpha1 - alphaI) 1.75e-14 (rule 1e-12); max abs Uc 6.9e-14 (rule 1e-10) | SWENSE check | pass |
| R2 (earlier attempts 14717, 14720 stopped; superseded) | max abs(alpha1 - alphaI) first > 1e-12 at t = 2.03 s (9.1e-8 there, 1.8e-7 max by t = 9.7 s); max abs Uc first > 1e-10 at t = 2.03 s (1.9e-8 there, 3.9e-8 max) | R1 rule for the whole run | **FAIL, stop** |

R0 detail: relative to max|p| = 4.66e4 Pa, the p difference is 3.1e-16; p_rgh is 7.4e-13 relative
to its max 9.85. The pc field itself peaks at 1.5e-11, so its difference is round-off of the p scale.
Continuity "sum local" also differs in the 22nd digit. Not investigated (no cause measured).
Solver runs: serial, login node, v2506, LD path with gcc 14 libstdc++ and CUDA 11.0.3 (for libmpi).
Build: both builds compiled from source in their own worktrees (before: 8848341; after: this branch).
Per-field comparison: doc/swense_analytic/water_check/cmp_fields.py.

Commits: 6dbff05 (waveTheoryForcing and seaLevel-independent Grid2GridCoeffs), d6a6116 (periodic air
basis, source theory), 065f6f2 (water checks).
R2b log: no fatal error; 3.9 M "alphaI has drifted" warnings by t = 9.7 s. R1 log: 374 k "alphaI has drifted" warnings (swenseIncident.C:501, W strictly between 0 and 1); not a rule, same in R2.
Open: R2 result; waveIsoFoamSWENSE and body cases wait for go-ahead. Grid2Grid pExcess units not changed.

## R2 full run (job 14736, fresh case R2full_streamFunction), revised rules

| Rule | Value | Rule | Result |
| --- | --- | --- | --- |
| SWENSE: max abs(alpha1 - alphaI), 6780 steps | 2.72e-7 (at t = 35.0-35.5 s) | 1e-6 | pass |
| SWENSE: max abs Uc, 6780 steps | 5.65e-8 m/s (same window) | 1e-6 m/s | pass |
| Growth from burst to burst | bursts (window max alpha) 1.6e-7 at 2.0 s, 1.8e-7 at 5.5 s, 2.7e-7 at 35 s; no monotone trend, many windows at 7e-8 to 2e-7 | no growth | pass (judgement) |
| Gauge amplitude, first harmonic, x = -8, -1, 0, 1, 8 (t 20-40 s) | rel. diff -0.037, -0.061, -0.064, -0.062, -0.037 % | 1% | pass |
| Gauge phase, same gauges | |diff| below 0.001 deg | 2 deg | pass |
| Gauge mean level, same gauges | diff 0.004 to 0.005 mm | 0.5 mm | pass |
| Wall time, serial, 40 s run, this machine (2317 s) | 38.6 min | HOS coarse: README says about 20 min on 8 cores (serial not measured) | not comparable as stated |

Notes: 0 fatal errors; the 'alphaI has drifted' warning is still emitted (about 2 400 per step;
log 2.2 GB). Solver: psiProj averages 231 iterations per step. Cells: 16 262.

## waveIsoFoamSWENSE (isoAdvector), streamFunction case, full 40 s (job 14860, case R2iso_full)

Case: tutorial waveIsoFoamSWENSE_coarse with the R2 waveProperties; fvSolution: cAlpha 1 and icAlpha 0
added to the alpha.water.* and alphaI.* dictionaries (the tutorial lacks them; the base class needs them).

| Rule | Value | Rule | Result |
| --- | --- | --- | --- |
| SWENSE: max abs(alpha1 - alphaI), 6780 steps | 0 | 1e-6 | pass |
| SWENSE: max abs Uc | 0 m/s | 1e-6 m/s | pass |
| Gauge amplitude, x = -8, -1, 0, 1, 8 (t 20-40 s) | -0.033, -0.032, -0.036, -0.033, -0.034 % | 1% | pass |
| Gauge phase | |diff| below 0.001 deg | 2 deg | pass |
| Gauge mean level | diff 0.001 to 0.002 mm | 0.5 mm | pass |
| Wall time, serial, 40 s (1798 s) | 30.0 min (VoF R2 full: 38.6 min) | HOS coarse: about 20 min on 8 cores (serial not measured) | not comparable as stated |

Fatal errors: 0. Drift warnings: 5 865 over the run (VoF: about 3.9 M by 9.7 s, so far fewer).

## stokesFifth, waveFoamSWENSE, 200 steps (case R1_stokesFifth)
phaseSpeed() = omega/K and waveNumber() = K added to stokesFifth and stokesSecond (headers only).
max|alpha1 - alphaI| = 1.5e-14, max|Uc| = 5.2e-14 (R1 rules 1e-12 / 1e-10): pass. Fatal errors: 0.
Air2Air: source theory, L = 5.261 m, retained 10 of 64 modes. stokesSecond: not run (per instruction).
