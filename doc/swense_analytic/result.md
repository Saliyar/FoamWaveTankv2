# SWENSE analytic incident field: result (stopped at R0)

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
| R0: HOS case, before (8848341) vs after, 10 steps, fields at t = 0.059 s | max abs diff 1.5e-11 (p), 2.5e-12 (pc), 7.3e-12 (p_rgh), 6.7e-16 (alpha), 0 (U) | absolute 1e-14 | **FAIL, stop** |
| R1: streamFunction case, 200 steps | not run | SWENSE check | not run |
| R2: streamFunction case, 40 s | not run | R1 + gauges | not run |

R0 detail: relative to max|p| = 4.66e4 Pa, the p difference is 3.1e-16; p_rgh is 7.4e-13 relative
to its max 9.85. The pc field itself peaks at 1.5e-11, so its difference is round-off of the p scale.
Continuity "sum local" also differs in the 22nd digit. Not investigated (no cause measured).
Solver runs: serial, login node, v2506, LD path with gcc 14 libstdc++ and CUDA 11.0.3 (for libmpi).
Build: both builds compiled from source in their own worktrees (before: 8848341; after: this branch).
Per-field comparison: doc/swense_analytic/water_check/cmp_fields.py.

Commits: 6dbff05 (waveTheoryForcing and seaLevel-independent Grid2GridCoeffs), d6a6116 (periodic air
basis, source theory), 065f6f2 (water checks).
Open: Phase 4 R1 and R2 and the body cases are not run. Grid2Grid pExcess units not changed.
