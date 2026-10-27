# SWENSE analytic incident field: result (stopped at first failed rule)

Reference wave (HOS regular case, modes_HOS_SWENSE_clean.dat, x = 0 gauge):
T = 1.82721 s (zero-crossing mean over 54 crossings), H = 2 a1 = 0.175838 m (a1 = 0.087919 m), h = 5 m,
zero mean current. Script: doc/swense_analytic/water_check/hos_gauge.py.

| Check | Value | Rule | Result |
| --- | --- | --- | --- |
| Air: periodic vs cosine basis (even input) | max abs dA 1.4e-18 | 1e-12 | pass |
| Air: A_t rotation vs central difference | relative 4.7e-6 | 1e-5 | pass |
| Air: kinematic mismatch at z = eta, orders 1 / 2 / 3 | 1.18e-2 / 1.56e-3 / 6.22e-4 | decreasing, order 3 < 1e-2 (revised) | pass |
| Water streamFunction (N 20): div U at 50 points | 1.3e-9 of k a omega | 1e-6 | pass |
| Water streamFunction (N 20): kinematic FSBC, 200 points | 1.12e-6 of a omega | 1e-6 | **FAIL** |
| Water streamFunction (N 20): dynamic FSBC | 5.7e-6 of g a | 1e-5 | pass |
| Water streamFunction (N 30, Niter 100, one increase per spec) | kinematic 1.116e-6 (no change) | 1e-6 | **FAIL, stop** |
| Water stokesFifth: div U | 1.3e-9 | 1e-6 | pass |
| Water stokesFifth: kinematic | 8.3e-5 | 1e-6 | FAIL |
| Water stokesFifth: dynamic | 1.87e-5 | 1e-5 | FAIL |
| Water stokesSecond: div U | 1.3e-9 | 1e-6 | pass |
| Water stokesSecond: kinematic | 1.70e-2 | 1e-6 | FAIL |
| Water stokesSecond: dynamic | 6.8e-3 | 1e-5 | FAIL |
| Phase 4 R0 (HOS regression, before vs after) | not run | identical fields | not run |
| Phase 4 R1 (streamFunction, 200 steps) | not run | SWENSE check | not run |
| Phase 4 R2 (streamFunction, 40 s) | not run | R1 + gauges | not run |

Check settings: t = T/4 for all water checks (the theory eta is zero for t < Tstart = 0, so a
central time difference at t = 0 is invalid); eta_t by central difference with tau = T/1e5; eta_x by
central difference with dx = 1e-4 m; divergence step 1e-4 m; dynamic check p/rho = pExcess - g eta at eta.
Stokes theories were run with the same H, T, h (stokesDrift 0, Tsoft 0).

Commits: 6dbff05 (waveTheoryForcing, seaLevel-independent Grid2GridCoeffs), d6a6116 (periodic air basis, source theory).
Open: the streamFunction kinematic residual is 1.1e-6 and does not change with N; the cause is not measured.
Build: scratch libraries; the worktree bashrc with OFPLUSBRANCH=1, GSL 2.7, gcc 14 libstdc++ at run time.
Open item: Grid2Grid pExcess units per your note (not changed).
