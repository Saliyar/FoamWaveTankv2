# Phase 3: unit checks (partial, stopped at first failed rule)

Build: libwaves2Foam and libswenseWaveTank compile and link against OpenFOAM v2506
(worktree bin/bashrc, OFPLUSBRANCH=1, GSL 2.7 from /services/tools/gsl/2.7, Grid2Grid
lib copied from the main tree). Output to scratch, nothing installed.

Air side, standalone (doc/swense_analytic/test_air.cpp, nondimensional, M = 32, L = 2 pi, c = 1):

| Check | Value | Rule | Result |
| --- | --- | --- | --- |
| Periodic vs cosine basis, even eta and even eta_t | max abs dA = 1.4e-18 | 1e-12 | pass |
| A_t: rotation rule vs central difference, tau = T/1000, orders 1-3 | max abs diff 6.2e-7, relative to max abs A_t = 4.7e-6 | relative 1e-5 | pass |
| Air kinematic mismatch at z = eta, order 1 | 1.18e-2 | decreasing | pass |
| Air kinematic mismatch at z = eta, order 2 | 1.56e-3 | decreasing | pass |
| Air kinematic mismatch at z = eta, order 3 | 6.22e-4 | decreasing | pass |
| Ratio order 2 / order 3 | 2.50 | above 3 | FAIL |

Water-side checks (divergence, kinematic FSBC, dynamic FSBC) were not run: they need the
OpenFOAM test application and a streamFunction case, which is not built yet.
Phase 4 is not started.
