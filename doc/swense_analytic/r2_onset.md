# R2 onset (cancelled run R2b, no new run)

Source: R2b log (t up to 9.70 s, cancelled by me). Cases: tutorial waveFoamSWENSE_coarse, identical
controlDict and fvSolution; waveProperties differs only in externalForcing, the streamFunction block
(N 30) and Air2Air (source theory, nModes 64, modesFile and Ha removed).

1. Warning source: src/swenseWaveTank/swenseIncident.C:499-503 (weights(), "W strictly between 0 and 1").
   Quantity: W = clamp((alphaI - 1e-6)/(1 - 2e-6), 0, 1). It warns when 0 < W < 1 and |z - eta| > guard_.
   guard_ = max(delta, 3 dz_interface) = max(0.024, 3 x 0.02417) = 0.0725 m (log: "guard d_g = 0.0725 m").
   The same test sets evalWater/evalAir at lines 484-485; FatalError at 487-492 if neither applies.
   The 3.9 M warnings are cells in the air side whose alphaI is between 0 and 1 more than d_g from eta.
   The first warning (t = 0.0118 s) is at x = -6.87 m, z = 0.074 m. The first step with
   max|alpha1 - alphaI| > 1e-12 (t = 2.0296 s, value 9.1e-8) has 7 199 warnings within 0.01 s, all at
   z = -0.1812 m, x from -2.43 m upward (listed in the run log).

2. First step with max|alpha1 - alphaI| > 1e-12: t = 2.0296 s, 9.06e-8. The log has no cell location
   for the max; the location is not logged. The warnings above are the only coordinates.
   Max|Uc| first exceeds 1e-10 at the same step (1.9e-8).

3. Per 0.5 s (maxima over the SWENSE check lines in each window; locations not logged):

   window (s) | n | max|alpha1-alphaI| | max|Uc|
   1.5-2.0 | 84 | 1.532e-14 | 2.953e-14
   2.0-2.5 | 85 | 1.558e-07 | 3.336e-08
   2.5-3.0 | 85 | 7.211e-14 | 5.896e-14
   3.0-3.5 | 85 | 4.657e-14 | 7.180e-14
   3.5-4.0 | 84 | 4.452e-14 | 8.475e-14
   4.0-4.5 | 85 | 3.608e-14 | 9.834e-14
   4.5-5.0 | 85 | 2.343e-14 | 1.185e-13
   5.0-5.5 | 85 | 2.554e-14 | 1.313e-13
   5.5-6.0 | 84 | 1.816e-07 | 3.905e-08
   6.0-6.5 | 85 | 3.508e-14 | 1.595e-13
   6.5-7.0 | 85 | 3.730e-14 | 1.744e-13
   7.0-7.5 | 85 | 5.196e-14 | 1.814e-13
   7.5-8.0 | 84 | 5.368e-14 | 2.081e-13
   8.0-8.5 | 85 | 4.319e-14 | 2.164e-13
   8.5-9.0 | 85 | 6.850e-14 | 2.420e-13
   9.0-9.5 | 85 | 9.193e-14 | 2.718e-13
   9.5-10.0 | 82 | 7.117e-14 | 2.916e-13

   The failures are two bursts (2.0-2.5 s and 5.5-6.0 s, about 1.6e-7 and 1.8e-7). Between them
   the check is at round-off (about 1e-14 alpha, 1e-13 Uc). The earlier statement of a steady growth
   to 1e-7 was wrong: this is a burst pattern.

4. Settings. Relaxation zones (waveProperties, same as tutorial): farfieldPatchNames (inlet outlet),
   farfieldDistances (8.76 8.76), blendingDistances (8.76 8.76), spatialExponent 3.5. Domain x from
   -10.11 to 10.11 m (blockMeshDict). Computed zone edges at |x| = 10.11 - 8.76 = 1.35 m (inlet zone
   x < -1.35, outlet zone x > 1.35); not checked against the relaxation code. The first burst site
   (x about -2.4 m, z = -0.181 m) lies inside the inlet zone.
   Column level correction: no levelCorrection* keys in either case. Defaults from
   swenseColumnLevel.H: levelCorrectionColumnLog false, levelCorrectionZMin -GREAT, no skip box.
   Identical in the tutorial, so not a difference between the two cases.
