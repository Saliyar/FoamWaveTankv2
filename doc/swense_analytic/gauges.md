# Gauge comparison: CFD vs streamFunction theory

Window t = 20-40 s (3390 samples, 10.95 periods). First harmonic and mean level fitted by least squares at omega = 3.43868 rad/s; theory sampled at the same gauge x and times (gaugeTheory on the same coefficients).

Rules: amplitude within 1%, phase within 2 deg, mean within 0.5 mm.

## VoF (waveFoamSWENSE, R2full)

| x (m) | amp CFD (m) | amp theory (m) | amp diff (%) | phase diff (deg) | mean diff (mm) | rules |
|---|---|---|---|---|---|---|
| -8 | 0.08754 | 0.08757 | -0.037 | -0.002 | +0.004 | pass |
| -1 | 0.08748 | 0.08754 | -0.061 | -0.001 | +0.005 | pass |
| 0 | 0.08747 | 0.08752 | -0.064 | -0.000 | +0.005 | pass |
| 1 | 0.08749 | 0.08754 | -0.062 | -0.001 | +0.005 | pass |
| 8 | 0.08752 | 0.08756 | -0.037 | +0.000 | +0.004 | pass |

## isoAdvector (waveIsoFoamSWENSE, R2iso_full)

| x (m) | amp CFD (m) | amp theory (m) | amp diff (%) | phase diff (deg) | mean diff (mm) | rules |
|---|---|---|---|---|---|---|
| -8 | 0.08754 | 0.08757 | -0.033 | -0.001 | +0.002 | pass |
| -1 | 0.08751 | 0.08754 | -0.032 | -0.000 | +0.002 | pass |
| 0 | 0.08749 | 0.08752 | -0.036 | +0.000 | +0.001 | pass |
| 1 | 0.08751 | 0.08754 | -0.033 | -0.001 | +0.001 | pass |
| 8 | 0.08753 | 0.08756 | -0.034 | +0.001 | +0.002 | pass |

