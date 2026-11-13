# 35°N circulation audit through day 200

[Watch the three-depth vorticity/tracer movie and read the diagnostic below it](https://npizzo.github.io/jupiter-ffr-breeze/#beta35-200).

The late overturning is more persistent, but the resolved momentum budget does not establish eddy-driven Ferrel cells. This read-only diagnostic ran on Unity (job 1668230, completed in 76 seconds); no simulation was modified. Only finished figures and the numerical summary were transferred for publication.

## Method

The same calculation as the earlier 100-day audit is applied to 1,201 three-hourly saved outputs from Earth days 50–200. Segments 1–10 belong to the original 100-day campaign; segments 11–20 are its checkpoint continuation. All 20 metadata records are complete, physical/diagnostic settings agree, saved times are consecutive, and all loaded fields are finite. The continuation configuration has a different hash because its stopping time changed. Segment metadata hashes are recorded in the summary.

Trapezoidal time averages use nine windows: each of days 50–100, 100–150 and 150–200, and their two 25-day halves. Zonal means and eddy covariances are read from saved diagnostics. The streamfunction satisfies rho v = -dPsi/dz and rho w = dPsi/dy, using the saved native mass-flux reconstruction. Coordinates: x east, y north, z up.

The centered resolved zonal-momentum diagnostic is

    residual = d(ubar)/dt - f vbar - E_h - E_v - A
    E_h = -d(<u'v'>)/dy
    E_v = -(1/rho) d(rho <u'w'>)/dz
    A = < -vbar d(ubar)/dy - wbar d(ubar)/dz >_time

This is not an exact discrete model budget. Viscosity and numerical stresses are omitted; the residual also includes diagnostic errors. RMS norms and centered-field pattern inner products are density-weighted over |y| <= 2,400 km and 10 <= z <= 75 km. Streamfunction pattern similarity is an unweighted normalized inner product on its native face/corner grid within that region. It is not a centered Pearson correlation or a significance test.

Figures use a separable Gaussian filter with reflected edge values, standard deviations 200 km in y and 5 km in z. All columns share the same physical color limits within each figure row. Metrics are unsmoothed except the explicitly labeled smoothed alignment and streamfunction comparison. Filtering can change alignment substantially and does not prove momentum-budget closure.

## Results

| Window (Earth days) | Residual / f-v RMS | Eddy alignment with -f-v | Smoothed alignment |
|---|---:|---:|---:|
| 50–100 | 2.347 | 0.0419 | 0.7095 |
| 100–150 | 1.123 | 0.2954 | 0.6483 |
| 150–200 | 1.235 | 0.2514 | 0.8326 |

The final two 25-day streamfunction patterns have similarity 0.924 (0.872 after visualization smoothing), compared with 0.529 (0.515 smoothed) in the original days 50–100 audit. The earlier residual ratios, unsmoothed eddy alignment and streamfunction maxima are reproduced within 1e-12 relative tolerance.

For days 150–200, RMS accelerations in m/s² are: Coriolis 2.642e-6, vertical eddy 3.530e-6, horizontal eddy 2.664e-7, mean advection 1.766e-6, tendency 1.143e-7 and residual 3.262e-6. The mean advection is not negligible. A better pattern match after smoothing must not be interpreted as a closed two-term balance.

The broad overturning remains concentrated near the bottom and channel walls, while zonal bands are nearly depth-independent. The observed persistence supports persistent Eulerian overturning; it does not establish that eddies drive thermally indirect cells analogous to Jupiter's. The mass-flux reconstruction error remains below 1.85e-14 and maximum native mass divergence below 9.0e-18 in the saved diagnostic units. These continuity checks do not validate the momentum balance.

## Limits and reproducibility

Still needed for a stronger physical conclusion: the complete native discrete momentum budget (including viscous/numerical contributions), a thermodynamic indirectness test, and sensitivity to sidewalls and the stable-cap/sponge/thermal-feedback treatment. Residual or Lagrangian circulation is not measured here. The saved jet-energy fraction includes any uniform zonal component; it is not evidence for inverse energy transfer.

Analysis source: `scripts/ferrel35_day200_audit.py` (executed on Unity as `scripts/ferrel35_day200_audit_v2.py`). Its SHA256 is `326996530321f7401f12292d6054a73933fe00cf56feb5afc5b76c35a1ef0ec8`. Figures and all window metrics: `gallery/beta35-full200/ferrel35_day200_{overturning.png,evolution.png,summary.json}`. The first audit attempt stopped before reading fields because it required identical whole-configuration hashes across the continuation; the successful version checks matching physical settings instead.
