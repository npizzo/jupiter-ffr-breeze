# 35°N circulation audit over Earth days 0–300

[Watch the three-depth vorticity/tracer movie](https://npizzo.github.io/jupiter-ffr-breeze/#beta35-300) and [read the analysis beneath it](https://npizzo.github.io/jupiter-ffr-breeze/#ferrel35-300-evidence).

Persistent Eulerian overturning continues through day 300, and the centered momentum diagnostic approaches balance, but eddy-driven Ferrel cells remain unestablished. All analysis and plotting ran on Unity; simulations were not modified. Analysis job 1709118 completed in 3 min 13 s. Only finished publication assets were transferred.

## Method and checks

The same diagnostic as `ferrel35_day200_audit.py` now reads all 2,401 three-hourly outputs, days 0–300, across 30 completed ten-day segments. It checks completeness, unique consecutive times, matching physical settings, finite loaded fields and unchanged positive reference density. Configuration hashes can differ because continuation stopping times differ. Segment metadata and source hashes are recorded in the summary.

Trapezoidal means cover six 50-day windows (0–50 through 250–300), their twelve 25-day halves, and the entire 0–300 interval: 19 windows in total. All 153 scalar metrics from the nine earlier day-200 windows reproduce within rtol=1e-12, atol=1e-14. The whole-run mean includes initial adjustment and is not treated as stationary.

Coordinates are right-handed: x east, y north, z up. The native mass-flux streamfunction satisfies rho*v = -dPsi/dz and rho*w = dPsi/dy, in kg/(m s). The centered momentum estimate is unchanged:

    residual = d(ubar)/dt - f vbar - E_h - E_v - A
    E_h = -d(<u'v'>)/dy
    E_v = -(1/rho) d(rho <u'w'>)/dz
    A = time_mean(-vbar d(ubar)/dy - wbar d(ubar)/dz)

RMS norms and normalized inner products of centered fields are density weighted over |y| <= 2,400 km and heights 10–75 km. Streamfunction similarity is an unweighted normalized inner product over the corresponding native face/corner region. These are not centered Pearson correlations or significance tests.

Visualization uses the same separable Gaussian filter as the earlier audit (standard deviations 200 km in y and 5 km in z, reflected edge values). Metrics are unsmoothed except explicitly labeled quantities. Evolution plots use one scale per field across all six windows; the final-halves figure and whole-run mean have separate shared scales. The latter includes spin-up. The whole-run plotting helper only repairs title layout, recomputes the same trapezoidal mean on Unity and checks it against the audit summary.

## Window comparison

| Earth days | Residual / f-v RMS | Eddy alignment with -f-v | Smoothed alignment |
|---|---:|---:|---:|
| 0–50 | 4.340 | 0.232 | 0.719 |
| 50–100 | 2.347 | 0.042 | 0.710 |
| 100–150 | 1.123 | 0.295 | 0.648 |
| 150–200 | 1.235 | 0.251 | 0.833 |
| 200–250 | 0.830 | 0.576 | 0.970 |
| 250–300 | 0.620 | 0.704 | 0.958 |
| 0–300 | 1.124 | 0.237 | 0.842 |

For days 250–300, RMS accelerations in m/s² are: tendency 6.835e-8, Coriolis 2.797e-6, horizontal eddy 2.274e-7, vertical eddy 2.921e-6, mean advection 1.316e-6 and residual 1.734e-6. Mean advection is appreciable, and the residual is still 62% of the Coriolis RMS. Higher alignment after smoothing does not establish a closed two-term balance.

The final 25-day streamfunction patterns have similarity 0.91297 (0.98020 smoothed). At days 150–200 those values were 0.92445 and 0.87201: repeatability is not improving monotonically at every scale. The interior RMS mean vertical velocity changes from 0.52597 m/s in days 150–200 to 0.28975 m/s in days 250–300. Pattern persistence therefore does not imply stationary amplitude. The zonal velocity remains nearly barotropic, and the strongest overturning remains concentrated near the bottom and sidewalls.

Maximum saved mass-flux reconstruction error is 2.665e-14 and maximum native mass divergence 1.796e-17 in the diagnostic's saved units. These checks concern continuity, not momentum closure.

## Interpretation and limits

This is evidence for persistent Eulerian overturning accompanying zonal organization. The approximate centered budget looks more consistent with eddy compensation late in the run, but it cannot identify the driving mechanism reliably on its own.

The separate day 200–210 diagnostic restart passes the exact native total momentum-budget check and exposes a substantial native-minus-centered advection remainder. This is an operator/interpolation/numerical-transport issue, not automatically a missing physical eddy stress. Its result applies to that ten-day branch; no exact native full-run budget is fabricated for days 0–300. The legacy residual cannot be interpreted as a measured physical force or rescued simply by extending the averaging window.

Thermodynamic indirectness, residual/Lagrangian circulation, resolution convergence and independence from the sidewalls, cap relaxation, sponge and thermal feedback remain untested here. A movie, zonal kinetic-energy fraction or persistent streamfunction does not by itself establish Jovian Ferrel cells or an inverse cascade.

Analysis source: `scripts/ferrel35_day300_audit.py`, SHA256 `0429456b1a352fc013d9f553ef7e04e71d5b037b45246bee3621302b716ae831`. Plot-layout helper: `scripts/ferrel35_day300_whole_plot.py`. Published assets: `gallery/cascade35-full300/ferrel35_day300_{late.png,evolution.png,whole_run.png,summary.json}`.
