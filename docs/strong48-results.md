# Stronger heating: completed convection–turbulence comparison

All six 48-hour simulations completed. Increasing bottom heating produces a large, sustained response over the second day: vigorous vertical columns, positive buoyancy work and much more kinetic energy than in the matched cooling-only controls. This establishes a useful convection-driven regime in the current model; it does not establish Jupiter-like folded filamentary morphology.

[Snapshots and data](../gallery/strong48/README.md) · [Experiment and parameter conventions](strong48.md)

## Measured response

Second-day averages cover hours 24–48. Vertical speed and vertical-vorticity RMS are volume averages over the domain. Column kinetic energy uses density-weighted face velocities. The final resolved upward heat flux is a maximum over height at 48 hours, not the net boundary flux or an exact discrete heat-budget closure.

| Case | Mean column KE [J/m²] | Mean vertical RMS [m/s] | Mean ζ RMS / f | Final maximum resolved upward heat flux [W/m²] |
|---|---:|---:|---:|---:|
| 30× heated | 179693 | 1.152 | 0.379 | 120.056 |
| 30× cooling only | 19541 | 0.122 | 0.045 | 0.136 |
| 100× heated | 503446 | 1.977 | 0.582 | 431.151 |
| 100× cooling only | 19554 | 0.122 | 0.045 | 0.048 |
| 100× heated, 256² | 539765 | 2.069 | 0.914 | 437.260 |
| 100× cooling only, 256² | 19555 | 0.123 | 0.045 | 0.014 |

At 30× heating (162 W/m²), second-day kinetic energy is 9.20 times its cooling-only control, and vertical RMS speed is 9.44 times the control. At 100× heating (540 W/m²), those ratios are 25.75 and 16.19 on the 128² horizontal grid, and 27.60 and 16.85 on the 256² grid. This is a much larger response than the +0.117% energy change in the 5.4 W/m² baseline.

Second-day diagnosed buoyancy work is about 19.7 W/m² at 30× and 72 W/m² at 100×, versus approximately zero in the controls. Peak absolute vertical speeds over the entire integrations reach 14.6 m/s at 30×, 20.3 m/s at 100× coarse, and 21.3 m/s at 100× fine. These peaks are distinct from the domain RMS values in the table. The heated runs develop an initial overshoot and then a substantially flatter second-day energy history; 48 hours is insufficient to establish statistical stationarity.

## Rotation, resolution and morphology

Jovian f = 2.88096 × 10⁻⁴ s⁻¹ remains fixed, with beta = 0 and no imposed shear. Initial eddy Rossby is 0.1 in every run. The evolved vertical-vorticity RMS divided by f reaches 0.38 at 30×, 0.58 at 100× coarse and 0.91 at 100× fine (second-day means). Thus the convection does not remain uniformly at the initial small Rossby number. The finest run has order-one vorticity relative to f, with larger local excursions; a strongly geostrophic regime is not demonstrated.

Halving horizontal spacing changes the 100× heated second-day vertical RMS speed by **+4.66%**, kinetic energy by **+7.21%**, and vertical-vorticity RMS by **+57.03%**. Final maximum resolved upward heat flux changes by **+1.42%**. Bulk transport is much less sensitive than vorticity, but two horizontal grids and one vertical grid do not establish convergence. Fine gradients and morphology remain unresolved concerns. Discrete sampling and normalization differ slightly across grids; the thermal noise realization also changes. Within each heated/control pair, the initial state and parameters match apart from bottom heat input.

Vertical snapshots show convection filling much of the initially neutral 0–80 km layer, with relatively little vertical motion in the cap. Horizontal maps near 20 km show numerous small vortices and short elongated features. They do not show an established large-scale folded-filament system. The slice spectra document redistribution of horizontal kinetic energy; they are not a spectral-transfer calculation or proof of an inverse cascade.

The mean potential-temperature contrast between 82.5 and 117.5 km remains within 1% of its initial value in all cases. This indicates that the broad stable cap survives, not that its full profile is unchanged: the cap RMS temperature-profile departure is about 0.25 K at 30× and 0.82 K at 100×, largely associated with prescribed cooling.

## Validation and cost

All six cases passed finite-output, completion, prescribed-source, potential-temperature-inventory and initial-turbulence amplitude checks. Initialization checks also cover discrete mass continuity and seed reproducibility. Analysis verifies identical saved initial velocity maps, energy profiles and mean thermal profiles within pairs, plus matching source hashes and physical parameters. Final relative potential-temperature inventory residuals are below 4 × 10⁻¹⁰; this is not a full total-energy closure test.

Unity job 66919897 completed in **42 minutes 7 seconds**, using one L40S and about 4.1 GiB peak host memory reported by Slurm. Coarse solver runs took 3.1–3.4 minutes each; fine runs took 13.4–13.6 minutes each. See [resources](../gallery/strong48/resources.json) and [machine-readable metrics](../gallery/strong48/summary.json). Raw hourly fields remain on Unity; compact figures and full five-minute diagnostic histories are published here.

## Next useful experiment

The stronger forcing achieves the intended large convective modification. Before adding beta or cyclonic shear, the useful next step is a longer, resolution-controlled comparison: refine the 30× case (which retains lower vorticity relative to f), assess vertical resolution as well, and follow whether the small structures organize over several more days. Additional seeds would test robustness. These are recommendations; no additional campaign was submitted with this report. The enhanced heat fluxes remain idealized sensitivities, not estimates of Jupiter's actual forcing.
