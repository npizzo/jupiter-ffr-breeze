# Neutral initial layer: completed results

All four 48-hour runs completed and passed the output, thermal-inventory, and initial-neutral-profile checks. Removing the initial mean instability exposes a weak response to bottom heating. The heated case grows late while the controls decay, but the short horizontal wavelengths require better resolution before interpreting the patterns physically. Sustained, developed Jovian convection or FFR morphology is not established.

[Snapshots and diagnostics](../gallery/neutral48/README.md) · [Design](neutral48.md) · [Full summary](../gallery/neutral48/summary.json)

## Measured response

Second-day averages cover hours 24–48. RMS vertical velocities are volume averages at cell centers, so they differ from local maxima in a snapshot.

| Case | Second-day w RMS [m/s] | Final w RMS [m/s] | Final maximum resolved upward heat flux [W/m²] | Run time [s] |
|---|---:|---:|---:|---:|
| Heated + cooled | 0.000471472 | 0.00127041 | 0.000342424 | 201.8 |
| Cooled only | 1.50477e-05 | 7.19606e-06 | 2.1658e-09 | 181.2 |
| Unforced | 1.50477e-05 | 7.19603e-06 | 2.16574e-09 | 182.3 |
| Heated + cooled, Nz=128 | 0.000347938 | 0.000943719 | 0.000190614 | 357.1 |

At the standard resolution, final vertical RMS speed is **1.270 mm/s** in the heated case versus **0.0072 mm/s** in the cooling-only control and **0.0072 mm/s** in the unforced control. The heated/control contrast is therefore clear in these simulations, despite the small absolute speeds. The strongest local vertical speed in the standard heated run reaches **1.365 cm/s**. The cooling-only and unforced runs test the decay of the seed disturbance rather than an exactly motionless initial state.

The heated case's final maximum resolved convective heat flux is only **0.0063%** of the imposed 5.4 W/m². This is a maximum over height, not a boundary flux or a complete heat-transport budget. Most of the prescribed input has not become resolved convective transport. The 48-hour experiment shows early growth, not thermal equilibrium or developed turbulence.

## Resolution limits

The horizontal grid spacing is 9.375 km. In the standard heated run, a one-dimensional spectrum of the final vertical-velocity section below 40 km peaks at a horizontal wavelength of **20.69 km**, only **2.21 cells**. About **50.2%** of the section's nonzero-mode power lies at wavelengths shorter than four cells. This is a substantial resolution warning. It does not establish that the structure is a numerical artifact, but its morphology and growth cannot yet be treated as resolved physical predictions.

Doubling vertical resolution changes second-day vertical RMS speed by **-26.20%**. It leaves the horizontal grid unchanged. The refined case also uses a different spatial realization of the seeded noise, so this comparison is a bulk sensitivity check, not a convergence measurement. Its section-spectrum power below four horizontal cells is **35.0%**. A horizontal-resolution check is needed before extrapolating these weak patterns toward FFRs.

## What stayed controlled

The mean lower-layer initial N² was zero; the temperature-noise amplitude parameter was 10⁻⁵ K. Standard-resolution cases shared exactly the same initial seed and reference state. The stable cap, tapered thermal diffusivity, viscosity, Jovian rotation, box, and dry anelastic model were retained. Beta and imposed shear stayed zero.

The upper source remained −5.4 W/m² in cooled cases and zero in the unforced case. All forced potential-temperature inventory residuals were below **1e-08** relative to integrated absolute sources. The unforced absolute drift was **6.4e-07 kg K/m²**. This verifies source accounting for density-weighted potential temperature, not exact total-energy closure.

The stable cap's temperature contrast between 82.5 and 117.5 km remained close to its initial value; per-case values are in the summary. Motion grows mainly below the cap in the heated-case sections. Near-cap vorticity is therefore not a direct picture of the strongest bottom-region motion.

## Snapshots and cost

The gallery contains vertical velocity and temperature anomaly at 6, 12, 24, and 48 hours for every case, plus final comparisons and horizontal vorticity. Velocity is interpolated from vertical faces to cell centers. Temperature anomaly is relative to the horizontal mean at each height. Final comparisons share color scales across cases; each time sequence uses fixed scales within that sequence. The snapshots are direct solver fields, not synthetic clouds or observed-Jupiter comparisons.

All four runs ran sequentially on one L40S. The measured solver costs are in the table; total job timing and host-memory use are in the [resource record](../gallery/neutral48/resources.json). These costs apply to this small, weak-flow experiment and do not predict a larger equilibrated FFR simulation.

The next useful step is a horizontal-resolution check, or a smaller-domain test that resolves the emerging structures at the same physical parameters, before adding beta or cyclonic shear. A separate [seeded-turbulence experiment](turbulence48.md) now tests how convection modifies existing motion. It retains the same model and does not resolve the neutral test’s horizontal-resolution concern.
