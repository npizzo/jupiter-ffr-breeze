# Revised convection controls: completed results

All six 48-hour cases and four short validation cases completed successfully on Unity. Correcting the upper heat sink and cap mixing fixes the earlier thermal-budget problems, but does not establish sustained bottom-driven convection or folded filamentary regions. The heated and unforced cases remain remarkably similar: the initial unstable stratification dominates this diagnostic window.

[Figures and data](../gallery/controlled48/README.md) · [Experiment design](controlled48.md)

## Quantitative comparison

“Second day” means the time average over hours 24–48, computed from five-minute diagnostics.

| Case | Second-day vertical RMS speed [m/s] | Final vertical RMS speed [m/s] | Cap contrast / initial | Measured run time [s] |
|---|---:|---:|---:|---:|
| Heated + cooled, 0.2 K | 0.17114 | 0.11011 | 1.0081 | 192.4 |
| Cooled only, 0.2 K | 0.17031 | 0.11017 | 1.0081 | 185.9 |
| Heated + cooled, 0.02 K | 0.19524 | 0.12264 | 1.0092 | 186.5 |
| Cooled only, 0.02 K | 0.19504 | 0.12317 | 1.0091 | 182.3 |
| Unforced, 0.02 K | 0.19504 | 0.12316 | 1.0092 | 183.6 |
| Heated + cooled, 0.02 K, doubled vertical resolution | 0.20216 | 0.12417 | 1.0106 | 365.7 |

The perturbation labels denote amplitude parameters, not actual initial RMS temperature fluctuations. Cap contrast is measured between fixed physical heights of 82.5 and 117.5 km, interpolated from the horizontal-mean profiles, so both grids use the same interval.

Bottom heating raises second-day vertical RMS speed by **0.487%** and **0.104%** relative to its matched cooling-only control at the larger and smaller perturbation amplitudes. The corresponding kinetic-energy changes are **1.05%** and **1.14%**. At the smaller amplitude, removing the upper cooling as well produces almost no further change. Both imposed thermal powers are therefore weak controls on the bulk motion during these 48 hours.

Local upward/downward speeds reach absolute maxima of 16.7–21.3 m/s during the initial burst, then decay. The final horizontal vorticity slices contain small vortices and elongated structures, but they have not been quantitatively matched to observed FFR morphology. These are dry dynamical fields, not cloud images.

## What the revision fixes

The upper input remains exactly −5.4 W/m² to floating-point precision in every cooled case. It no longer supplies the roughly +56 W/m² seen late in the previous relaxed-profile experiment. The unforced case has zero prescribed thermal input.

The revised cap retains its initial temperature contrast to about 1%, versus approximately 49% retention in the previous small-perturbation heated case. Its final mean-profile RMS departure is only 0.020–0.022 K, versus 1.22 K previously. The minimum estimated diffusive heat flux is about −3 W/m² instead of −491 W/m². This latter quantity is internal redistribution, not an additional external heat source. The chosen diffusivity is still an idealized sensitivity parameter, not a calibrated Jovian closure.

The density-weighted potential-temperature inventory agrees with the integrated sources to a relative residual below 1.3 × 10⁻⁸ in the forced cases. The unforced case drifts by only 7.7 × 10⁻⁷ kg K/m² in absolute inventory. Equal and opposite enthalpy powers do not cancel in potential-temperature units because the Exner function varies with height. These checks do not establish a closed discrete total-energy budget.

## Resolution and computational cost

Doubling vertical resolution from 64 to 128 levels changes second-day vertical RMS speed by +3.55% and kinetic energy by −2.87%. The initial noise has a different spatial realization on the changed grid. This is a bulk sensitivity check, not proof of grid convergence; horizontal resolution was unchanged, and thin boundary layers may remain unresolved.

Unity job 66917135 completed in **24 minutes 19 seconds** on one L40S GPU, including validation, initialization, all six runs, and output. The standard cases each took about 3.0–3.2 minutes inside the solver; the refined case took 6.1 minutes. These are measured costs of these small experiments, not projections for a larger or equilibrated Jovian simulation. The interrupted A16 allocation is retained separately as partial output and is excluded from the comparisons.

All cases used the same solver source SHA-256, recorded in `summary.json`. The model stayed dry anelastic on the 55°N f-plane, with beta = 0 and no imposed shear. No alternative dynamical model was tested here.

## Implication for the next experiment

Before adding beta or shear, a useful next diagnostic would start from a nearly neutral lower layer and compare imposed heating/cooling against an unforced control. That would remove the current large initial instability as an energy reservoir. Its forcing, diffusivity, rotation, and runtime should be selected through explicit nondimensional parameters and measured onset/transport behavior. Merely increasing forcing for speed would change the dynamical regime and needs to be documented as a sensitivity experiment. This next experiment has **not** been configured or submitted.
