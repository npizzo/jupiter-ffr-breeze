# Output gallery

Figures below are generated from saved Breeze fields. They are preliminary numerical tests, not a reproduction of Jupiter's folded filamentary regions.

## Larger domain and maintained cap: validation and launch

[Campaign design](../docs/large40.md). New results will include three vorticity videos, near 20, 60 and 100 km above the model bottom. The older entries below predate this requirement and have not been reconstructed at missing depths.

## Twenty-day runs: interim thermal finding

Three long runs completed. Part of the cap loses stability late in the integration, complicating inverse-cascade interpretation; beta changes final-day bulk energy by about 1%. Horizontal-resolution comparisons are still underway.

[Interim figure and interpretation](regional20-interim/README.md) · [Campaign design](../docs/regional20.md)

![Thermal evolution](regional20-interim/thermal_evolution.png)

## Five-day 30× convection: completed, with videos

Sustained convection, matched controls, and a vertical-resolution comparison.

[Vorticity videos, snapshots and data](convection120/README.md) · [Interpretation](../docs/convection120-results.md)

[![128-level vorticity movie poster](convection120/z128_vorticity.png)](convection120/z128_vorticity.mp4)

## Stronger convection: six completed runs

30× and 100× heating produce vigorous convective columns and a large change relative to cooling-only controls. The finer 100× run gives similar bulk vertical speeds but substantially stronger vorticity.

[Snapshots and data](strong48/README.md) · [Interpretation](../docs/strong48-results.md)

![Stronger-convection snapshots](strong48/f100_h256_heated_snapshots.png)

## Seeded turbulence: three completed thermal-forcing controls

Identical initial Rossby-0.1 velocity fields isolate the incremental effect of bottom heating. Changes in second-day bulk energy and vertical speed are small over 48 hours.

[Snapshots, spectra and data](turbulence48/README.md) · [Interpretation](../docs/turbulence48-results.md)

![Heated-case vorticity evolution](turbulence48/heated_vorticity_snapshots.png)

## Neutral initial layer: four completed runs

Heating produces weak late growth while the controls decay. Short horizontal wavelengths require further resolution checks before interpreting the patterns.

[Snapshots and data](neutral48/README.md) · [Interpretation](../docs/neutral48-results.md)

![Heated-run snapshots](neutral48/neutral48_heated_snapshots.png)

## Revised forcing and cap mixing: six completed runs

Fixed upper cooling and reduced thermal diffusion preserve the stable cap. Matched controls show that initial instability still dominates; heating-driven FFRs are not established.

[Latest gallery and data](controlled48/README.md) · [Interpretation](../docs/controlled48-results.md)

![Revised convection histories](controlled48/evolution.png)

## Controlled 48-hour convection experiment

Four runs isolate bottom heating and initial-perturbation amplitude. All completed; the comparisons identify forcing and diffusion limitations in the initial setup.

[Comparison gallery and data](convection48/README.md) · [Interpretation](../docs/convection48-results.md)

![Paired convection histories](convection48/evolution.png)

## GPU convection pilot

128 × 128 × 64 cells; four simulated hours on one NVIDIA L40S. Convective columns form below the stable cap. This is a transient dry-convection test; folded-filament morphology is not established.

[Animation and run parameters](pilot/README.md) · [Validation and cost](../docs/results.md)

![GPU pilot final state](pilot/latest.png)

## CPU smoke test

16 × 16 × 24 cells; 60 seconds of simulated time. This checks construction, stepping, and output. The initial random perturbations remain visible at this short duration.

[Animation and run parameters](cpu-smoke/README.md)

![CPU smoke-test final state](cpu-smoke/latest.png)

The scientific target is sustained Jupiter-like FFR morphology.
