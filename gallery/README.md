# Output gallery

Figures below are generated from saved Breeze fields. They are preliminary numerical tests, not a reproduction of Jupiter's folded filamentary regions.

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
