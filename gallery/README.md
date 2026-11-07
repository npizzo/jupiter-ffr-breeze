# Output gallery

Figures below are generated from saved Breeze fields. They are preliminary numerical tests, not a reproduction of Jupiter's folded filamentary regions.

## GPU convection pilot

128 × 128 × 64 cells; four simulated hours on one NVIDIA L40S. Convective columns form below the stable cap. This is a transient dry-convection test; folded-filament morphology is not established.

[Animation and run parameters](pilot/README.md) · [Validation and cost](../docs/results.md)

![GPU pilot final state](pilot/latest.png)

## CPU smoke test

16 × 16 × 24 cells; 60 seconds of simulated time. This checks construction, stepping, and output. The initial random perturbations remain visible at this short duration.

[Animation and run parameters](cpu-smoke/README.md)

![CPU smoke-test final state](cpu-smoke/latest.png)

The scientific target is sustained Jupiter-like FFR morphology.
