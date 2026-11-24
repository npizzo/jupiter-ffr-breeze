# Short compressible energy-budget resolution test

Authorized 24 November 2026. Three branches from the same preserved native checkpoint of clean compressible job 2214011_0, at t=200. Duration 10 time units. Existing runs and source remain untouched.

| Case | Grid | Maximum outer timestep |
|---|---|---|
| Baseline | 256 x 256 x 80 | 1/120 |
| Refined | 384 x 384 x 120 | 1/180 |
| Smaller timestep | 256 x 256 x 80 | 1/180 |

All use the unchanged dry three-dimensional compressible potential-temperature model, 4 x 4 x 1 periodic horizontal box, impermeable top/bottom, full rotation at 67.5 N, Omega=3, beta=0, bottom heat flux 0.03125, local cap restoration and top conductive flux, WENO9 and no closure. No jet, damping, energy correction, new perturbation, or change to reference-state construction. This is a short budget experiment, not a new LSV campaign.

Preparation runs on Unity CPU: restore native donor, preserve Float64 prognostic fields, verify exact same-grid reinstallation, conservatively remap density, potential-temperature density, dye density and periodic horizontal momentum dual volumes using minmod piecewise-linear reconstruction. Vertical face momentum is linearly interpolated with exact zero at impermeable boundaries. All conservation errors, total energy jump, kinetic energy change and initial fine-grid Mach are recorded. Do not force energy equality or subtract the evolving physical residual. Reject transfer if kinetic energy changes by 10% or more, or the total-energy jump exceeds 10% of the 10-unit bottom heat input. Actual production transfer must pass before GPU submission.

Each branch starts its own energy ledger at its actual mapped initial energy. Power from local cap cooling, conductive top loss and bottom heating is accumulated at every step using trapezoidal integration. An independent one-unit power sampling replicates the parent diagnostic's integration cadence. Energies, temperature/density positivity, mass and Mach are recorded every 0.1 units; three-depth vorticity/tracer frames every 0.5; rolling checkpoint every 5. The entire energy ledger serializes with native checkpoints.

Interpretation: compare residual slopes over the whole interval and excluding the first 1-2 units (grid-transfer adjustment), with the absolute error normalized by common bottom input. Fine versus smaller-dt separates spatial effects at matched maximum dt; smaller-dt versus baseline isolates time-step effects. Inspect actual dt histories. Spatial refinement also changes the resolution of the imposed bottom flux and cap profile; a result implicates the spatial discretization broadly, not WENO alone. A ten-unit result is a local sensitivity check, not convergence of long-term vortex statistics.

Mandatory full-driver CPU parsing/advance/output/readback/restart and remap tests precede submission; exact fingerprints include source, scripts, settings, environment, pinned Julia/Breeze, prepared fields, and checkpoint provenance. Each allocation repeats small GPU checks before production. Three concurrent array slots; each allocation at most 90 minutes including startup, production capped at 80 minutes with graceful checkpoint margin; no automatic continuation. No changes to shared dependencies.

Version 2 explicitly refreshes diagnostic density, temperature and velocity after native checkpoint restoration, before extracting the donor state. The first CPU-only preparation stopped at the exact same-grid audit; no GPU jobs were released from version 1.

Version 3 repairs a GPU-only diagnostic failure: Statistics.mean on a nested GPU slice attempted forbidden scalar access. A preallocated GPU top-plane buffer and GPU sum now perform the same mean. Scalar access remains disabled. The previous array stopped in small GPU checks before any full-size branch; total failed allocation time was about 14 minutes. The first allocation must pass the full small GPU/restart/output exercise before releasing the remaining two cases; it then continues directly into the baseline branch. Physics, state transfer and comparison settings are unchanged.
