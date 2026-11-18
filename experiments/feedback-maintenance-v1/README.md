# Does direct signed reinforcement maintain the large vortices?

Authorized November 18, 2026. All analysis, simulation and rendering take place on Unity. This is an isolated two-case intervention on the completed 512² polynomial experiment, not a new convection or beta-plane campaign.

Both cases restart the identical mature positive-f checkpoint at dimensionless t=400 and run to t=450. The periodic domain is 2π×2π. Ordinary kinematic viscosity remains 0.0001; the original linear instability forcing has multiplier a k² in 12≤|k|≤16, excluding kx=0, with a=0.0002. Polynomial coefficients remain b=0.002, f=1, c=0.0002. There is no beta, additional drag, hyperviscosity, moisture or changed initial velocity. The cubic vorticity source is the previously authorized saturation hypothesis, not velocity drag.

The full source is Q=P0(bf ζ²−c ζ³), with P0 removing spatial mean. The ablation changes only the quadratic output:

    full:              P0(bf ζ²) − P0(c ζ³)
    no_low_quadratic:   (I−P_L) P0(bf ζ²) − P0(c ζ³)
    P_L:               0 < |k| ≤ 2

The cubic source remains unchanged at every scale. Both cases retain identical numerical operators and adaptive timestep criteria. The full control calls the original forcing path; the ablation adds one cached FFT per Runge–Kutta stage. No old source file or checkpoint is modified.

This tests maintenance of an existing state. Failure to maintain it after the intervention would establish reliance on direct large-scale quadratic feedback in this model. It would not prove that direct convective forcing of large modes is unphysical, nor establish whether the state could form from noise. A single pair is not a statistical robustness or convergence test.

## Existing-output analysis

The completed cubic-only and ±f polynomial cases are analyzed over t=0–400. Original source-work scalars are sampled every 0.5, and 81 saved vorticity fields per case are inspected every 5 time units. Native staggered inverse-curl reconstruction independently recovers the recorded energy and quadratic/cubic work within field-storage precision. The advection reconstruction is separately checked against the native Oceananigans operator in the maintenance validation.

At t=400 the positive-f case has E=0.062205, E(0<|k|≤2)=0.044749 and RMS speed 0.35272. The box-crossing time is 17.81, so the 50-unit continuation covers about 2.8 initial crossing times. The major cyclonic core count is 2 under the declared fixed threshold (twice original RMS vorticity) and periodic connectivity; this is descriptive, not an objective LSV criterion.

Over t=350–400, the measured low-mode energy change is +0.011115. Trapezoidal contributions are +0.026600 from quadratic feedback, −0.015239 from cubic feedback, −0.000754 from viscosity and −0.000523 from sparsely sampled reconstructed advection. The unassigned difference is +0.001032 and is retained as sampling/reconstruction error. This is evidence for substantial direct quadratic reinforcement, not an exact time-integrated budget or a claim of no advective transfer at every instant.

## Validation and bounded execution

The release gate checks original-model equivalence; independent source curl and momentum work; removal only of the intended quadratic modes; native advection/diffusion and centered-flux reconstruction; serialized restart; timestep refinement; finite state; output scheduling; and measured execution speed. A separate short production-driver check verifies actual frames, scalars and checkpoint files before production is released.

Production uses one CPU and 4 GB per case, one-hour Slurm limits and 55-minute internal limits. The timestep ceiling remains 0.0025 with the original CFL/reaction guards. Target 450 is a hard campaign bound; no automatic extension or parameter tuning is authorized by this protocol. Checkpoints every 5, vorticity frames every 0.5, instantaneous native low-mode work every 0.25.

Trapezoidal accumulated powers are compared with independent endpoint energy change. Their residual is explicitly retained and is not described as exact every-stage accounting. Recording native terms at sparse times does not remove temporal quadrature error.

The Unity renderer will make one synchronized two-panel vorticity movie over t=400–450, with one fixed shared asinh scale spanning all extrema. These horizontal 2D models have no tracer or vertical layers. Completed media will be validated before gallery publication.

[Earlier polynomial movies](https://npizzo.github.io/jupiter-ffr-breeze/gallery/viscous2d-feedback512-v1/) · [Contrasting stress experiment](https://npizzo.github.io/jupiter-ffr-breeze/gallery/stress2d-v1/)

[Adversarial literature review and proposed next mechanism test](literature-review.md).

## Execution record

Unity validation 1956155 and actual-driver check 1956236 passed. Full control 1956335 and quadratic low-mode removal 1956337 were submitted, with rendering 1956341 dependent on both completing successfully. Each production case has one CPU and a one-hour wall limit. These identifiers record submission, not a claim that the target has completed. See [execution.json](execution.json) and [measured validation](validation.toml).
