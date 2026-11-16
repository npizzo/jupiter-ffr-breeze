# Long 512² forced/unforced runs, without analysis products

Authorized November 16, 2026. The user requested the higher-resolution forced/unforced long runs and explicitly said no comparison movies or diagnostics are needed. This campaign has two cases only; it does not restart the stochastic-stirring cases.

Assumptions stated before launch:

- Horizontal 2D Oceananigans in the same periodic 2π square, f=1, beta=0.
- 512² grid, dt=0.0025, ordinary native Laplacian viscosity ν=10⁻⁴, exactly zero drag.
- Identical refined initial vorticity realization in both cases. The actual initial velocity RMS is 0.1994510449469 after preserving the coarse vorticity Fourier coefficients; no renormalization is introduced.
- Forced case: a=2×10⁻⁴ times k²ζ in the 12–16 band, excluding k_x=0. Control: a=0. No stochastic pulses.
- Target model time t=200, divided into five checkpointed segments ending at 40,80,120,160,200. All segments retain the same parameters. This is dimensionless model time, not planetary days.
- Output only checkpoints, raw vorticity every 0.5 time unit, configuration, basic progress and completion status. No spectra, vorticity-budget reports, comparison movies or diagnostic-analysis jobs are created by this campaign.

## Driver and restart verification

`LinearTerms` applies exactly the existing instability and native Laplacian-viscosity tendencies with drag zero. It combines their linear multipliers and omits the expensive native-budget observation callbacks. A short Unity check compares both cases with the existing audited driver to relative velocity error <10⁻¹⁰, verifies a split/restarted integration against an uninterrupted integration at the same tolerance, and measures runtime after compilation. This is an implementation/restart check, not a scientific parameter survey. Both integrations retain finite-state and CFL<0.5 guards.

The existing short 512² resolution check was still running when this longer campaign was authorized. These runs do not imply a resolution-convergence claim and are not gated on that check's movie or analysis. Only the new driver-equivalence/restart check gates production. Later segments require the preceding segment to complete successfully; a partial segment preserves its checkpoint and stops its dependents. No automatic timestep, forcing, dissipation or grid change is allowed.

One CPU and 16 GB per job, at most one active segment per case. Each job has a six-hour scheduler ceiling and 5.5-hour internal limit. The two case chains can run concurrently. All computation and raw output stay on Unity. A segment's initial saved frame duplicates the preceding segment's endpoint and is labeled with the same global iteration/time.
