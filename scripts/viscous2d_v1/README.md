# Viscosity-only 2D cascade and LSV experiment

User-authorized November 16, 2026. This supersedes the weaker-forcing campaign that retained cubic drag. Importing van Kan-style instability forcing did not authorize importing its nonlinear dissipation. No linear or nonlinear drag is used here.

## Assumptions disclosed before launch

- Horizontal incompressible 2D Oceananigans, periodic 2π square, constant f=1, beta=0. Uniform f does not select vorticity sign in this system.
- Grid 256². Ordinary Laplacian viscosity ν=10⁻⁴ only. The compatible native-grid Laplacian is used; there is no hyperviscosity, drag, cooling or hidden energy sink.
- New isotropic random initial vorticity in 12≤|k|≤16, seed 62041, normalized to RMS velocity 0.2. All four runs use the identical source file. Initial turbulence is not copied from a drag-maintained equilibrium.
- Four cases: unforced decay; instability forcing a=2×10⁻⁴; random-sign stochastic stirring; cyclonic-core stochastic stirring. Forcing acts in 12≤|k|≤16, with k_x=0 excluded as in the existing forcing implementation. Unlike the forcing, the initial annulus includes k_x=0.
- Stochastic pulses every 0.2 model time, same locations in both polarity cases. Expected energy input matches instability work at the common initial state. The measured value is 0.0014841159507395108; actual cross work is recorded. Polarity refers to the core: compensating opposite-sign surroundings preserve zero circulation.
- Dimensionless time 0–200, dt=0.0025. No mapping to Earth/Jupiter days. One realization per case.
- Production is gated on the actual unforced model showing negative native advective flux and increasing absolute large-scale energy, with acceptable spectral tails and closed native budgets. No claim of a broad stationary inertial range or resolution convergence.

The initial integral-wavelength turnover time is 2.2950; E/input is 13.4760. Each stochastic pulse is 12.18% of the initial RMS velocity. Net linear instability growth, including viscous damping on this grid, ranges from 0.01450 to 0.02593 over the forcing band. These values are measured on Unity from the actual initializer.

## Exact gate

The validation job generates the common initial condition and runs the production driver for t=8 unforced and t=2 in each forcing case. The unforced check requires:

1. Mean Π(K)<0 at K=4,6,8 over t≥4, with ≥75% of samples negative at each cutoff.
2. Increase in E(|k|≤4) of at least 1% of the initial total energy. This excludes an interpretation based only on a growing energy fraction during decay.
3. Less than 1% of enstrophy at |k|≥N/3 throughout saved probe states, as a basic spectral-tail screen, not a convergence test.
4. Zero injection and zero drag work, integrated kinetic-energy closure to relative 10⁻⁸, native stepwise vorticity and third-moment checks, and CFL<0.5.

Forcing probes also check matched initial fields, matched random positions and quadratic input, exact pulse work and the requested sign laws. The phase-cycling arithmetic is tested before these probes; production diagnostics sample every 81 steps while pulses occur every 80. Complete groups of 80 samples visit every discrete pulse phase once. No threshold is weakened after seeing a failure.

Production requires PASS.json plus an afterok dependency. Four one-CPU runs, at most two concurrent, each with a six-hour scheduler ceiling and a 5.5-hour internal stop. Partial output is retained. No automatic continuation, grid escalation or dissipation change. Existing 3D/Breeze and beta-plane experiments are unaffected.

## What answers the question

Report absolute energy at box scales, native spectral energy flux, spectra, full and low-pass vorticity skewness, integrated input/viscous work and synchronized vorticity evolution. An unforced flow can itself make large vortices through transient inverse transfer. The comparison therefore asks whether weak forcing enhances or maintains a box-scale vortex and whether imposed small-scale polarity changes its sign statistics. A box-scale vortex in one frame, or a large low-k fraction after most energy has decayed, is insufficient.

All four cases are rendered on Unity as a 2×2 movie with a shared fixed asinh vorticity scale spanning full extrema. These are actual horizontal 2D fields, with no vertical layers or tracer. Output every 0.5 model time alternates between post-pulse and mid-cycle states in the stochastic cases; the inherited run configuration's `movie_sampling` text incorrectly says post-pulse only. The 81-step diagnostic series, rather than the movie frames, supplies phase-balanced statistics. The renderer documents the actual movie cadence. No raw fields are downloaded.

## Connection to rotating convection

A depth-averaged equation in reduced rapidly rotating convection has forcing −〈J(ψ′,ζ′)〉 from depth-dependent motions. Temperature and vertical velocity affect this correlation; they are not represented by directly substituting a random temperature into the 2D vorticity equation. The stochastic and instability terms here are surrogate forcings of the horizontal flow. They omit convection/large-scale-flow feedback and cannot derive the physical origin of an f-dependent sign asymmetry.

- [Reduced convection, equation 2.9](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/on-the-inverse-cascade-and-flow-speed-scaling-behaviour-in-rapidly-rotating-rayleighbenard-convection/65EA7F62ED610264B8588035016CFBF3)
- [Guervilly, Hughes & Jones: large-scale vortices](https://arxiv.org/abs/1403.7442)
- [van Kan et al.: forcing and nonlinear-damping effects](https://arxiv.org/abs/2206.12355)
