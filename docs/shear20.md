# Cyclonic shear, passive tracer, and thermal control

[Watch the gallery](https://npizzo.github.io/jupiter-ffr-breeze/) · [Submitted jobs](shear20-launch.json)

This new campaign uses a separate driver, `scripts/run_shear.jl`. It does not change the running large40 f/beta pair. Four 20-Earth-day experiments have been submitted in ten-day checkpointed segments, gated on thermal/restart validation and full-grid GPU probes. Submission is not completion; see the launch record for job IDs.

| Case | Grid | Bottom heating | Initial cyclonic shear |
|---|---|---|---|
| shear20_heated | 1024² × 64 | 162 W/m² | peak relative vorticity 0.05 f |
| shear20_control | 1024² × 64 | zero | identical to heated |
| shear20_fine | 1536² × 64 | 162 W/m² | same physical shear |
| shear20_baseline | 1024² × 64 | 162 W/m² | zero |

All cases use a 6000 × 6000 × 120 km beta-plane channel centered at 55°N: f = 2.88095774e-4 s⁻¹ and beta = 2.8216700143e-12 m⁻¹ s⁻¹. Initial turbulence has eddy Rossby number 0.1 and wavelength 300 km. The lower 80 km is initially neutral; the stable cap, free-slip walls, diffusivities and random seeds follow the large40 configuration. Horizontal spacing is 5.859 km or 3.906 km. Vertical spacing remains 1.875 km, so this is a horizontal-resolution comparison, not a complete convergence study.

The barotropic shear is U(y) = -Us sin(πy/Ly), with Us = 0.05 f Ly/π ≈ 27.5 m/s. It satisfies the free-slip sidewall condition and has cyclonic vorticity throughout the channel interior. An initial kinematic pressure anomaly balances its Coriolis acceleration, including beta. There is no sustained shear forcing; the flow can redistribute its initial shear energy.

## Heat accumulation

Refining a grid cannot remove a net external heat input. The new campaign therefore changes the thermal control independently of resolution. A local relaxation maintains the cap reference profile. Every 600 seconds, a controller measures that relaxation's column enthalpy input, Pcap, and the column quantity E = integral(rho cp Exner theta dz). It adds a distributed potential-temperature tendency -alpha theta, with alpha = (Fbottom + Pcap)/E. The external enthalpy source is zero at each controller update; measured residuals between updates are integrated and checked.

This is an idealized feedback sink, not Jovian radiative transfer. The controller can heat if the other sources produce net cooling. Its uniform fractional theta tendency preserves the logarithmic vertical theta gradient by itself; the full simulation still has a cap-stability guard. It controls imposed enthalpy input, not exact discrete total energy, and does not establish thermal equilibrium by construction. Heating-off controls also use the same feedback rule and can therefore experience compensating tendencies from cap relaxation.

The validation suite checks balanced rest, uninterrupted versus restarted fields, controller checkpoint state, tracer passivity and inventory, a 40-day small-box cap/heat test, then both full grids on an H100. Production requires all prerequisites to pass. A probe must project less than 18 wall hours per ten simulated days; production requests 24 hours per segment. These are scheduling gates, not guaranteed mature-flow runtimes.

## Tracer and movies

Three smooth passive patches are initialized identically at every height. Breeze stores their prognostic field as rho C; output divides by reference density to recover dimensionless concentration. It has no dynamical feedback and is not cloud opacity or a JIRAM observable. Its diffusivity is 100 m²/s.

New movies have three depth columns, with vorticity above tracer concentration, on a single synchronized timeline. The requested heights are 20, 60 and 100 km; the actual centers are stored in metadata. Each vorticity depth uses a time-independent scale, while concentration shares a scale across depths. Movies retain every saved frame.

The first unbounded-WENO trial failed the tracer range check (minimum concentration -0.0357). Its dependent production jobs were automatically cancelled. The replacement uses bounded WENO5 for water and dye, SSPRK3, a summed advective CFL target of 0.2, and the extra halo cells required by the limiter. Checks were tightened to numerical-tolerance bounds, rather than clipping concentrations or relaxing acceptance thresholds. The failed attempts and replacement job IDs are retained in the launch record.

Successful morphology alone will not demonstrate an inverse energy cascade. Assess tracer stretching, folding and persistence together with shear evolution, cap stability, transfer diagnostics and resolution sensitivity.

## Validation completed — November 8, 2026

The 40-day heat-budget/tracer test and both full-grid probes passed. [Watch the validation movie and read the measured budgets](https://npizzo.github.io/jupiter-ffr-breeze/#heat-budget). Post-warmup probe costs project 6.39 wall hours per ten simulated days on 1024², and 15.33 hours on 1536². These two-hour probes are preliminary runtime estimates. The first production job is now waiting for resources; it is no longer blocked on validation.
