# Controlled convection experiment

The objective is to determine whether the present Breeze setup sustains convection through its imposed bottom heat flux, as a first step toward Jupiter-like FFRs. There is no comparison with another solver.

## Cases

All four runs use 128 × 128 × 64 cells, a 1,200 × 1,200 × 120-km box, the pilot's weakly unstable lower layer and stable cap, and seed 42. Duration is 172,800 s (48 Earth hours, about 4.8 Jovian rotations). Bottom heating is the only difference within each pair.

| Case | Bottom enthalpy flux | Initial perturbation parameter |
|---|---:|---:|
| convection48_on_p020 | 5.4 W/m² | 0.2 K |
| convection48_off_p020 | 0 W/m² | 0.2 K |
| convection48_on_p002 | 5.4 W/m² | 0.02 K |
| convection48_off_p002 | 0 W/m² | 0.02 K |

The random field is mean-subtracted independently at every height and tapered vertically. The quoted parameters multiply that field; they are not its RMS or maximum temperature amplitude.

## Rotation and geometry

These are traditional f-plane runs with f = 2Ω sin(55°) = 2.88096 × 10⁻⁴ s⁻¹ and Ω = 1.7585 × 10⁻⁴ s⁻¹. Beta is zero; horizontal components of the planetary rotation vector are omitted. The side boundaries are periodic. A future beta-plane experiment requires a consistent meridional geometry and boundaries; simply adding a linear f(y) across a doubly periodic box would introduce a seam.

The initial convective Rossby number sqrt(−N²)/f ≈ 1.1 is not a strongly rotation-dominated limit. The measured relative-vorticity RMS divided by f is recorded separately as an evolving flow diagnostic.

## Diagnostics

Every 300 simulated seconds, `diagnostics.csv` records face-weighted column kinetic energy, cell-centered vertical/horizontal RMS speed, maximum vertical speed, density-weighted potential-temperature content, imposed bottom potential-temperature flux, signed upper-layer relaxation rates, sponge work, approximate buoyancy work, relative-vorticity RMS/f, and horizontal potential-temperature variance. Fields and profiles are saved hourly.

The conservative thermal diagnostic is Θ = integral rho_ref × horizontal_mean(theta) dz, in kg K/m². Its expected tendency is the imposed rho-theta boundary flux plus the diagnosed relaxation source. The bottom flux equals 5.4/(cp × Exner_bottom) in heated cases. Integrating the sampled relaxation rate with the trapezoidal rule introduces a quadrature error; this is reported as a residual, not silently interpreted as a solver conservation error. Θ is not total thermal energy. The signed relaxation enthalpy rate is reported separately in W/m².

Resolved convective enthalpy-flux profiles use rho_ref cp Exner × mean(w′ theta′) at cell centers. Diffusive enthalpy-flux profiles use −rho_ref cp Exner κ d(mean(theta))/dz with finite differences. These are diagnostic approximations, not the exact numerical transport flux. Kinetic energy, buoyancy work, and sponge work do not constitute a closed discrete kinetic-energy budget: numerical/viscous losses and collocation errors are not independently resolved.

The upper relaxation can heat or cool; it is not assumed to be a cooling-only sink. Explicit diffusion can redistribute heat from the stable cap even in the bottom-heating-off controls. Both effects must be examined when attributing the flow to bottom heating.

## Execution

Submit `sbatch slurm/convection48.sbatch` from the Unity project directory. One GPU runs the four cases sequentially in one Julia process, reusing compilation. Fresh output is written to `data/convection48-v1`; existing nonempty outputs are never overwritten. The campaign first runs hydrostatic-rest and heating-only checks of the new diagnostics and stops on a failed check.

The first 48 hours are a diagnostic window, not an assumed thermal equilibrium. Interpretation will compare late-time behavior and forcing budgets, including dependence on the initial perturbations. A visually interesting transient alone is insufficient evidence for sustained FFR dynamics.
