# Controlled 48-hour convection results

All four runs completed on Unity and passed finite-state and diagnostic checks. They show an initial convective burst followed by strong decay. This setup has not demonstrated sustained bottom-heated FFR dynamics.

## Paired comparison

| Initial perturbation parameter [K] | Bottom heating | Mean vertical RMS speed, hours 24–48 [m/s] | Vertical RMS speed at 48 h [m/s] | Mean RMS vorticity/f, hours 24–48 |
|---:|---|---:|---:|---:|
| 0.2 | on | 0.08052 | 0.03546 | 0.1246 |
| 0.2 | off | 0.07949 | 0.03420 | 0.1238 |
| 0.02 | on | 0.08773 | 0.03739 | 0.1076 |
| 0.02 | off | 0.08735 | 0.03693 | 0.1068 |

Turning bottom heating on increases the second-day mean vertical RMS speed by **1.30%** for the 0.2-K perturbation parameter and **0.44%** for the 0.02-K parameter. Mean kinetic-energy changes are approximately 0.99% and 1.31%, respectively. These are matched deterministic runs with one seed, not an ensemble estimate of a small heating response. The smaller perturbation mainly delays the convective burst; both amplitudes subsequently decay.

## Why the current forcing needs revision

The upper Newtonian relaxation acts as a **heat source**, adding about **56 W/m²** on average during hours 24–48 in all four cases. Bottom heating supplies only 5.4 W/m² when enabled. The relaxation replenishes heat lost from the initial stable cap as constant thermal diffusion smooths its temperature gradient. This is not the intended radiative heat sink.

The explicit thermal diffusivity is 1,000 m²/s throughout the domain. At 48 hours, the diagnostic downward diffusive heat flux reaches about **491 W/m²** within the cap. This is internal redistribution, not an additional global energy source. The mean profile shows the stable region extending downward substantially, while the remaining vertical motion becomes weak and confined to lower levels. In the heated 0.2-K case, the largest resolved upward convective heat flux at 48 hours is only about **0.25 W/m²**.

The heating-off controls evolve almost identically to the heated cases. Together with the initial superadiabatic profile, this supports attributing the early strong convection mainly to stored initial instability. The 48-hour window does not establish a final thermal equilibrium and does not show that Breeze is incapable of sustained convection or FFRs.

## Rotation

Every run uses a **traditional f-plane**, not a beta plane: Ω = 1.7585 × 10⁻⁴ s⁻¹, latitude 55°N, f = 2.88096 × 10⁻⁴ s⁻¹, and beta = 0. The rotation period is about 9.93 hours. The initial convective Rossby estimate is 1.10. Measured volume-RMS relative vorticity/f falls to roughly 0.11–0.12 on average during the second day, indicating stronger rotational constraint on the remaining flow; this does not make the initial burst an asymptotic rapidly rotating regime.

## Validation and cost

The new diagnostics passed separate hydrostatic-rest and heating-only tests before the campaign. Matched pairs share the grid, initial stratification, seed, rotation, transport coefficients, and solver source hash. The density-weighted potential-temperature inventory agrees with integrated boundary and relaxation sources to within **4.4 × 10⁻⁷ relative** across the cases. That residual includes five-minute source-integration quadrature error. The inventory is not total energy, and the kinetic-energy diagnostic is not a closed discrete energy budget.

All four 48-hour runs plus short validation cases completed in **14 min 48 s** on one NVIDIA L40S. Individual measured run intervals were about 179–188 s, including diagnostics and output. Slurm recorded 3.55 GiB maximum host RSS. A direct GPU-memory sample was 849 MiB; sampled values can miss brief peaks and are not guaranteed capacity requirements.

## Next model change

Replace the reference-profile relaxation with an explicitly controlled upper heat sink and account for its net input/output. Reduce excessive stable-layer mixing using a justified diffusivity profile or stability-dependent closure, with resolution checks. Then repeat the source-control experiment before adding beta, larger domains, cyclonic shear, or moist storms. Keep Jupiter’s rotation fixed. These changes address the experiment’s identified forcing and transport limitations; morphology alone should not be used to tune them.

[Comparison gallery](../gallery/convection48/README.md) · [Experiment design](convection48.md) · [Numerical summary](../gallery/convection48/summary.json).
