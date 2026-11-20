# Bounded anelastic convection pilot

This is a new, exploratory Breeze configuration for finding an affordable convecting parent state before testing cyclone/anticyclone amplification. It is not a reproduction of Käpylä et al. or Chan & Mayr, and does not establish an anticyclonic regime in advance. The older long compressible campaigns remain held or cancelled.

## Disclosed physical assumptions

- Dry, three-dimensional anelastic dynamics; constant positive `f=0.5`; no beta, imposed shear, vortex seed, moisture or momentum drag.
- Periodic horizontal box `2 × 2 × 1`, `128 × 128 × 64` cells, Float64. Impermeable, stress-free top and bottom.
- Adiabatic hydrostatic reference: `g=0.8`, `cp=1`, `R=0.4`, base temperature and density both 1. The reference temperature is `1−0.8z` and density `(1−0.8z)^1.5`, giving approximately 11.18:1 bottom/top density contrast. The anelastic constraint is divergence of reference-density-weighted velocity = 0.
- Bottom upward potential-temperature-density flux `9.375e-7`; zero top diffusive thermal flux. In the upper 20%, smooth potential-temperature relaxation toward a stable target removes heat. This is a penetrative-convection idealization, not a fixed-temperature top boundary or the Chan–Mayr thermal prescription.
- Native density-weighted Newtonian viscosity `nu=1e-5`, potential-temperature diffusivity `kappa=2.5e-5`, WENO5 transport. No subgrid model. Viscous heating is omitted: a conservative potential-temperature balance is checked, **not a closed total-energy budget**.
- Initial velocity is zero. The lower layer has a 3% unstable potential-temperature contrast and small, horizontally mean-free thermal noise. The upper layer's stable target increases potential temperature by 1%. The initial unstable profile contains available potential energy, so startup growth must not be mistaken for a statistically sustained convective state.
- Passive dye is initialized smoothly at concentrations 0.05–0.15. Three physical heights nearest `z=(0.2,0.6,0.9)` are saved for synchronized vorticity/tracer movies.
- Time is nondimensional; these are not Earth days or a dimensional Jupiter calibration. Rotation dominance and convection scale must be measured rather than inferred from the value of f alone.

## Cost and implementation gates

CPU preflight job 2029786 passed resting-atmosphere, bottom-flux conservation, mass projection and short heated-integration checks. The uniform-heating potential-temperature budget error was below `1.1e-9` of the imposed input at the sampled times. These are implementation checks, not turbulence validation.

One GPU allocation is limited to **one hour total**, including GPU compilation/checks and parent development. It uses an FP64-capable Unity GPU; the actual device is recorded. It first checks rest and serialized restart agreement, then times the actual 128²×64 grid through t=5. It stops if the remaining budget cannot plausibly cover another 100 time units with a 50% timing margin. The parent has an upper target t=2000, but the wall-clock cap takes precedence. There is no automatic extension, production release or vortex-seeded branch.

Runtime checks stop on nonfinite fields, appreciable wall penetration, mass-flux divergence above 1e-8, potential-temperature departures above 10%, Mach number above 0.3, or sampled thermal-budget discrepancy above 5% of bottom input after t=5. The sampled thermal-budget check uses a trapezoidal source integral; a failure requires distinguishing sampling error from a discretization defect before changing tolerances.

## Scientific decision after the pilot

Review the measured convective flux and thermal drift, depth and gradient Rossby numbers, Reynolds number, accumulated turnover count, vertical velocity scale and horizontal resolution, cap stability, and turbulent throughput. A stable run is not automatically an adequate parent state. In particular:

1. Sustained convective transport must replace release of the initial unstable profile. A short growing transient is insufficient.
2. Convective motions must be resolved. The recorded gradient-derived wavelength in grid cells is an initial screen, not a convergence certificate; inspect spectra/fields before releasing branches.
3. The stable layer must remain stable, and reference-state departures must remain small enough for anelastic interpretation.
4. The measured cost per turnover must make a matched response experiment affordable. If this budget cannot provide a useful state, report that outcome before committing more resources.

If these checks succeed, the proposed next experiment is a matched positive/negative compensated vortex pair in separate copies of the same parent, plus an unseeded control, with equal seed size/energy, balanced initialization, and no sustained seed forcing. Compare convection's contribution to vortex growth and maintenance over measured turnover times. This would test finite-amplitude support; it would not by itself prove spontaneous sign selection. A later seed-amplitude reduction and unseeded confirmation would be required. Those branches are not submitted by this pilot.

All scientific computation, diagnostics and rendering remain on Unity. Source/configuration and compact provenance are public; raw fields and checkpoints are not.
