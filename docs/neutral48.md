# Neutral-layer convection test

The previous controlled runs remained dominated by their initially unstable lower layer. This experiment sets the initial horizontal-mean lower-layer buoyancy frequency squared to exactly zero, while retaining the stable cap above 80 km. Zero-horizontal-mean temperature noise with amplitude parameter 10⁻⁵ K seeds symmetry breaking. That parameter is 2,000 times smaller than the previous 0.02 K case; the actual RMS fluctuation is smaller still. The initial noise remains a small perturbation energy source, so a matched unforced control is essential.

All other physical settings are retained: dry anelastic Breeze, 55°N Jovian f-plane, beta = 0, no imposed shear, 1,200 × 1,200 × 120 km, viscosity 1,000 m²/s, and thermal diffusivity tapered from 1,000 to 1 m²/s across 60–80 km. Changing the neutral reference profile also changes its consistent hydrostatic pressure and density profiles. The physical heat flux remains 5.4 W/m²; it has not been accelerated by increasing the forcing.

| Case | Bottom input [W/m²] | Upper input [W/m²] | Grid |
|---|---:|---:|---|
| neutral48_heated | +5.4 | −5.4 | 128 × 128 × 64 |
| neutral48_cooled | 0 | −5.4 | 128 × 128 × 64 |
| neutral48_unforced | 0 | 0 | 128 × 128 × 64 |
| neutral48_refined | +5.4 | −5.4 | 128 × 128 × 128 |

All runs last 48 Earth hours. The standard cases have identical initial conditions. The refined run checks sensitivity to resolving the heating layer, but has a different spatial realization of the same seeded-noise distribution. It is not a complete convergence study.

The existing finite-field and thermal-source checks apply to every run, together with verification that the initial lower-layer mean potential temperature is uniform. The comparison will examine vertical velocity, kinetic energy, resolved heat transport, the emerging mean stratification, and dependence on forcing and resolution. A difference from the unforced run alone does not demonstrate developed convection: absolute amplitude, growth or decay, vertical extent, and heat transport matter.

The 48-hour duration is an initial diagnostic window, not a claim that the weak forcing reaches equilibrium or even resolvable onset on this grid. If no sustained motion develops, interpret that in the context of rotation, imposed diffusivity, resolved heating-layer thickness, and duration; it would not establish that Breeze cannot simulate Jovian convection.

Submit `sbatch slurm/neutral48.sbatch`. Four runs execute sequentially on one L40S, with four CPUs, 8 GB host memory, and a one-hour limit. Outputs go to the new directory `data/neutral48-v1`. Previous results are preserved. No new beta, shear, or model choice is introduced.

Requested snapshot output includes vertical velocity and departures from the horizontal-mean potential temperature at 6, 12, 24, and 48 hours for every case. Final comparisons use common color scales across cases. Time sequences use fixed scales within each sequence; color bars retain actual physical units.
