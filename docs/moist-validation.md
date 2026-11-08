# Moist-convection validation in hydrogen–helium

[Gallery](https://npizzo.github.io/jupiter-ffr-breeze/) · [Launch record](shear20-launch.json)

These are implementation-validation experiments, separate from the dry FFR campaign. Breeze is pinned to 0479d9c8ce006adcb1e3314c324b84980ad4871b. Its mixed-phase saturation adjustment represents water vapor, liquid and ice with latent heating. The background gas uses molar mass 0.0023 kg/mol, cp = 13000 J/(kg K), and gravity 25 m/s². Water vapor is heavier than this carrier gas: importing an Earth light-vapor buoyancy approximation would be incorrect.

## Tests

1. **Thermodynamic identities.** Independent gas-mixture density and condensate-loading checks; a triple-point saturation check; and 36 combinations of pressure, temperature and water fraction. Check phase positivity, total water, saturation and independently reconstructed static energy before and after adjustment. These tests do not simply compare a routine to itself.
2. **Prescribed ascent.** Dry and water-bearing states evaluated along a specified pressure/height path, conserving the defined static energy and allowing reversible phase partition. This is a thermodynamic path test, not a simulated rising parcel or an independent moist-adiabatic benchmark.
3. **Resting moist atmosphere.** A 32³ periodic-horizontal, impermeable-top/bottom anelastic box, integrated for one hour. Check small vertical velocity, finite fields, phase closure and water inventory.
4. **Matched dry and moist plumes.** An 80 × 80 × 40 km box at 64³; a Gaussian bubble centered at 8 km with 5 km width and 2 K dry-equivalent buoyancy. The moist bubble has peak total-water mass fraction 0.020001. Initial temperature compensates for the mixture gas constant so initial density buoyancy matches the dry case, provided there is no initial condensate; that initial condensate is explicitly recorded. Bottom reference conditions are 5 bar and 300 K. Rotation uses the Jovian f at 55°N. There is no external heat forcing during this one-hour experiment.
5. **Resolution sensitivity.** Repeat the moist plume at 96³. Compare the saved cloud-water history, vertical motion and structure with 64³. Passing inventory and positivity checks at both grids does not by itself demonstrate convergence.

All dynamic cases check total-water and passive-tracer inventories, finite fields and phase closure. Plumes must produce condensate; the rest test must remain quiet. Three-depth output is saved every two minutes near 4, 12 and 24 km, with actual heights in metadata. Movies show vorticity, passive concentration and, for moist plumes, liquid-plus-ice concentration in aligned rows. Concentration is not infrared radiance.

## Current validation status and numerical correction

The first thermodynamic matrix passed 220 checks and the prescribed-ascent suite passed 406 checks. The initial dynamic trial failed when unbounded scalar advection produced a negative total-water fraction (-3.92e-7) in the moist plume. This is a validation failure, not accepted moist-convection output.

The replacement uses the pinned library's bounds-preserving WENO5 on water and dye, SSPRK3, a summed advective CFL target of 0.2 and enlarged halos required by that scheme. Negative values are not manually clipped and the original water tolerance remains in force. The revised dynamic suite completed successfully (job 66957583). All four dynamic cases passed their stated checks. This clears these implementation checks, not the broader validation required for a moist FFR production run.

This validates selected implementation properties, not an observed Jovian storm. There is no precipitation removal, ammonia chemistry, entrainment calibration, measured Jovian moisture profile, or radiative transfer here. Strong moist-plume speeds will also require checking the anelastic approximation's applicability. Follow-on validation should add an independently specified moist benchmark and quantitative refinement criteria before interpreting Jovian FFR moisture feedback.

## Reproduce on Unity

Use the existing pinned environment and submit `slurm/moist-validation.sbatch` with a new, empty output root. It calls `scripts/run_moist_validation.jl`, including the 0D tests, and writes `PASS.toml` only after every requested dynamic case passes. The current root and job ID are recorded in `docs/shear20-launch.json`.

## Completed first comparison

[Watch synchronized three-depth movies](https://npizzo.github.io/jupiter-ffr-breeze/#moist-validation). Peak domain-wide vertical RMS is 6.76 m/s in the dry 64³ case, 7.99 m/s in the moist 64³ case and 8.03 m/s in the moist 96³ case. Peak condensate mass fraction increases from 0.01537 to 0.01605 with refinement. Water and tracer relative inventory drift stays below 4e-14 in the recorded histories.

The maximum sampled vertical speed changes from 96.1 to 142.6 m/s between the moist grids (about 48%), although peak RMS changes by only 0.48%. Local extremes are therefore not resolution independent. This short, periodic plume test also requires domain/boundary sensitivity checks before physical interpretation. It does not establish convergence or Jupiter FFR formation.
