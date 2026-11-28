# Marcus continuation and horizontal-resolution check

Authorized28November2026. Source/model/reference copied from the completed four-day CV1 reconstruction. No physical change. Baseline384x384x256 resumes its exact day4 native checkpoint and stops at day12 for review, with day40 the intended subsequent milestone. Fine512x512x256 starts from the same analytic initial condition, runs a1800s full-grid pilot, then resumes to day2 only if it passes and projects under16GPUhours. This changes horizontal resolution only.

Same dry anelastic model, tilted rotation23S, beta0, prescribed background shear, reference stratification, impermeable lids,31km vertical and7000km meridional absorption, tau20000s, WENO9, no closure or heating. No automatic extension to40or120days.

Derived fields refresh after checkpoint restore before diagnostics. Native timestep callback is restored. CPU and GPU real-driver tests compare continuous and restarted trajectories, validate all added outputs and enforce source fingerprints. Live scheduler deadline controls checkpoint margin. Sequential chains restart only normal wall stops; failures stop without retry. Full-grid fine pilot precedes production. Existing jobs/source unchanged.

Additional outputs: horizontal velocity on the same three saved planes and fixed central x-z vorticity/theta/T/w sections. Original absolute vorticity/tracer retained. No extra full-volume series. Existing fixed-ellipse circulation remains explicitly a fixed-aperture diagnostic. Later analysis should supplement it with moving-core morphology and track sponge power; do not claim an exact energy budget.

Measured baseline cost8h27m/fourdays gives ~17additionalA100hours to day12 and~76to day40. Fine two-day estimate~10hours is provisional. Review coherence, hollowing, thermal structure, circulation and boundary absorption before further extension. These are initialized-vortex persistence experiments, not spontaneous sign selection.
