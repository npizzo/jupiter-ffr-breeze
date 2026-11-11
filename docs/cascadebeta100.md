# Matched 100-day heated beta-plane experiment

Requested November 11, 2026: repeat the completed heated f-plane experiment from its original initial state, introducing only Jupiter's planetary beta. Computation, diagnostic analysis, rendering and media validation take place on Unity.

The new `cascadebeta100_heated` configuration matches `cascade100_f_heated` for every parameter except the case name, rotation-model selector and beta. It retains the 6000 × 6000 × 120 km channel, 256 × 256 × 64 cells, random 1200-km eddies with seed 31415 and initial Rossby number 0.1, thermal perturbation seed 42, uniform 162 W/m² bottom heating, dry anelastic physics, stable cap, thermal controller and diffusivities. No imposed jet or shear is added. The new run starts at day zero, not from the final f-plane vortex.

At 55°N, `f(y) = f0 + beta*y`, with central `f0 = 2.88095774e-4 s^-1` and `beta = 2.8216700143325537e-12 m^-1 s^-1`, calculated as `2*Omega*cos(latitude)/radius` using the existing Jupiter constants. The channel remains periodic in x and bounded by free-slip meridional walls. The existing f0-based pressure initialization is retained to match the reference initial state; beta introduces its own initial adjustment.

An isolated copy of the original model driver, initializer, thermal controller and transfer diagnostics is hash-checked against the original passing validation record. GPU validation checks checkpoint/restart agreement, runs the full production grid for one day, checks cap stability, tracer conservation and thermal accounting, and compares the actual saved initial fields against the completed f-plane run. A successful throughput test and passing validation record are required for production. Each production segment checks source, configuration and environment hashes.

Production uses ten successive ten-day segments on one L40S GPU, each with four CPU cores, 32 GB RAM and a six-hour limit. Every segment depends on successful completion and diagnostics of its predecessor. Existing experiments are not modified. Outputs remain under `data/cascadebeta100-v1` on Unity.

After day 100, a dependent CPU job renders the two heated cases together, plus separate movies with shared scales. All 801 saved times are included. Three synchronized depth columns show vorticity above passive tracer, near 20, 60 and 100 km. Vorticity uses a fixed symmetric asinh scale at each depth shared across the two cases, normalized by the same central f0; tracer uses the fixed asinh width 0.01 and physical concentration ticks. Full decoding and frame/time checks are required before publication to the [GitHub Pages gallery](https://npizzo.github.io/jupiter-ffr-breeze/).

The comparison tests the effect of beta conditional on this coarse, strongly heated dry model. It does not establish resolution convergence or an inverse cascade from morphology. Existing resolved-flow flux and Reynolds-stress-work diagnostics are retained; their budget does not close all beta, numerical or boundary contributions.

[Submission and job dependencies](cascadebeta100-launch.json). Submission is not a passing validation result or a completed simulation.
