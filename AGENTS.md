# Simulation work and reporting

The user requires every new simulation-results presentation to include **vorticity-evolution videos at three different depths**. Use actual saved fields; do not infer missing depths from old single-level output. The large-domain campaign saves planes nearest 20, 60 and 100 km above the model bottom, and records the actual cell-center heights. Include all three videos in gallery result entries and user-facing results; mark partial runs and validation runs clearly. Status-only messages without results need not contain videos.

Preserve running experiments: introduce a new driver/configuration for a new campaign, and never change source files being used by active jobs. Keep provenance, source budgets, cap-stability checks and restart validation. Do not infer inverse energy transfer from morphology alone. Raw model data and credentials do not belong in the public repository.

The public gallery landing page is https://npizzo.github.io/jupiter-ffr-breeze/ and is served from this repository (main branch, root index.html, GitHub Pages). Link this page for watching movies, not GitHub blob/README pages or externally hosted sites. New movies must combine three depths in one synchronized 3-column frame, with vorticity on the top row and passive-tracer concentration directly below. Moist validation can add a condensate row. Keep real heights, time stamps, fixed scales and provenance.

Tracer movies use a labeled asinh color scale with linear width 0.01 and physical concentration ticks, fixed across time and all depths, so diluted tracer is visible. Preserve unchanged data and explicit scale provenance; do not silently rescale each frame.

The balanced-heating validation has an explicitly labeled contrast variant requested by the user: tracer asinh centered at C=0.08 with width 0.005, range 0–1, fixed across all times and depths. Production shear keeps the original tracer scale. The contrast variants use symmetric asinh vorticity with width 0.1 in ζ/f0 and full per-depth extrema; retain original renderings and label normalization in provenance. Never use automatic per-frame limits.

## Compute and data locality

Run simulations, diagnostics, field exports, analysis, movie rendering, and media validation on Unity. Do not download model fields or run these computations locally unless strictly necessary. Prefer publication directly from Unity when supported; transfer only the minimal finished publication assets when required by the available GitHub upload mechanism. Lightweight repository editing and browser playback checks are allowed as necessary publication work. This user instruction applies to all future runs and gallery updates.
