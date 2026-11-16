# Simulation work and reporting

The user requires every new simulation-results presentation to include **vorticity-evolution videos at three different depths**. Use actual saved fields; do not infer missing depths from old single-level output. The large-domain campaign saves planes nearest 20, 60 and 100 km above the model bottom, and records the actual cell-center heights. Include all three videos in gallery result entries and user-facing results; mark partial runs and validation runs clearly. Status-only messages without results need not contain videos.

Preserve running experiments: introduce a new driver/configuration for a new campaign, and never change source files being used by active jobs. Keep provenance, source budgets, cap-stability checks and restart validation. Do not infer inverse energy transfer from morphology alone. Raw model data and credentials do not belong in the public repository.

The public gallery landing page is https://npizzo.github.io/jupiter-ffr-breeze/ and is served from this repository (main branch, root index.html, GitHub Pages). Link this page for watching movies, not GitHub blob/README pages or externally hosted sites. New movies must combine three depths in one synchronized 3-column frame, with vorticity on the top row and passive-tracer concentration directly below. Moist validation can add a condensate row. Keep real heights, time stamps, fixed scales and provenance.

Tracer movies use a labeled asinh color scale with linear width 0.01 and physical concentration ticks, fixed across time and all depths, so diluted tracer is visible. Preserve unchanged data and explicit scale provenance; do not silently rescale each frame.

The balanced-heating validation has an explicitly labeled contrast variant requested by the user: tracer asinh centered at C=0.08 with width 0.005, range 0–1, fixed across all times and depths. Production shear keeps the original tracer scale. The contrast variants use symmetric asinh vorticity with width 0.1 in ζ/f0 and full per-depth extrema; retain original renderings and label normalization in provenance. Never use automatic per-frame limits.

## Compute and data locality

Run simulations, diagnostics, field exports, analysis, movie rendering, and media validation on Unity. Do not download model fields or run these computations locally unless strictly necessary. Prefer publication directly from Unity when supported; transfer only the minimal finished publication assets when required by the available GitHub upload mechanism. Lightweight repository editing and browser playback checks are allowed as necessary publication work. This user instruction applies to all future runs and gallery updates.

## Historical single-height exception

On November 11, 2026, the user explicitly requested clearly labeled single-height movies for older campaigns that saved only one vorticity level. Show the actual saved height, state that other depths and tracer were not saved, and never fabricate missing panels. This exception applies to historical output only; retain the synchronized three-depth vorticity/tracer layout for new simulation output.

## Full 100-day tracer contrast variant

On November 11, 2026 the user requested clearer late-time tracer gradients in the completed 100-day movies. The labeled contrast variant uses unchanged absolute concentration, BrBG colors, fixed center C=0.0793082, and two-sided asinh width 0.001 with the full concentration range. The same mapping applies to both cases, all depths and all times; no framewise normalization, mean subtraction or clipping. Preserve original movies and original vorticity normalization. Do not modify frozen model or renderer source used by queued campaigns; add separate rendering variants when needed.


## Explicit assumptions and viscosity-only 2D tests

Before launching simulations, explicitly state the assumptions to the user in short bullets: equations and dimensionality, geometry, initial conditions, forcing, dissipation, resolution, duration, controls and any release gates. Do not infer authorization for one physical term from authorization for another. This disclosure does not itself require a new permission question when the user has already authorized the experiment.

For the current minimal 2D cascade/LSV experiments, use ordinary Laplacian viscosity only. The user explicitly rejected cubic drag and did not authorize importing van Kan-style dissipation along with the instability forcing. Do not add linear drag, nonlinear drag or hyperviscosity without an explicit change in the user's instructions. Compare unforced decay, weak instability forcing, weak random-sign stirring and weak f-aligned-core stirring from identical initial conditions. Verify transient unforced inverse energy transfer using native flux and absolute low-wavenumber energy, not morphology or energy fractions alone. Keep this separate from beta-plane Ferrel-cell work.


## Long 512² two-case campaign

On November 16 the user authorized the higher-resolution forced/unforced long runs and said no comparison movies or diagnostics are needed. For `viscous2d_long512_v1`, run only the instability-forced and unforced cases through t=200, with ordinary viscosity and no drag. Do not submit movie or diagnostic-analysis jobs for this campaign. Retain checkpoints and raw vorticity; basic finite-state/CFL and implementation/restart checks are retained. This explicit two-case request supersedes the earlier four-case plan for this campaign.
