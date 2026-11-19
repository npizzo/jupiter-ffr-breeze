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

On November 16 the user authorized the higher-resolution forced/unforced long runs. They clarified that they DO want the 512² movies, just no comparison to 256². For `viscous2d_long512_v1`, run only the instability-forced and unforced cases through t=200, with ordinary viscosity and no drag. Render individual vorticity movies and a synchronized forced/unforced movie on Unity and publish them to the GitHub gallery. No additional scientific diagnostics are requested; basic finite-state/CFL, implementation/restart and media checks remain. These genuinely horizontal 2D models have one plane and no tracer: label this clearly rather than fabricating depths. This explicit two-case request supersedes the earlier four-case plan for this campaign. Preserve frozen simulation source and original configuration metadata; the original output-description string predates the movie clarification.

## Efficient continuation to t=400

The user subsequently authorized both cases to continue from their own t=200 checkpoints to t=400 with a substantially more efficient implementation. `viscous2d_fast512_v1` preserves all physical parameters and dt=0.0025, uses native viscosity, cached in-place forcing FFTs and specialized array loops. The Unity release check must pass equation equivalence, serialized restart and measured speed/allocation criteria before production. Do not resume the allocation-heavy prototype for further long runs. Preserve the original t≤200 output and source. Continue the same 2D vorticity movie convention, without 256² comparisons or additional scientific-analysis jobs.



## Explicit changes before full campaigns

Before submitting or releasing a production campaign, give the user concise bullets stating every change from the named reference case: domain, grid/spacing, initial velocity RMS and eddy wavelength, Rossby definition/value, seed/mode selection, forcing, boundary/cap treatment, duration and total compute limit. Distinguish unchanged controls from changed parameters and numerical stability checks from evidence of physical fidelity. Do not silently resize the domain or seed scale to save cost, and do not expand a single-case authorization into a parameter sweep. On November 18 the user authorized one no-shear 50-day test restoring the 6000-km box and 1200-km eddies at Ro=0.1 with approximately 5.50 m/s initial RMS, allowing 128 x 128 x 64 resolution for speed; other current physical settings and odd-mode seed selection remain unchanged. The superseded Ro=0.3 setup was held before allocation and must not also be run.


## Authorized parallel shear extension

Later on November 18, the user explicitly authorized adding the cyclonic and anticyclonic initial-shear cases concurrently with the revised fixed-Ro no-shear test. This supersedes the earlier single-case restriction for this campaign only. Preserve the 6000-km box, 1200-km random eddies, random RMS about 5.50 m/s, random Ro=0.1, 128 x 128 x 64 grid and all other current control physics. Add only initial background shear with peak zeta/f=+0.05 or -0.05; disclose that this means peak background speed about 27.5 m/s in the wider domain, separate from random RMS. The three independent jobs each have a one-day stability/cost gate and three-GPU-hour limit, nine GPU-hours combined. No automatic extension or further parameter sweep is authorized. Render the three synchronized three-depth vorticity/tracer movies on Unity after all cases finish; publish only verified completed assets to GitHub.


## Authorized 256-squared shear resolution reruns

On November 18 the user authorized rerunning the three fixed-Ro large-box cases under the same physical conditions at higher resolution, as soon as possible. `strong256_v1` changes nx and ny from 128 to 256, retaining nz=64 and all physical settings. Explicitly disclosed operational changes: bulk diagnostics and cap-stability checks every 3 simulated hours instead of 30 minutes; retain 3-hour movie frames, daily full fields, 2-day checkpoints, 10-minute thermal controller, CFL and NaN checks. Allow one A100, H100 or L40S per independently queued job. A one-day timing/stability gate precedes each continuation to day 50; limit six GPU-hours per case, eighteen total, with no automatic expansion or retry. A Unity CPU renderer follows all three completions. These initial gates do not establish late-time spatial fidelity: compare the near-grid stripe wavelengths and sign statistics against the completed coarse cases before treating the result as converged.


## Authorized A100/H100 checkpoint continuation

After all three `strong256_v1` L40S runs completed day 1 but stopped at the cost gate, the user explicitly authorized continuing the retained checkpoints on A100/H100. `strong256_resume_v1` changes no physical or numerical configuration. Advance day 1 to day 2 as an actual-hardware timing gate, verify restart fields, then continue to day 50 only when the forecast with 25% margin fits. The new allocation is 21000 seconds per case; including the prior 568/564/553 seconds, this stays below the existing six-GPU-hour per-case and eighteen-GPU-hour total limits. Preserve previous day-1 output; do not regenerate initial conditions or duplicate completed integration. Replace the cancelled movie dependency with a job assembling days 0-1, 1-2 and 2-50. Further failed gates must retain output for review rather than trigger automatic retries.


## Direct continuation supersedes the intermediate timing gate

The user explicitly requested running these cases straight to day 50 rather than staging another timing gate. `strong256_direct_v1` supersedes the unstarted `strong256_resume_v1` jobs, which were held and cancelled with zero GPU time. Resume the original saved day-1 checkpoints directly to day 50 on A100/H100, with no day-2 forecast-based stop. All physical settings, resolution, output cadence, normal stability checks and periodic checkpoints remain unchanged. Keep the existing 5h50m continuation allocation (prior plus new stays within six GPU-hours per case), and reconnect the dependent three-depth movie renderer to the direct jobs. Do not reinstate the intermediate runtime gate without an explicit reason and user direction.
