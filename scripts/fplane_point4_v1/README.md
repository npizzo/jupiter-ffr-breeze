# F-plane depth-mean-vorticity sign-selection tests

These are hypothesis tests requested after the 35°N horizontal-resolution audit. All simulation, saved-field diagnostics and rendering run on Unity. Existing completed data and frozen running sources are unchanged. No claim of a third invariant or Ferrel-cell detection follows from vortex morphology or a fitted curve alone.

## Bounded initial campaign

| Cases | Interval (Earth days) | Purpose |
|---|---:|---|
| early_replay_heated, early_replay_control | 0–5 | Recreate the original random-turbulence initial conditions and original timestep ceiling to connect early skewness with the actual completed day-100 endpoints. |
| late_heated, late_control | 90–100 | Restart the actual completed cascade100 checkpoints to save full 3D velocities at matched daily times, with the original 30 s maximum timestep. |
| rest_q0, rest_q5p4, rest_q54, rest_q162, rest_q540 | 0–5 | Start at rest with identical weak physical-mode temperature noise; vary only bottom heat flux. |
| rest_q162_flipf | 0–5 | Reverse f at fixed magnitude; same weak-noise physical initial condition. A cyclone then has negative vorticity. |
| rest_q162_fine | 0–5 | 512²×64 companion to the 256²×64 162 W/m² case; same physical noise modes and 15 s timestep ceiling. |
| late_restart_control, late_restart_reversed | 100–110 | Restart the completed day-100 heated state unchanged or with only its density-weighted barotropic horizontal velocity reversed. |

The from-rest cases isolate convective injection; they cannot alone predict the existing turbulence-seeded condensate. The additional early replays supply the matching initial conditions for that comparison.

The short runs are an onset and source-budget screen. Five days need not be enough for a condensate, and ten days need not erase a reversed mature vortex. An absence within these windows is not evidence against eventual sign selection. The flux sweep is a sensitivity study, not four calibrated Jovian atmospheres. The thermal cap, volume heat-balance correction and upper vertical-velocity sponge remain in all cases, including zero bottom flux.

## Data and conservation checks

Original native Float64 staggered velocity, centered potential temperature and dye mass density are saved every 6 hours (daily for long replay/reversal pairs), with reference density, column mass, actual time, grid spacing and f. Three-depth vorticity/tracer movies retain 3-hour sampling. The start and end states are always saved. No source is reconstructed from images. Scientific files stay on Unity; only final media, compact reports and source code are published.

The reversal transforms horizontal momentum by subtracting twice its density-weighted depth-mean velocity times reference density. It preserves vertical momentum, thermodynamics, tracer and the baroclinic horizontal component. The preflight verifies native mass divergence and wall fluxes, and records the original pressure as an initial guess for subsequent projection. A matched unchanged restart separates normal evolution from the reversal response.

Six GPU preflights cover both day-90 restarts, actual day-100 reversal, 540 W/m², negative f, and 512². Analytic CPU checks on Unity validate the channel Poisson solver, signs, moment normalization, quadratic coefficient conversion and zero baroclinic sources for a barotropic field. Production requires both gates. Sources are hashed at the gate; frozen code changes invalidate it.

## Diagnostics and interpretation

- Z is the reference-density-weighted vertical mean of vertical vorticity. Full-channel, interior (outer 10% removed at each wall), wall-strip, zonal and nonzonal statistics are reported separately.
- Streamfunction uses a periodic-x/channel-y Poisson solve. Equal wall values define zero transport; the measured uniform transport is restored as a linear streamfunction component, then its spatial mean is set to zero. The nonzonal inversion is reported separately.
- Quadratic and identifiable logistic fits give R², raw c/b and standardized curvature, pooled across equal-time samples and separately at each snapshot. Nonlinear fitting is deterministically thinned to at most 3,000 points; all retained pooled samples score the fit. Snapshot variation supplies a temporal sensitivity check. These are descriptive fits, not independent-pixel uncertainty estimates or evidence of equilibrium.
- In the proposed two-level form, occupancy a and offset μ are structurally degenerate. Only the combined logistic intercept is fitted. An occupancy inferred from the mean is labeled as such, and compared cautiously with explicitly thresholded positive/cyclonic core areas. Upper-vorticity-level occupancy is not automatically cyclonic occupancy when f is negative.
- Raw and central skewness, PDFs and sign fractions are recorded for Z and for 3D vorticity. Fu & Sun's S≈2.5 Ro relation concerns volumetric 3D vorticity in rotating Boussinesq convection, not automatically Z or an anelastic capped channel. Both uniform-volume and reference-density-weighted 3D results are retained.
- Smoothing widths are Gaussian FWHM 50, 100 and 200 km, with odd curl extension at free-slip y walls and periodic x. The 50 km filter spans only 2.13 coarse cells and is explicitly flagged as unresolved under a four-cell criterion. The 200 km split means Gaussian-low-pass plus its exact complementary remainder, not a sharp Fourier cut.
- V = −2〈ζ′δ′〉 and H = −〈u′·∇hζ′〉 + 〈ωh·∇h w〉 are evaluated after full 3D differentiation, with density weighting. Signs, spatial integrals, convection-conditioned values and 200 km core-ring values are measured, not imposed. The through-origin V regression against 2(f+Z)〈δ′²〉 gives an effective τ; it is not assumed to be a grid-dissipation time.
- The contribution to the third moment is 3〈Z² source〉. Finite-time storage is compared with the sampled centered source reconstruction; all unclosed remainders are retained. They include viscosity, native WENO differences, temporal sampling and omitted terms and cannot be relabeled as physical dissipation. Daily or 6-hour source snapshots may undersample fast convection.

The comparison therefore tests whether the observations are consistent with a skewed injection mechanism. It cannot establish a conserved third invariant for the forced, dissipative, three-dimensional atmosphere.

## Source reference

Fu & Sun (2024), *Vorticity skewness of finite-amplitude rapidly rotating Rayleigh–Bénard convection*, JFM 997 A3, [doi:10.1017/jfm.2024.571](https://doi.org/10.1017/jfm.2024.571).
