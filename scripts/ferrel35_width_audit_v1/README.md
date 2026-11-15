# Matched channel-width audit, days 50–100

Read-only analysis on Unity of the completed 35°N narrow and wider channels. No simulation is modified or restarted. Run `tools/point4-analysis-env/bin/python scripts/ferrel35_width_audit_v1/analyze.py` from the Unity project root; the output directory must not already exist.

The channels are 6,000 × 6,000 × 120 km and 6,000 × 12,000 × 120 km, respectively, on 256 × 256 × 64 and 256 × 512 × 64 grids. The analysis checks equal grid spacing and matched forcing, rotation at shared coordinates, mixing, reference atmosphere, cap, and initialization parameters. The realized random channel modes differ, so this is one matched-statistics realization at each width, not an identical-field pair.

## Predefined comparison

- Analyze all 401 actual three-hour snapshots from Earth days 50 through 100, with trapezoidal time averages.
- Compare the same central physical region: |y| ≤ 2,400 km and 10 ≤ z ≤ 75 km. Use reference-density-weighted norms, without spatial smoothing.
- Separately measure north and south wall strips within 600 km of each wall. Plot RMS profiles against distance from each wall, using fixed 200-km bins and the same height range. Do not combine the two hemispheric sides of the channel.
- Compare days 50–75 with 75–100 for descriptive persistence. These correlations are not confidence intervals and do not establish equilibrium.
- Plot full-width and common-interior mass streamfunctions. The native mass-flux definition is ρv = −∂z Ψ and ρw = ∂y Ψ. The displayed maps use fixed 200-km meridional and 5-km vertical Gaussian smoothing; quantitative statistics are unsmoothed. The displayed limits are shared between the widths.

## Momentum accounting

Both runs saved centered zonal velocities and eddy covariances. The analysis uses a common conservative centered-flux decomposition: horizontal and vertical mean momentum flux divergence and horizontal and vertical eddy momentum flux divergence. Centered derivatives are second order inside, first order at the outermost cells, matching the wider-run observer.

Only the wider run saved every-step, correctly weighted SSPRK3 native momentum integrals. Sum their actual increments across segment boundaries, taking differences within a segment when a requested window begins at day 75. Check complete coverage and independently compare the result with endpoint zonal momentum changes. Separate:

1. Exact native time-mean Coriolis, diffusion, external forcing and advection.
2. Three-hour sampled centered mean and eddy advection.
3. Sampled native advection minus the centered decomposition (reconstruction difference).
4. Exact integrated native advection minus sampled native advection (sampling difference).
5. Independent total closure residual and exact-versus-sampled Coriolis difference.

Neither reconstruction nor sampling difference is identified as a physical eddy stress. The narrow run's unresolved residual also includes native diffusion and forcing, which cannot be recovered exactly from its saved centered fields. No replacement simulation is launched to fill that gap.

Report eddy-force alignment with −f v, magnitude, persistence and resolved shear conversion. Regional shear conversion is not regional force work: flux through region edges matters. A closed total budget is necessary but insufficient for a Ferrel-like mechanism; robust interior eddy driving, boundary sensitivity, and resolution remain separate questions.

## Checks and outputs

The script reproduces the previously published narrow-run day 50–100 velocity RMS values, independently reproduces saved wider-run centered terms at every snapshot, checks common coordinates and physical parameters, and verifies every-step native integral coverage and closure. It writes four figures, a numerical report, and a checksum/validation manifest. All processing and figure rendering stay on Unity. Publication reuses the existing synchronized three-depth vorticity/tracer movies.

Moving the walls also changes allowed channel modes and the total f range. A single pair cannot uniquely identify causality or establish statistical convergence. The cap, bottom boundary, dry dynamics and thermal feedback remain unchanged. No residual/Lagrangian circulation, ammonia chemistry, or thermodynamic indirectness claim is made by this audit.
