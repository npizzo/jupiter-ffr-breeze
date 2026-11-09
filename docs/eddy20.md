# Larger random eddies: matched heating and rotation experiments

Four new 20-Earth-day experiments use the same 6000 × 6000 × 120 km channel and 1024 × 1024 × 64 grid as the shear campaign. Each runs in two ten-day segments. These are new initial-value experiments, not restarts of the existing 300-km seed cases.

| Case | Rotation | Bottom heat flux |
| --- | --- | --- |
| eddy20_f_heated | f-plane | 162 W/m² |
| eddy20_f_control | f-plane | 0 W/m² |
| eddy20_beta_heated | beta-plane | 162 W/m² |
| eddy20_beta_control | beta-plane | 0 W/m² |

The central Coriolis parameter is 2.880957739764384e-4 s⁻¹ (55°N). Beta is zero or 2.8216700143325537e-12 m⁻¹ s⁻¹. All four retain the same imposed initial cyclonic shear, with peak ζ/f₀ = 0.05 and maximum speed 27.51 m/s. It is not maintained by a shear force. Cap relaxation, heat-budget feedback, scalar transport, diffusivities, tracer patches and thermal perturbations are unchanged from the shear campaign. Heating-off means no bottom heat input; cap relaxation and heat-budget feedback remain enabled.

## Initial random flow

Random Fourier/sine modes lie in the annulus 0.75 k₀ ≤ |k| ≤ 1.25 k₀, where k₀ = 2π/(1200 km). Thus the permitted wavelengths are 960–1600 km, subject to discrete channel modes. There is no random zonal-mean mode. Seed 31415 is used in every case.

A corner streamfunction generates staggered velocities u = ψ_y and v = −ψ_x. The random horizontal flow is identical at every height, including the stable cap. Initial w is zero. This gives a discrete mass-flux divergence of roundoff size for a reference density depending only on height, and exactly zero normal flow at the meridional walls and top/bottom.

The face-weighted random horizontal RMS speed is normalized to 5.5022 m/s, giving Ro = U/(f₀ L) = 0.1 with L = 1/k₀ = 191 km. This normalization excludes the separately imposed cyclonic shear. Compared with the earlier 300-km seeds at the same Rossby number, random speed is four times larger and speed-squared is sixteen times larger. Vertical structure also changes. Consequently the old/new comparison does not isolate seed scale; a future equal-speed comparison can do that.

The random kinematic-pressure estimate is −f₀ ψ interpolated to cell centers, added to the existing shear pressure. It matches the linear f₀ Coriolis balance under centered staggered interpolation. The model pressure remains diagnostic and evolves through its pressure solve. Nonlinear acceleration and beta-plane adjustment are retained; this is not an exact steady equilibrium. The maximum fractional variation of f from the channel center is about 2.94% for beta-plane cases.

## Validation and queue gates

1. CPU initializer checks: actual prescribed RMS/Rossby number, divergence, wall velocities, full-depth coherence, linear f₀ pressure consistency, reproducibility, independent random seed, and identical fields across the four configurations. Repeat constraint checks at production horizontal resolution with two vertical levels.
2. L40S validation: a checkpoint/restart comparison and four 24-hour simulations at 128 × 128 × 32 over the full physical domain. Check finite fields, tracer bounds and conservation, cap stability, thermal source accounting, initial-field matching, checkpoints and three-depth output.
3. H100 production-grid probe: 1024 × 1024 × 64, two simulated hours. Require all output checks and an extrapolated ten-day runtime below 18 hours before allowing 24-hour production allocations.
4. Eight H100 production segments, chained with successful-completion dependencies. The order is f heated, f control, beta heated, beta control through day ten, then the same order through day twenty. Any failed gate prevents dependent production jobs from starting.

The implementation is isolated in `scripts/eddy1200/`, with new configs and submission scripts. Existing campaigns are not edited or cancelled. Runtime driver, initializer and heat-controller hashes are recorded, and continuations reject changed source/configuration.

## Output and interpretation

Save every three hours at the actual planes nearest 20, 60 and 100 km above the bottom. Gallery movies will show all three depths simultaneously, vorticity over tracer, with fixed and labeled color scales. Depth-mean velocities and Reynolds stresses are retained as in the existing campaign.

This comparison asks whether convection deforms, fragments, stretches or sustains pre-existing larger eddies and their tracer patterns, and how beta changes that response. Seeded large structures by themselves do not establish an inverse cascade or Jupiter folded-filament reproduction. The model remains a dry, idealized anelastic atmosphere.

Submission IDs and gate status are recorded in [eddy20-launch.json](eddy20-launch.json).
