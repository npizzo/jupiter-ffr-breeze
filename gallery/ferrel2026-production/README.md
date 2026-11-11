# Ferrel transient mechanism screening

Completed production publication assets; generated and decoded on Unity.

[Watch on the GitHub gallery](https://npizzo.github.io/jupiter-ffr-breeze/).

| Case | Grid | Duration (Earth h) | Native sampled closure steps |
|---|---|---:|---:|
| Unheated | 512 × 512 × 96 | 48 | 32 |
| Heated · 162 W/m² | 512 × 512 × 96 | 48 | 43 |
| Heated · boundary sensitivity | 512 × 512 × 96 | 48 | 43 |
| Axisymmetric · unheated | 512 × 512 × 96 | 6 | 11 |
| Axisymmetric · heated | 512 × 512 × 96 | 6 | 11 |

These outputs do not yet establish statistically steady Ferrel cells. The first 12 hours are conservatively excluded as startup, without claiming a measured equilibration time; six-hour controls/refinement probe this early window only.

Movies use actual saved vorticity/tracer at three heights. Time-independent asinh color scales are shared across cases, with a common tracer range and separate full-range vorticity limits for each depth. No automatic per-frame rescaling, field clipping or inferred depths.

See `report.json` for numeric diagnostics, exact averaging windows, limitations, input hashes and media checksums. Native step closure and approximate centered covariance interpretation are reported separately.

- Dry anelastic beta-channel transient mechanism screening; no steady or deep Jovian-cell claim.
- The first 12 Earth hours are excluded as a conservative startup window, not a demonstrated equilibration timescale. Six-hour pilots and controls lie entirely within that window.
- Native momentum closure samples complete SSPRK3 steps; it is not a continuous long-run integrated momentum budget.
- Centered covariance stresses, force work and mean-flow transfer are approximate diagnostics, not exact WENO fluxes or a closed kinetic-energy budget.
- Positive covariance times mean shear feeds zonal-mean kinetic energy. This is not by itself proof of an inverse cascade.
- A single seed, short duration, cap relaxation, sponge and global thermal feedback limit physical attribution.
- The boundary case changes both cap relaxation and vertical-velocity damping; it cannot isolate their individual effects.
- No residual/TEM circulation is inferred in the nearly neutral layer; mean and eddy tracer fluxes are reported directly.
- Psi is overturning transport per unit zonal length in kg/(m s), not a planet-integrated kg/s streamfunction.
- Temporal variability and coherence are descriptive; short windows and a single seed do not establish statistical significance.
