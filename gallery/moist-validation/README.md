# Moist implementation validation

**[Watch all movies on the GitHub Pages gallery](https://npizzo.github.io/jupiter-ffr-breeze/#moist-validation)**

One hour in an 80 × 80 × 40 km rotating hydrogen–helium box. Each movie contains three depth columns (near 4, 12 and 24 km), with vorticity on top and passive concentration below. Moist movies add liquid-plus-ice concentration in a third row. Heights are above the bottom. All 31 saved frames are retained. These are plume-validation experiments, not FFR simulations or cloud-radiance predictions.

- [Dry 64³ movie](dry_plume_asinh.mp4)
- [Moist 64³ movie](moist_plume_asinh.mp4)
- [Moist 96³ movie](moist_plume_fine_asinh.mp4)

The revised transport scheme passes water/tracer conservation, positivity, phase closure, rest and condensate-formation checks. Peak vertical RMS changes by 0.48% between moist grids, but maximum sampled vertical speed changes by about 48%. Local extremes are not converged.

[Method and limitations](../../docs/moist-validation.md) · [Numerical summary](summary.json) · [Validation record](PASS.toml). Per-case metadata, sampled histories and movie color-scale/time records are included here.

Tracer panels now use a fixed nonlinear asinh color scale (linear width 0.01, limits 0–1), with actual concentration tick labels. It is shared across all depths and times. The original fields and diagnostic results are unchanged.
