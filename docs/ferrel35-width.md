# 35°N Ferrel-cell diagnostics and channel-width sensitivity

Authorized November 13, 2026. Computation, diagnostics and rendering run on Unity. Existing production experiments and their source files are unchanged.

This experiment tests whether the persistent overturning in the 6,000-km channel is tied to its sidewalls, while improving the momentum accounting needed to assess eddy driving. It does not assume the observed circulation is a Ferrel cell.

## Two experiments

| Parameter | Diagnostic restart | Wider channel |
|---|---|---|
| Initial state | Actual 35°N heated checkpoint at Earth day 200 | Fresh random flow with matched physical statistics |
| Time interval | Days 200–210 | Days 0–200 |
| Domain, east × north × height | 6,000 × 6,000 × 120 km | 6,000 × 12,000 × 120 km |
| Grid | 256 × 256 × 64 | 256 × 512 × 64 |
| Horizontal spacing | 23.4375 km | 23.4375 km |
| Bottom heating | 162 W/m² | 162 W/m² |
| Rotation | Jupiter beta plane centered at 35°N | Same central f0 and beta |
| Initial eddy wavelength / RMS speed | Inherited evolved state | Approximately 1,200 km / 5.50 m/s |

The cap, thermal feedback, vertical sponge, viscosity, diffusivities, dry anelastic equations and free-slip impermeable boundaries remain unchanged. The larger channel has twice as many meridional cells, not coarser cells. Its random seed number is retained, but the realized channel modes differ; it is a matched-statistics sensitivity, not an identical-field comparison. No checkpoint is stretched.

## Momentum diagnostics

The isolated driver is copied from the existing 35°N code, adding observers. Native advection, Coriolis, diffusion and external zonal forcing are evaluated at all three SSPRK3 stages, checked against the native RHS and integrated with the correct stage weights on every complete timestep. Independent momentum changes are compared with those integrated terms. Pressure-gradient cancellation in the periodic zonal direction is also checked. Detailed stage files are saved sparsely; cumulative arrays are saved with each three-hour frame. This improves on the earlier Ferrel pilot's sampled-step closure test.

Centered mean/eddy flux estimates, heat and tracer transport, mean kinetic-energy conversion, and the native-minus-centered advection remainder are retained. The remainder must not be labeled automatically as physical eddy transport. The day200–210 review separately reports exact time-integrated native terms, the three-hour sampled mean/eddy estimates, and the sampling discrepancy on their common interval. A restored callback may omit the duplicate checkpoint-time frame; the review uses the actual saved interval rather than assuming it begins at day200.

A closed total momentum budget is a numerical prerequisite, not evidence that eddies drive thermally indirect circulation. Boundaries and cap treatment remain possible causes; the wider run retains the same zonal extent and changes allowed meridional modes as well as wall distance. Further seeds and resolution tests would be needed for general conclusions.

## Validation and release

1. Unity CPU continuous-versus-restarted test on a 64×128×32 wider grid, including cumulative momentum integrals and ordinary conservation, cap and three-depth output checks.
2. Unity GPU repeats the restart tests, then compares an instrumented restart from the actual day200 checkpoint against the original driver over 1,200 seconds. Physical fields must agree within the existing restart tolerances.
3. A full 256×512×64 six-hour GPU probe checks output, cap, conservation and runtime. Estimated ten-day production time must be under eight hours for the twelve-hour allocation.
4. The full day200–210 diagnostic branch must complete with every-step and cumulative closure checks and a numerical review before wider production starts.

The final v3 driver reuses the CPU test from v2 only after verifying that all ten tested numerical/physical files and both configurations are byte-identical, and that the runtime environment is unchanged. The v3 changes concern campaign paths and correct treatment of the snapshot sampling interval in the later review. Unused earlier GPU validation submissions were retired before running; no existing physical simulation was canceled.

A failed dependency cancels descendants rather than allowing unvalidated production. Each wide segment verifies source, configuration and runtime hashes and requires the previous segment to complete. GPU jobs request one GPU with at least 48 GB VRAM in the URI GPU partition; the actual device is recorded. The existing verified CUDA 12.2 environment is reused without modification.

## Queue and output

[Job identifiers and dependencies](../campaigns/ferrel35-width-v3/jobs.tsv) and [publication-time status](../campaigns/ferrel35-width-v3/status.json) are snapshots, not a live scheduler.

- Source: `scripts/ferrel35_width_v3/`; configuration: `configs/ferrel35-width-v3/`.
- Diagnostic branch: `data/ferrel35-width-v3/diagnostic/segment_1` on Unity.
- Wide run: `data/ferrel35-width-v3/wide/segment_1` through `segment_20` on Unity.
- Movie jobs follow the diagnostic branch and wide day200, respectively. They combine vorticity above tracer, with three simultaneous depth columns, fixed physical scales and explicit timestamps. The wide movie preserves the channel's 2:1 north–south/east–west aspect ratio. The day200 diagnostic starting frame comes from the actual original day200 output.
- The same improved BrBG tracer contrast and per-depth fixed asinh vorticity normalization are retained. Finished movies must pass full decoding and frame-count checks before gallery publication.

[Current gallery](https://npizzo.github.io/jupiter-ffr-breeze/#beta35-200). No new production results are claimed by this launch record.
