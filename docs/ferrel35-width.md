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

The main comparison should use the same central physical region in both channels (|y| ≤ 2,400 km, heights 10–75 km), alongside separately labeled whole-channel and distance-from-wall diagnostics. Doubling width increases the total range of f(y), but f(y) at each shared central coordinate is unchanged. Using the same fractional width for the primary diagnostic would compare different physical regions and should be avoided.

## Momentum diagnostics

The isolated driver is copied from the existing 35°N code, adding observers. Native advection, Coriolis, diffusion and external zonal forcing are evaluated at all three SSPRK3 stages, checked against the native RHS and integrated with the correct stage weights on every complete timestep. Independent momentum changes are compared with those integrated terms. Pressure-gradient cancellation in the periodic zonal direction is also checked. Detailed stage files are saved sparsely; cumulative arrays are saved with each three-hour frame. This improves on the earlier Ferrel pilot's sampled-step closure test.

Centered mean/eddy flux estimates, heat and tracer transport, mean kinetic-energy conversion, and the native-minus-centered advection remainder are retained. The remainder must not be labeled automatically as physical eddy transport. The day200–210 review separately reports exact time-integrated native terms, the three-hour sampled mean/eddy estimates, and the sampling discrepancy on their common interval. A restored callback may omit the duplicate checkpoint-time frame; the review uses the actual saved interval rather than assuming it begins at day200.

A closed total momentum budget is a numerical prerequisite, not evidence that eddies drive thermally indirect circulation. Boundaries and cap treatment remain possible causes; the wider run retains the same zonal extent and changes allowed meridional modes as well as wall distance. Further seeds and resolution tests would be needed for general conclusions.

## Validation and release

1. Unity CPU continuous-versus-restarted test on a 64×128×32 wider grid, including cumulative momentum integrals and ordinary conservation, cap and three-depth output checks.
2. Unity GPU repeats the restart tests, then compares an instrumented restart from the actual day200 checkpoint against the original driver over 1,200 seconds. Physical fields must agree within the existing restart tolerances.
3. A full 256×512×64 six-hour GPU probe checks output, cap, conservation and runtime. Estimated ten-day production time must be under eight hours for the twelve-hour allocation.
4. The full day200–210 diagnostic branch must complete with every-step and cumulative closure checks and a numerical review before wider production starts.

The v4 observer explicitly initializes its first stage after a checkpoint restore. The first CPU test caught a missing initial restart timestep in the diagnostic accumulator (the flow, tracer and cap checks passed). That validation failure blocked the original dependent chain before any new production ran. The repair adds checks that the budget starts at the exact restart time, ends at the exact segment endpoint, and covers every intervening timestep without gaps. CPU and GPU restart tests are repeated for the repaired version. No existing physical simulation was canceled.

A failed dependency cancels descendants rather than allowing unvalidated production. Each wide segment verifies source, configuration and runtime hashes and requires the previous segment to complete. GPU jobs request one GPU with at least 48 GB VRAM in the URI GPU partition; the actual device is recorded. The existing verified CUDA 12.2 environment is reused without modification.

## Queue and output

[Job identifiers and dependencies](../campaigns/ferrel35-width-v4/jobs.tsv) and [publication-time status](../campaigns/ferrel35-width-v4/status.json) are snapshots, not a live scheduler.

- Source: `scripts/ferrel35_width_v4/`; configuration: `configs/ferrel35-width-v4/`.
- Diagnostic branch: `data/ferrel35-width-v4/diagnostic/segment_1` on Unity.
- Wide run: `data/ferrel35-width-v4/wide/segment_1` through `segment_20` on Unity.
- Movie jobs follow the diagnostic branch and wide day200, respectively. They combine vorticity above tracer, with three simultaneous depth columns, fixed physical scales and explicit timestamps. The wide movie preserves the channel's 2:1 north–south/east–west aspect ratio. The day200 diagnostic starting frame comes from the actual original day200 output.
- The same improved BrBG tracer contrast and per-depth fixed asinh vorticity normalization are retained. Finished movies must pass full decoding and frame-count checks before gallery publication.

[Current gallery](https://npizzo.github.io/jupiter-ffr-breeze/#beta35-200). No new production results are claimed by this launch record.

## Launch verification

CPU job 1693660 and GPU preflight job 1694154 passed. The GPU test ran on an NVIDIA L40S with CUDA 12.2; all 21 restart field/integral comparisons and all 20 evolved-state comparisons passed, along with the six-hour full-size probe. The review script was also exercised successfully on both the evolved-state and wider-grid preflight outputs. Measured post-startup cost projects approximately 3.24 GPU-hours per ten simulated days in the wider channel; evolving velocities and queue delays can increase the total elapsed time.

The diagnostic job is 1694155; wider segments are 1694157–1694176, followed by movie job 1694177. The diagnostic movie job is 1694156. The earlier GPU attempt pointed to a stale preflight directory and stopped before a model run; the replacement job uses the correct v4 gate. Current and superseded submission records are retained.
