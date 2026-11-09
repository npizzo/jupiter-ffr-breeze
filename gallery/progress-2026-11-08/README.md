# Production and heat-budget progress — November 8, 2026

**[Watch the synchronized movies](https://npizzo.github.io/jupiter-ffr-breeze/)**

## First ten days of the large f-plane run

The first 10-day segment completed on a 1024 × 1024 × 64 grid in a 6000 × 6000 × 120 km channel. Slurm elapsed time was 2 h 19 min. This is the first quarter of a planned 40-day experiment, not a completed 40-day campaign. The matching beta-plane run is still in progress; the snapshot below records the check time.

[Three-depth movie](large40_f_first10.mp4) shows all 81 saved times, with three depth columns near 20, 60 and 100 km. This earlier campaign did not carry a passive tracer. Its movie therefore shows vorticity only; no concentration has been inferred.

The protected cap core retained at least 99.974% of its initial N². Final domain-wide vertical RMS velocity was 1.225 m/s. This older thermal setup still has net external heating (166.05 W/m² at day 10); it has not been changed while running. No inverse-cascade or FFR-reproduction claim follows from this movie alone.

## New heat-budget and tracer test: forty days complete

[Synchronized flow and tracer movie](shear40_budget_asinh.mp4) has three depth columns, vorticity above concentration, and all 321 saved times. This is the 600-km, 64³ validation box, not the new 6000-km shear production run.

The 40-day mean residual external enthalpy input was **0.00173 W/m²**, compared with 162 W/m² imposed at the bottom. The maximum sampled imbalance just before a controller update was **0.0951 W/m²**. The protected cap retained at least **99.991%** of its initial N². The tracer stayed nonnegative and the source-budget, inventory, passivity and restart checks passed.

The distributed feedback controls imposed enthalpy input; these results do not establish exact total-energy conservation or a realistic Jovian radiative equilibrium. See [the controller design](../../docs/shear20.md).

## Full-grid checks and scheduling

Both 1024² × 64 and 1536² × 64 probes completed on an 80-GB H100, including output and checkpoints. Initial post-warmup throughput projects **6.39 and 15.33 wall hours per ten simulated days**, respectively. Each probe covers only two simulated hours; later convection can change the timestep and cost. Production requests 24 hours per segment.

As of November 8 at 8:56 p.m. Eastern, the first shear production job was eligible and waiting for resources. The remaining heated, heating-off, refined and no-shear segments depend on successful predecessors. The original large-domain beta run had reached 4.125 Earth days. These are timestamped observations, not a live queue display.

[Numerical summary](summary.json) · [Campaign job IDs](../../docs/shear20-launch.json)

Movie data are lossless subsets of the original saved vorticity/tracer planes; no spatial or temporal resampling is applied before rendering. MP4 encoding is lossy. Per-movie JSON records the fixed color scales and saved times; metadata records actual heights and source hashes. Diagnostic CSVs here are sampled at the three-hour movie cadence. Summary extrema and integrated budgets use the full diagnostic histories, retained with the raw fields on Unity.

## Tracer display correction

The revised movie uses an asinh color normalization with linear width 0.01 and fixed concentration limits 0–1, shared across all times and depths. Colorbar labels are actual concentrations. This reveals diluted tracer that was compressed into the purple end of the earlier linear scale. The final full-domain concentration range is 0.01888–0.12975. The fields are unchanged and no values are clipped. Late-time spatial uniformity can still reflect real mixing. The [earlier linear-scale render](shear40_budget.mp4) is retained for comparison.
