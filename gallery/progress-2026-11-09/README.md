# First ten-day beta and heated-shear results — November 9, 2026

**[Watch the gallery](https://npizzo.github.io/jupiter-ffr-breeze/)**

The original 6000-km f-plane and beta-plane experiments have both completed their first ten-day segments. The new heated cyclonic-shear/tracer experiment has also completed its first ten days. All are 1024 × 1024 × 64 grids in 6000 × 6000 × 120 km channels; these are completed segments of longer campaigns, not completed 40-day or 20-day campaigns.

## Matched f-plane/beta-plane movie

[Synchronized comparison](f_beta_first10.mp4): f-plane on the top row, beta-plane below; three columns at 19.69, 59.06 and 100.31 km above the bottom. Every saved time is included. Each depth uses one fixed vorticity scale shared by both cases throughout the movie. These earlier runs have no passive tracer.

The final-day mean kinetic energy differs by approximately 0.17% between f and beta. Their final-day vertical RMS speeds are 1.22585 and 1.22555 m/s. This is a small bulk difference over ten days; it is not an ensemble estimate or evidence that beta cannot organize the flow over longer times. Both protected cap cores retained more than 99.97% of their initial N². The older campaigns still receive about 165.97 W/m² net external heating and were not modified while running.

The Slurm jobs took 2 h 19 min each. Their next segments continue toward forty simulated days per case.

## Heated cyclonic shear with tracer

[Synchronized shear/tracer movie](shear_heated_first10.mp4): vorticity above passive concentration, with the same three depth columns and all 81 saved times. Concentration uses the labeled asinh color scale, fixed across time and depths, with actual concentration tick labels and no data clipping.

At day ten the tracer shows elongated, wavy ribbons, particularly in the stable cap. This is a visual description, not a quantitative match to observed Jovian folded filamentary regions. The run completed in 7 h 18 min of Slurm time. Its protected cap retained at least 99.992% of initial N². The ten-day mean residual external heat input was 0.00637 W/m² with 162 W/m² entering at the bottom; the largest sampled imbalance immediately before a controller update was 0.0935 W/m². Tracer mass drift was less than 5.4e-13 relative, and concentration remained nonnegative.

Final-day vertical RMS was 1.176 m/s. The much larger horizontal RMS (19.46 m/s) includes the deliberately imposed initial shear; it must not be interpreted as convection-generated turbulence. The heating-off control and the no-shear case are not yet available. Likewise, comparisons with the older f/beta runs confound shear with the new thermal feedback and scalar-transport treatment.

These numerical checks do not establish a realistic Jovian radiative balance, an inverse cascade, or reproduction of FFR morphology.

## Queue snapshot

At approximately 5:44 a.m. Eastern on November 9, none of the campaign jobs had failed. Both active campaigns were waiting for GPUs. Unity estimated 12:44 p.m. November 9 for the next large-domain segment, and 12:50 a.m. November 10 for the heating-off shear control. These scheduling estimates can move. The refined 1536² production case remains dependent on completion of the control; only its full-grid probe has run.

[Numerical summary](summary.json) · [Large-domain jobs](../../docs/large40-launch.json) · [Shear jobs](../../docs/shear20-launch.json)

The movie inputs are lossless subsets of saved model planes with no spatial or temporal resampling before rendering. MP4 encoding is lossy. Metadata and movie manifests record actual heights, fixed scales, source hashes and times. Public diagnostic CSVs are sampled at the three-hour movie cadence; summary statistics use the full diagnostic histories retained on Unity.
