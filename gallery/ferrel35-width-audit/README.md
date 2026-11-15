# Channel-width comparison: 35°N, days 50–100

[Read the analysis and watch the three-depth movies in the GitHub Pages gallery](https://npizzo.github.io/jupiter-ffr-breeze/#ferrel35-width-audit).

All computation used saved output on Unity; no simulation was launched or modified. Job 1788246 completed successfully. The initial short analysis job stopped because the published reference report was absent; the exact archived report was restored and the check passed on the rerun.

The strongest overturning and meridional-flow profiles peak within 0–200 km of each wall in both widths. Wide-channel wall-strip streamfunction RMS is about 10–13 times the common-interior value. Interior streamfunction RMS decreases from 169.1 to 109.4 kg m⁻¹ s⁻¹ while zonal-velocity RMS increases from 1.839 to 3.178 m/s. The interior contains a persistent broad pattern, but its unsmoothed streamfunction cosine between the two 25-day windows is only 0.205 (0.726 after fixed smoothing).

The wider-run native momentum budget closes with a central residual RMS of 1.24e-18 m/s². Its sampled centered decomposition has a reconstruction remainder of 1.22e-5 m/s², versus resolved eddy acceleration 2.39e-6 m/s². Eddy-force alignment with negative Coriolis is encouraging at broad scales (cosine 0.955 smoothed, 0.467 raw), but is not by itself a closed physical attribution. The older narrow run lacks every-step native momentum integrals.

One realization per width; changing width changes allowed modes and total f range. Wall association does not uniquely prove causation. No boundary-independent Jovian Ferrel-cell identification is claimed.

- [Numerical report](summary.json)
- [Validation and checksums](PASS.json)
- [Method and source](../../scripts/ferrel35_width_audit_v1/README.md)
- [Circulation](circulation.png), [interior persistence](interior_persistence.png), [wall profiles](wall_distance.png), [momentum balance](momentum_balance.png)

Figures use labeled 200-km/5-km smoothing; primary norms are unsmoothed and density weighted. Raw model fields remain on Unity. The existing wider and narrow three-depth movies are reused, without rerendering or changing color scales.
