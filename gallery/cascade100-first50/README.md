# Large-eddy f-plane comparison: days 0–50

[Watch the three movies on the gallery landing page](https://npizzo.github.io/jupiter-ffr-breeze/#cascade50).

Both matched cases completed fifty Earth days. This is the first half of the planned hundred-day campaign, not a completed hundred-day result.

| Movie | Layout | Saved frames |
|---|---|---|
| [Heating on/off comparison](comparison_0_50days.mp4) | Four rows: heated vorticity/tracer, then unheated vorticity/tracer | 401 |
| [Heating on](heated_0_50days.mp4) | Vorticity above tracer | 401 |
| [Heating off](unheated_0_50days.mp4) | Vorticity above tracer | 401 |

Every movie has three synchronized columns at **19.6875, 59.0625 and 100.3125 km above the bottom**. Real saved times span days 0–50 at three-hour intervals. Playback is 12 frames/s (about 33.42 seconds). Restart boundaries are joined without duplicate frames or temporal interpolation. Posters show day 50.

## Matched configuration

- Dry anelastic hydrogen–helium atmosphere, Jovian f-plane at 55°N: f = 2.88095774 × 10⁻⁴ s⁻¹, beta = 0.
- Domain 6,000 × 6,000 × 120 km; periodic zonal direction and free-slip meridional walls.
- Grid 256 × 256 × 64; horizontal spacing 23.4375 km and vertical spacing 1.875 km.
- Identical full-depth random horizontal eddies, configured wavelength band 960–1600 km centered on 1200 km; initial RMS speed 5.5022 m/s and Rossby number 0.1. No imposed shear or jet, and no sustained random momentum forcing.
- Bottom heat flux 162 versus 0 W/m². Both retain the stable cap and enthalpy-feedback controller.
- Passive tracer is dimensionless concentration, not cloud condensate. Neither run includes moisture.

See the [campaign design and diagnostic limitations](../../docs/cascade100.md).

## Fixed display scales

Vorticity is ζ/f₀ on a symmetric asinh scale with linear width 0.1. Per-depth limits are ±1.47202134, ±1.93124425 and ±0.36857685, respectively. These include extrema across both cases and all fifty days and are identical in all three movies. They can differ from earlier ten-day movies because the time window is longer. The unheated flow is therefore paler where its amplitude is smaller; it is not rescaled independently.

Tracer uses a fixed asinh scale with linear width 0.01, range 0–1 and physical concentration ticks, shared across cases, depths and times. Fields are unchanged; there is no automatic per-frame contrast adjustment. Display sampling is nearest-neighbor.

## Verification and provenance

Rendering and media validation ran on Unity in job **67036876**, dependent on successful unheated day-50 job **66992558**. Heated day-50 job **66992557** also completed. All three movies passed frame-count and duration checks on Unity; hashes and file sizes are in [provenance.json](provenance.json).

Movie manifests record every frame time, depth, normalization and renderer hash:
[comparison](comparison_0_50days.json), [heated](heated_0_50days.json), [unheated](unheated_0_50days.json).

Segment metadata:
| Segment | Heating on | Heating off |
|---|---|---|
| 1 | [metadata](heated_segment1_metadata.toml) | [metadata](control_segment1_metadata.toml) |
| 2 | [metadata](heated_segment2_metadata.toml) | [metadata](control_segment2_metadata.toml) |
| 3 | [metadata](heated_segment3_metadata.toml) | [metadata](control_segment3_metadata.toml) |
| 4 | [metadata](heated_segment4_metadata.toml) | [metadata](control_segment4_metadata.toml) |
| 5 | [metadata](heated_segment5_metadata.toml) | [metadata](control_segment5_metadata.toml) |

Published metadata replaces the Unity absolute repository prefix in checkpoint paths with repository-relative paths. Scientific parameters and recorded source/configuration hashes are unchanged. Raw fields and checkpoints remain on Unity; only finished publication assets were transferred for GitHub upload.

To reproduce on Unity after both fifth segments and their diagnostics complete:

```sh
sbatch slurm/cascade50-gallery.sbatch
```

The renderer is [prepare_cascade50_gallery.py](../../scripts/prepare_cascade50_gallery.py), using [render_cascade50_movies.py](../../scripts/render_cascade50_movies.py).

## Interpretation limits

These coarse runs screen large-scale evolution. Fine convective structures and filament thickness are not demonstrated to be resolution-converged. Changes in vortex size or morphology alone do not establish an inverse cascade. Use the scale-transfer diagnostics described in the campaign design before drawing that conclusion. The heating-off case is decaying turbulence with the same cap/controller and dissipation, not a statistically stationary forced reference.
