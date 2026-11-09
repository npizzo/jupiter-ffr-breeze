# Completed comparisons and large-eddy validation — November 9, 2026

Watch the synchronized movies on the [gallery landing page](https://npizzo.github.io/jupiter-ffr-breeze/). This entry covers completed ten-day heating pairs, both original rotation cases through twenty days, and the one-day / two-hour validation experiments for the new large-eddy campaign. These are completed segments or tests, not completed forty-day or twenty-day production campaigns.

## Matched heating-on/off comparison

`shear_on_off_10days.mp4` contains all 81 saved times, three depth columns, and four rows: heating-on vorticity, heating-on concentration, heating-off vorticity, heating-off concentration. Both cases use the same physical times, domain, prescribed cyclonic shear, random seed and initial fields. Initial saved vorticity, concentration, mean temperature profile and depth coordinates are checked for exact equality during rendering preparation.

The grid is 1024 × 1024 × 64 over 6000 × 6000 × 120 km. This pair uses the original 300-km random seeds, not the new 1200-km initialization. Bottom heat flux is 162 or 0 W/m². Cap relaxation and the global thermal-budget feedback remain active in both cases.

| Last-day diagnostic | Heating on | Heating off |
| --- | ---: | ---: |
| Volume RMS vertical velocity (m/s) | 1.175789 | 0.005817 |
| Volume RMS horizontal velocity (m/s) | 19.461032 | 19.454469 |
| Mean residual external heat input over ten days (W/m²) | 0.00636569 | 0.000412282 |
| Maximum relative tracer-inventory drift | 5.36e-13 | 3.43e-13 |

Heating sustains vertical motions while the heating-off case decays. Horizontal speed and most kinetic energy include the imposed shear, so they must not be interpreted as convection-generated turbulence. Morphology alone does not establish an inverse cascade or reproduce Jupiter's folded filamentary regions.

## Twenty-day f-plane/beta-plane comparison

`f_beta_20days.mp4` includes all 161 saved times from days 0–20. Each case joins its two completed ten-day segments; the second segment begins with frame 81, so no boundary frame is duplicated or omitted. The top row is f-plane and the bottom row is beta-plane, with three real depth columns. These older runs did not save passive tracer.

Final-day mean kinetic energies differ by 0.0154% (beta relative to f). This bulk diagnostic alone does not measure spectral transfer or jet organization. These older runs also retain their original thermal treatment, so a comparison with the newer shear campaign changes more than rotation or shear alone.

## New 1200-km eddy validation

- `large_eddy_f_validation.mp4`: heating-on/off f-plane tests through one Earth day, 128 × 128 × 32.
- `large_eddy_beta_validation.mp4`: the matching beta-plane tests through one Earth day, 128 × 128 × 32.
- `large_eddy_grid_probe.mp4`: the 1024 × 1024 × 64 beta-plane heated probe through two simulated hours; three actual saved frames shown at one frame per second, without temporal interpolation.

All retain the full 6000-km physical domain. The random flow has wavelengths around 1200 km, Ro = 0.1 and full-depth coherence. Initializer, restart, wall, finite-field, tracer, cap and thermal-budget checks passed before production was released to the queue. The short full-grid probe estimates 6.35 wall hours per ten simulated days; actual production timing may differ. The validation movies are not twenty-day production results.

## Color scales and provenance

Vorticity uses an asinh scale with linear width 0.1 in ζ/f₀. In each comparison, each depth shares the full extrema of both cases over all saved times. Concentration uses a shared fixed asinh scale with width 0.01 and limits 0–1. Colorbar labels are physical values; no saved field values are clipped. No per-frame normalization or temporal interpolation is applied. Standard video encoding is lossy; the figures are for visualization, and numerical analysis uses the saved fields.

`summary.json` records diagnostics computed from the complete histories. Metadata, rendering-export records, sampled diagnostics and per-movie manifests accompany the movies. The join records identify the two original rotation segments. Raw simulation volumes and restart files remain on Unity and are not committed to the repository.

## Re-rendering

Use `scripts/export_video_fields.py` on each completed segment. For a twenty-day case, place its two exported folders under a common parent named `segment_1` and `segment_2`; the renderer's collector reads them in order and requires every expected timestamp. Then run:

```sh
python scripts/compare_flow_tracer_movies.py render/f20 render/beta20 --labels f-plane beta-plane --title '6000-km domain | first 20 Earth days | same color scales' --fps 12 --output f_beta_20days.mp4
python scripts/compare_flow_tracer_movies.py render/shear_heated10 render/shear_control10 --labels 'Heating ON' 'Heating OFF' --title 'Cyclonic shear | matched first 10 Earth days | same color scales' --tracer --fps 8 --output shear_on_off_10days.mp4
```

The validation pairs use the same comparison command with `--tracer --fps 2`; the full-grid probe uses `flow_tracer_video.py --vorticity-linear-width 0.1 --fps 1`. Movie manifests record normalization, frame times and renderer hashes.

## Additional completed f-plane segment: day 30

[Watch the three-depth movie](https://npizzo.github.io/jupiter-ffr-breeze/#f-plane30). Job **66946546** completed day 20–30 successfully (2:17:31 elapsed). `f_plane_30days.mp4` covers **all 241 saved times from day 0 to day 30** at 12 frames/s, with no missing or repeated restart-boundary frames. The three synchronized columns use actual heights 19.6875, 59.0625 and 100.3125 km above the bottom.

Symmetric asinh vorticity scales use width 0.1 in ζ/f₀ and full extrema over all 241 frames, fixed separately for each depth. Limits may differ from the twenty-day f/beta comparison. This older run saved no passive tracer and retains the original cap relaxation, not the later balanced-heating feedback. It continues to gain heat overall. The matching beta-plane day 20–30 segment remains queued; the campaign target remains forty days per case. A longer movie alone does not establish an inverse cascade or jets.

Records: `f30_provenance.json`, `f30_segment3_metadata.toml`, `f30_segment3_export.json`, `f30_segment3_sampled_diagnostics.csv`, `f_plane_30days.json`, and `f30_media_check.json`. Earlier segments retain their original records in this directory.

Reproduction, with rendering exports joined chronologically in `work/gallery-update-f-first30`:

```sh
python scripts/flow_tracer_video.py work/gallery-update-f-first30 --vorticity-only --vorticity-linear-width 0.1 --label 'f-plane | 6000-km domain | completed through day 30' --fps 12 --crf 28 --output gallery/completed-2026-11-09/f_plane_30days.mp4
```
