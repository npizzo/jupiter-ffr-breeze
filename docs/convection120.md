# Five-day 30× convection and resolution experiment

All four simulations completed as Unity array **66932767**, one L40S at a time. Each starts from the seeded initial state and runs for **120 Earth hours**; these are fresh integrations, not restarts of the 48-hour output. Output is isolated in `data/convection120-v1`.

| Array task | Case | Bottom heat input [W/m²] | Upper heat input [W/m²] | Grid |
|---|---|---:|---:|---|
| 1 | heated, z64 | 162 | −162 | 256 × 256 × 64 |
| 2 | cooling only, z64 | 0 | −162 | 256 × 256 × 64 |
| 3 | heated, z128 | 162 | −162 | 256 × 256 × 128 |
| 4 | cooling only, z128 | 0 | −162 | 256 × 256 × 128 |

The 30× heating level produced a large response in the completed campaign with lower vorticity relative to f than 100×. This experiment asks whether that response persists and organizes over five days, and how sensitive it is to vertical resolution. Compared with the original 30× run, horizontal spacing halves from 9.375 to 4.6875 km. The two new pairs halve vertical spacing from 1.875 to 0.9375 km. Compare the first 48 hours against the original case for horizontal sensitivity, and the full new integrations for vertical sensitivity and longer evolution. This is not a full three-grid convergence study.

The model, 1,200 × 1,200 × 120-km box, initially neutral lower layer, stable cap above 80 km, viscosity, tapered diffusivity, seeds and initial eddy Rossby 0.1 stay fixed. Jovian f remains 2.88096 × 10⁻⁴ s⁻¹, with beta = 0 and no imposed shear. Across resolutions, the initial velocity modes and phases are shared but sampled and normalized on each discrete grid; the tiny thermal noise realization changes with the grid. Heated/control configurations at each resolution differ only in bottom heat input and case name. The 162 W/m² is an idealized enhanced-forcing choice, not an inferred Jovian flux.

## Execution and validation

`sbatch slurm/convection120.sbatch` submits an array with concurrency one, four CPUs, 24 GB host memory and a three-hour limit per case. This bounds simultaneous use to one GPU and preserves independent case outputs if another case fails. The solver remains pinned to the validated source and dependency versions. This run does not implement restart checkpoints; use a new output root for any rerun.

Every case runs the existing initial mass-continuity, amplitude, reproducibility, finite-output, source and thermal-inventory checks. Completion status and checks must be confirmed before interpreting results. Save full five-minute diagnostic histories and hourly field snapshots. Report means separately for hours 24–48 and 96–120; do not silently reuse the previous analyzer's default hours 24–end as a last-day average.

Compare column kinetic energy, vertical speed, vorticity/f, buoyancy work, heat transport, stable-cap evolution and spatial spectra against the matching cooling-only control. Assess whether bulk histories settle and whether small-scale vorticity remains resolution-sensitive. Longer-lived patterns alone do not establish FFR morphology or an inverse cascade.

## Vorticity videos for this and future galleries

The gallery output standard now includes **MP4 vorticity-evolution videos**, alongside still snapshots and metrics. Use every hourly saved frame, annotate actual simulated time, and keep one symmetric ζ/f color scale fixed across all frames and both members of a pair. Do not interpolate frames or rescale each frame. Each video includes a final-frame poster and a JSON record of source cases, times, frame rate and color limits. Map height is labeled at its actual cell-center position (near 20 km); this differs slightly with vertical resolution.

Install system FFmpeg with the libx264 encoder, in addition to the Python requirements. After retrieving every frame from completed runs:

```sh
python scripts/vorticity_video.py \
  data/convection120-v1/convection120_f030_heated_h256_z64 \
  data/convection120-v1/convection120_f030_cooled_h256_z64 \
  --output gallery/convection120/z64_vorticity.mp4

python scripts/vorticity_video.py \
  data/convection120-v1/convection120_f030_heated_h256_z128 \
  data/convection120-v1/convection120_f030_cooled_h256_z128 \
  --output gallery/convection120/z128_vorticity.mp4
```

At six frames per second, 121 hourly snapshots yield about 20 seconds per video. The renderer rejects missing frames and mismatched comparison grids/times by default. `--allow-sparse` is for explicitly labeled local previews, not completed gallery movies. Include clickable MP4 links and poster images in the gallery so the files remain usable even where Markdown does not play video inline. Render and inspect the movies before publication; no completed five-day results or movies are claimed at submission.

Completed results: [interpretation](convection120-results.md) and [videos and snapshots](../gallery/convection120/README.md). Public MP4s are compact 800-pixel H.264 viewing copies (CRF 32) of the 1200-pixel renders; no frames are dropped. Full renders remain on Unity.
