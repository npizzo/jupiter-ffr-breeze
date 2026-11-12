# 55°N beta-plane over heated f-plane: complete days 0–100

[Watch the synchronized comparison on the GitHub gallery](https://npizzo.github.io/jupiter-ffr-breeze/#beta55-f100).

Four rows by three depths, in this order: heated beta-plane vorticity, beta-plane passive tracer, heated f-plane vorticity, f-plane passive tracer. Actual heights are 19.6875, 59.0625 and 100.3125 km above the bottom. All 801 three-hourly saved times are included, with no time interpolation, for a 66.75-second movie at 12 fps.

Both dry runs use the same initial eddy pattern and speed, heating of 162 W/m², 6,000 × 6,000 × 120 km domain and 256 × 256 × 64 grid. The initial saved vorticity and tracer slices were verified identical. Segment parameters match except campaign name, stop time and rotation/beta configuration. The beta-plane has beta=2.8216700143325537e-12 m⁻¹ s⁻¹; the reference has beta=0. Both have f₀=2.880957739764384e-4 s⁻¹ (see original metadata for the exact stored value).

Vorticity uses symmetric asinh width 0.1 in ζ/f₀, with one full-extrema range per height shared by both cases over the entire movie. Positive vorticity is cyclonic in the northern hemisphere. Tracer uses the requested improved BrBG contrast mapping centered at C=0.0793082, two-sided asinh width 0.001, and the full physical concentration range, fixed across both cases, all heights and all times. The fields are unchanged; no framewise normalization, subtraction or clipping. Exact limits and scales are in the movie JSON. The poster shows day 100.

The final 55°N segment (job 1447348) completed successfully, including its stable-cap, budget, checkpoint, tracer and three-depth checks. Rendering job 1498827 ran entirely on Unity and verified all 801 encoded frames by full decoding. The original queued rendering job 1447349 was preserved. This separate variant implements the requested row order and improved tracer colors.

PASS.json records sizes, SHA-256 and Git blob hashes for the finished publication assets. The original segment metadata and rendering provenance are included; the raw model fields stay on Unity. QA stills show days 0, 50 and 100. Rendering code: scripts/beta55_f100_contrast; submission: slurm/beta55-f100-contrast.sbatch.

These coarse dry-model movies do not establish inverse energy transfer, jet equilibration, resolved filament widths or correspondence between passive dye and observed clouds. Interpretation requires the saved diagnostics and resolution checks.
