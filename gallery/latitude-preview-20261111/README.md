# Heated beta-plane previews

Watch the movies on the [GitHub Pages gallery](https://npizzo.github.io/jupiter-ffr-breeze/#latitude-previews).

These are partial campaigns: 55°N days 0–40 (321 frames) and 35°N days 0–20 (161 frames), at three-hourly intervals. Both movies show three real saved heights, with vorticity above passive tracer. Rendered and fully decoded on Unity; only finished publication assets were transferred.

The tracer mapping matches the revised full 100-day f-plane contrast movies: fixed BrBG, center 0.0793082, two-sided asinh width 0.001, full physical concentration range. Vorticity uses asinh width 0.1 in ζ/f₀, full per-depth extrema fixed over each movie. Vorticity ranges are not shared across latitude movies. See each movie JSON for exact limits, saved times and heights. No field values were changed.

Each latitude directory contains the movie, poster, provenance, frozen source inventory, original segment metadata, three QA stills and PASS.json with checksums and full-decode results. The capture excludes active files younger than 30 seconds and verifies stable file size/time and finite saved slices. The inventories are frozen before rendering; simulations continue independently.

Rendering jobs: 1470461 (55°N) and 1470462 (35°N, after successful day-20 job 1459468). Scientific model files and running chains were not modified. Renderer: scripts/latitude_preview_20261111; submission: slurm/latitude-preview-20261111.sbatch.

Dry idealized convection, 162 W/m² bottom heating, 6,000 × 6,000 × 120 km domain, 256 × 256 × 64 cells, and approximately 1,200-km initial random eddies. Morphology alone does not establish inverse energy transfer or resolution convergence.
