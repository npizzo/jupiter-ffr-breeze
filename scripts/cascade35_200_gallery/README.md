# Heated 35°N beta-plane movie, days 0–200

All 1601 three-hourly outputs, at 12 fps, in a synchronized two-row, three-column movie. Vorticity is above passive-tracer concentration, at heights 19.6875, 59.0625 and 100.3125 km above the bottom. Days 0–100 come from the original campaign; days 100–200 come from its checkpoint continuation. No simulation or raw fields are changed or downloaded.

The Unity-only preparer checks completed segment metadata, unchanged physical configuration, diagnostic completion, full chronological coverage, finite three-depth fields, and full movie decoding. Fixed symmetric asinh vorticity scales use the full extrema at each depth over all 200 days. Tracer retains the BrBG scale centered at C=0.0793082, width 0.001, full physical range, fixed across all times and depths. No temporal interpolation, field clipping or per-frame normalization. These scales may differ from earlier shorter movies because the full vorticity extrema are taken over 200 days.

Submit `sbatch slurm/beta35-full200-gallery.sbatch` on Unity. Output: `rendered/beta35-full200`. Existing output directories are refused. All scientific computation and rendering remain on Unity.
