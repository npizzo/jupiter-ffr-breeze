# Heated beta-plane latitude comparison, days 0–100

Render the completed 35°N run above the completed 55°N run. Four rows: 35°N vorticity, 35°N tracer, 55°N vorticity, 55°N tracer. Three columns: 19.6875, 59.0625 and 100.3125 km above the bottom. Every saved three-hourly field is included (801 frames at 12 fps).

Both runs use a 6000 × 6000 × 120 km channel, 256 × 256 × 64 grid, 162 W/m² bottom heating, stable-cap relaxation, thermal-budget feedback and free-slip impermeable walls. Initial dimensional eddy velocities, temperature and tracer match to rtol=1e-12, atol=1e-13. Latitude changes both f and beta; equal dimensional speed therefore gives different initial Rossby numbers (approximately 0.143 at 35°N and 0.1 at 55°N). This is not an isolated beta-only sensitivity experiment.

Vorticity is shown as ζ/f₀ using each latitude's central Coriolis parameter. One symmetric asinh color scale per depth is shared by both cases and fixed across all times. Tracer retains the improved BrBG, two-sided asinh mapping centered at C=0.0793082 with width 0.001; limits cover the full physical range and are shared across cases, depths and times. No framewise normalization, subtraction or clipping.

All field checks, rendering, encoding and decoding occur on Unity. `prepare.py` checks completed segment metadata, initial saved fields, transfer-diagnostic completion, synchronized times, layout, movie frame count and full decode before writing PASS.json. Original simulation code and original movies are unchanged.

Submit from the Unity repository root with `sbatch slurm/latitude100-contrast.sbatch`. Output: `rendered/latitude100-contrast`. The render script refuses an existing output directory.
