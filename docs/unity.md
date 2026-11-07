# Unity workflow

Use a separate project directory and Julia depot. Do not modify environments of existing simulations.

Create `local.env` (ignored by Git) containing absolute paths:

```sh
export FFR_ROOT=/your/project/jupiter-ffr-breeze
export JULIA_BIN=/your/runtime/julia-1.12/bin/julia
```

From the repository directory, create `logs/`, then submit `sbatch slurm/prepare.sbatch`. After successful dependency preparation and CPU validation, submit `sbatch slurm/pilot.sbatch`. Check the current Unity partition/constraint names and your allocation before use. The GPU script requests one GPU, four CPUs and 8 GB host memory for at most 30 minutes. It first runs a short GPU check, then the four-hour physical-time pilot; no retry or larger run is automatic.

Inspect `sacct` MaxRSS, elapsed time, GPU model, and actual usage before allocating the next case. Metadata `wall_seconds` includes simulation initialization, first-step compilation, and saved output, but excludes process startup and model construction. The `post_warmup_*` fields measure the interval after the first nonzero-time output, including output overhead. Use this interval for preliminary throughput estimates; it is not a controlled solver benchmark. Retain job logs outside Git. Transfer only selected raw run directories for analysis and publish the generated figures and metadata.

The runtime download must come from the official Julia distribution and match its published SHA-256. Keep the generated Manifest.toml with the project for reproducibility.

The controlled 48-hour campaign uses `sbatch slurm/convection48.sbatch`. It first validates the added diagnostics, then runs four cases sequentially within one GPU allocation. See [experiment design](convection48.md). Analyze downloaded output with `python scripts/analyze_convection.py data/convection48-v1 --output gallery/convection48`.

For the revised upper-sink and stable-cap experiment, use `sbatch slurm/controlled48.sbatch`. It requests one GPU and 45 minutes for validation plus six sequential 48-hour cases, including one vertically refined run. See [the design](controlled48.md).
