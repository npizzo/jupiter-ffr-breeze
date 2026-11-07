#!/bin/bash
#SBATCH --job-name=ffr-prepare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/prepare-%j.out
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
"$JULIA_BIN" --project=. scripts/setup.jl
"$JULIA_BIN" --project=. scripts/run.jl configs/rest.toml data/cpu-rest CPU
"$JULIA_BIN" --project=. scripts/check_output.jl data/cpu-rest
"$JULIA_BIN" --project=. scripts/run.jl configs/smoke.toml data/cpu-smoke CPU
"$JULIA_BIN" --project=. scripts/check_output.jl data/cpu-smoke
