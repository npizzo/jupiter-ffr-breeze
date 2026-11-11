#!/bin/bash
#SBATCH --job-name=ffr-seedv1-unit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/ferrel2026-seedtest-v1-unit-%j.out
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?This validation requires a Unity allocation}"
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-ferrel2026-gpu:$FFR_ROOT/depot"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
test "$#" -eq 2
"$JULIA_BIN" --startup-file=no --project=environments/ferrel2026_gpu scripts/ferrel2026_seedtest_v1/launch.jl "$2" CPU scripts/ferrel2026_seedtest_v1/validate_initialization.jl "$@"
