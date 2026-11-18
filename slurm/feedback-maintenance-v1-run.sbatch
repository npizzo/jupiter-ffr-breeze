#!/bin/bash
#SBATCH --job-name=ffr-maint
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/feedback-maint-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
test -f "$1/DRIVER_PASS.json"
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/feedback_maintenance_v1/run.jl "$2" data/viscous2d-feedback512-v1-1925368/biased_plus/to_400/checkpoint.jld2 450 "$1/PASS.toml" "$3"
