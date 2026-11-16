#!/bin/bash
#SBATCH --job-name=ffr-2d-512-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/viscous2d-refine-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
out="data/viscous2d-refine-v1-${SLURM_JOB_ID}"
"$JULIA_BIN" --project=. scripts/viscous2d_refine_v1/initialize.jl data/viscous2d-v1-check-1846256/initial.jld2 "$out"
"$JULIA_BIN" --project=. scripts/viscous2d_refine_v1/run.jl "$out/unforced" unforced "$out/initial.jld2" probe
tools/point4-analysis-env/bin/python scripts/viscous2d_refine_v1/analyze.py data/viscous2d-v1-check-1846256/unforced "$out/unforced" "rendered/viscous2d-refine-v1-${SLURM_JOB_ID}"
