#!/bin/bash
#SBATCH --job-name=ffr-2d-visc-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=logs/viscous2d-v1-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
out="data/viscous2d-v1-check-${SLURM_JOB_ID}"
"$JULIA_BIN" --project=. scripts/viscous2d_v1/initialize.jl "$out"
for case_name in unforced forced random_sign cyclonic_sign; do
 "$JULIA_BIN" --project=. scripts/viscous2d_v1/run.jl "$out/$case_name" "$case_name" "$out/initial.jld2" probe
done
tools/point4-analysis-env/bin/python scripts/viscous2d_v1/verify.py "$out"
