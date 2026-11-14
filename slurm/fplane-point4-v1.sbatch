#!/bin/bash
#SBATCH --job-name=ffr-point4
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/fplane-point4-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
if [ "$1" = array ]; then
 mapfile -t cases < scripts/fplane_point4_v1/case-order.txt
 "$JULIA_BIN" --project=. scripts/fplane_point4_v1/cuda_launch.jl campaign "$2" "${cases[$SLURM_ARRAY_TASK_ID]}"
else
 "$JULIA_BIN" --project=. scripts/fplane_point4_v1/cuda_launch.jl "$@"
fi
if [ "$1" = preflight ]; then
 tools/gallery-env/bin/python scripts/fplane_point4_v1/test_diagnostics.py "$2/preflight"
 cp "$2/preflight/GPU_PASS.toml" "$2/preflight/PASS.toml"
fi
