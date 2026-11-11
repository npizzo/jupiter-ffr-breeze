#!/bin/bash
#SBATCH --job-name=ffr-ferrgpu-quick
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ferrel2026-gpu-quick-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/ferrel2026_gpu/gpu_job.sh "$1" scripts/ferrel2026_v3/preflight.jl "$1" quick
