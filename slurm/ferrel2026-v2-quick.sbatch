#!/bin/bash
#SBATCH --job-name=ffr-ferr26v2-quick
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ferrel26v2-quick-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/ferrel2026_v2/gpu_job.sh scripts/ferrel2026_v2/preflight.jl "$1" quick
