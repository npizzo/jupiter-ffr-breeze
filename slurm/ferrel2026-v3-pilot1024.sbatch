#!/bin/bash
#SBATCH --job-name=ffr-ferr26v3-p1024
#SBATCH --partition=uri-gpu
#SBATCH --gpus=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel26v3-p1024-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/ferrel2026_v3/gpu_job.sh scripts/ferrel2026_v3/preflight.jl "$1" 1024
