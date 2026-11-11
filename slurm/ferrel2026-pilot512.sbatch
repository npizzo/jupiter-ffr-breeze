#!/bin/bash
#SBATCH --job-name=ffr-ferr26-p512
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel26-p512-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/ferrel2026/gpu_job.sh scripts/ferrel2026/preflight.jl "$1" 512
