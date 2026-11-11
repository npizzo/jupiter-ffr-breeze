#!/bin/bash
#SBATCH --job-name=ffr-conv26-p1024
#SBATCH --partition=uri-gpu
#SBATCH --gpus=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/convection26-p1024-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/convection2026/gpu_job.sh scripts/convection2026/preflight.jl "$1" 1024
