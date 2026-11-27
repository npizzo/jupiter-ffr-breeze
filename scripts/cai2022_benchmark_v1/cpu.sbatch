#!/bin/bash
#SBATCH --job-name=cai-benchmark-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/cai-benchmark-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cai2022_benchmark_v1/environment.sh
python scripts/cai2022_benchmark_v1/gate.py preflight "runs/cai2022-benchmark-v1/cpu-$SLURM_JOB_ID"
python scripts/cai2022_benchmark_v1/release.py validation "runs/cai2022-benchmark-v1/cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"
