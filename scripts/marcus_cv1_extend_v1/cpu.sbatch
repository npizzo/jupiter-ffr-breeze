#!/bin/bash
#SBATCH --job-name=marcus-extend-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/marcus-extend-cpu-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_extend_v1/environment.sh
python scripts/marcus_cv1_extend_v1/gate.py preflight "runs/marcus-cv1-extend-v1/cpu-$SLURM_JOB_ID"

python scripts/marcus_cv1_extend_v1/chain.py selftest baseline
python scripts/marcus_cv1_extend_v1/chain.py selftest fine
python scripts/marcus_cv1_extend_v1/allocation.py selftest
