#!/bin/bash
#SBATCH --job-name=cai-B0-assessment
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/cai-B0-assessment-%j.out
set -euo pipefail
source ./local.env
source scripts/cai2022_benchmark_v2/environment.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/cai2022_calibration_v2/assess.py "$1"
if test -f runs/cai2022-benchmark-v2/B0_CALIBRATION_PASS.json; then
    python scripts/cai2022_benchmark_v2/release.py production "$2" B7 "$3"
fi
