#!/bin/bash
#SBATCH --job-name=ffr-totalE-env
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=logs/chan-total-energy-v1-setup-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy_v1/environment.sh
"$JULIA_BIN" --startup-file=no scripts/chan_total_energy_v1/setup.jl
