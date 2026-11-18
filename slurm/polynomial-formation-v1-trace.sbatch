#!/bin/bash
#SBATCH --job-name=ffr-poly-check-trace
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:04:00
#SBATCH --output=logs/polynomial-formation-trace-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. -e 'println("TRACE_JULIA_STARTED");flush(stdout);@time include("scripts/polynomial_formation_v1/model.jl");println("TRACE_MODEL_INCLUDED");flush(stdout);@time d=load("data/viscous2d-refine-v1-1848240/initial.jld2");println("TRACE_INITIAL_LOADED");flush(stdout);@time m,b=setup_noadvection(d["u"],d["v"]);println("TRACE_NOADV_MODEL_BUILT");flush(stdout);@time q=native_advection_curl(m,b.base.op);println("TRACE_ZERO_ADV=",maximum(abs,q));flush(stdout)'
