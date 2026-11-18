#!/bin/bash
#SBATCH --job-name=ffr-maint-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/feedback-maint-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --startup-file=no -e 'function check(ex); ex isa Expr || return; ex.head in (:error,:incomplete) && error(ex); foreach(check,ex.args); end; foreach(p->check(Meta.parseall(read(p,String))),ARGS); println("PARSE_PASS")' scripts/feedback_maintenance_v1/model.jl scripts/feedback_maintenance_v1/run.jl scripts/feedback_maintenance_v1/check.jl
"$JULIA_BIN" --project=. scripts/feedback_maintenance_v1/check.jl data/viscous2d-feedback512-v1-1925368/biased_plus/to_400/checkpoint.jld2 "data/feedback-maintenance-check-${SLURM_JOB_ID}"
