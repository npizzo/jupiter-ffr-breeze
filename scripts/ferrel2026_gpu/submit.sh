#!/bin/bash
# New queues only; optionally attach to the separately submitted setup job.
set -euo pipefail
root=${1:?new production root required}
preflight=${2:?runtime preflight root required}
setup_job=${3:-}
if [[ -e "$root" ]]; then echo 'Production root exists; refusing duplicate submission' >&2; exit 1; fi
if [[ -n "$setup_job" ]]; then
    [[ "$setup_job" =~ ^[0-9]+$ ]] || { echo 'Setup job must be an integer Slurm ID' >&2; exit 1; }
    mkdir -p "$preflight"
else
    if [[ -e "$preflight" ]]; then echo 'Use a fresh preflight root or provide its existing setup job ID' >&2; exit 1; fi
    mkdir -p "$preflight"
fi
mkdir -p "$root" logs
record="$root/submitted-jobs.tsv"
printf 'stage\tcase\tdays\tjob\tafterok\n' > "$record"
if [[ -z "$setup_job" ]]; then
    setup_job=$(sbatch --parsable slurm/ferrel2026-gpu-setup.sbatch "$preflight")
fi
printf 'runtime_setup\tall\t0\t%s\t\n' "$setup_job" >> "$record"
unit=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$setup_job" slurm/ferrel2026-gpu-unit.sbatch "$preflight")
printf 'unit\tall\t0\t%s\t%s\n' "$unit" "$setup_job" >> "$record"
cpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/ferrel2026-gpu-cpu.sbatch "$preflight")
printf 'cpu64\tclosure+restart\t0.0014\t%s\t%s\n' "$cpu" "$unit" >> "$record"
quick=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$cpu" slurm/ferrel2026-gpu-quick.sbatch "$preflight")
printf 'quick96\tvalidation\t0.014\t%s\t%s\n' "$quick" "$cpu" >> "$record"
pilot=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$quick" slurm/ferrel2026-gpu-pilot512.sbatch "$preflight")
printf 'pilot512\tmatched+axisymmetric\t0.25\t%s\t%s\n' "$pilot" "$quick" >> "$record"
fine=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$pilot" slurm/ferrel2026-gpu-pilot1024.sbatch "$preflight")
printf 'pilot1024\tb162\t0.25\t%s\t%s\n' "$fine" "$pilot" >> "$record"
previous=$pilot
for case in b0 b162 boundary162; do
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/ferrel2026-gpu.sbatch "$case" "$root" "$preflight")
    printf 'screen512\t%s\t2\t%s\t%s\n' "$case" "$job" "$previous" >> "$record"
    previous=$job
done
cat "$record"
