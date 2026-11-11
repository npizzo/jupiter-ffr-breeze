#!/bin/bash
set -euo pipefail
root=${1:?new production root required}
preflight=${2:?new preflight root required}
if [[ -e "$root" || -e "$preflight" ]]; then echo 'Use fresh roots; refusing duplicate submission' >&2; exit 1; fi
mkdir -p "$root" "$preflight" logs
record="$root/submitted-jobs.tsv"
printf 'stage\theat_W_m2\tend_day\tjob\tafterok\n' > "$record"
unit=$(sbatch --parsable slurm/convection2026-unit.sbatch "$preflight")
printf 'unit\tall\t0\t%s\t\n' "$unit" >> "$record"
pilot=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/convection2026-pilot512.sbatch "$preflight")
printf 'pilot512\tall\t0.25\t%s\t%s\n' "$pilot" "$unit" >> "$record"
fine=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$pilot" slurm/convection2026-pilot1024.sbatch "$preflight")
printf 'pilot1024\t162\t0.25\t%s\t%s\n' "$fine" "$pilot" >> "$record"
previous=$pilot
# Serialize production to avoid consuming three GPU allocations simultaneously.
# Fine pilot is an independent refinement check, not proof of convergence.
for segment in $(seq 1 5); do
  for heat in 0 54 162; do
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/convection2026.sbatch "$heat" "$segment" "$root" "$preflight")
    printf 'production512\t%s\t%s\t%s\t%s\n' "$heat" "$((2*segment))" "$job" "$previous" >> "$record"
    previous=$job
  done
done
cat "$record"
