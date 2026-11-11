#!/bin/bash
set -euo pipefail
root=${1:?new production root required}
preflight=${2:?new preflight root required}
if [[ -e "$root" || -e "$preflight" ]]; then echo 'Use fresh roots; refusing duplicate submission' >&2; exit 1; fi
mkdir -p "$root" "$preflight" logs
record="$root/submitted-jobs.tsv"
printf 'stage\tcase\tdays\tjob\tafterok\n' > "$record"
unit=$(sbatch --parsable slurm/ferrel2026-v3-unit.sbatch "$preflight")
printf 'unit\tall\t0\t%s\t\n' "$unit" >> "$record"
cpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/ferrel2026-v3-cpu.sbatch "$preflight")
printf 'cpu64\tclosure+restart\t0.0014\t%s\t%s\n' "$cpu" "$unit" >> "$record"
quick=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$cpu" slurm/ferrel2026-v3-quick.sbatch "$preflight")
printf 'quick96\tvalidation\t0.014\t%s\t%s\n' "$quick" "$cpu" >> "$record"
pilot=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$quick" slurm/ferrel2026-v3-pilot512.sbatch "$preflight")
printf 'pilot512\tmatched+axisymmetric\t0.25\t%s\t%s\n' "$pilot" "$quick" >> "$record"
fine=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$pilot" slurm/ferrel2026-v3-pilot1024.sbatch "$preflight")
printf 'pilot1024\tb162\t0.25\t%s\t%s\n' "$fine" "$pilot" >> "$record"
previous=$pilot
for case in b0 b162 boundary162; do
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/ferrel2026-v3.sbatch "$case" "$root" "$preflight")
    printf 'screen512\t%s\t2\t%s\t%s\n' "$case" "$job" "$previous" >> "$record"
    previous=$job
done
cat "$record"
