#!/bin/bash
set -euo pipefail
root=${1:?new production root required}
preflight=${2:?new preflight root required}
if [[ -e "$root" || -e "$preflight" ]]; then echo 'Use fresh output roots; refusing duplicate submission' >&2; exit 1; fi
mkdir -p "$root" "$preflight" logs
record="$root/submitted-jobs.tsv"
printf 'stage\tcase\tend_day\tjob\tafterok\n' > "$record"
unit=$(sbatch --parsable slurm/cascade100-unit.sbatch "$preflight")
printf 'unit\tboth\t0\t%s\t\n' "$unit" >> "$record"
previous=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/cascade100-preflight.sbatch "$preflight")
printf 'validation\tboth\t1\t%s\t%s\n' "$previous" "$unit" >> "$record"
for index in $(seq 1 20); do
  job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/cascade100.sbatch "$index" "$root" "$preflight")
  if (( index % 2 == 1 )); then case_name=cascade100_f_heated; else case_name=cascade100_f_control; fi
  day=$(( ( (index - 1) / 2 + 1 ) * 10 ))
  printf 'production\t%s\t%s\t%s\t%s\n' "$case_name" "$day" "$job" "$previous" >> "$record"
  previous=$job
done
cat "$record"
