#!/bin/bash
set -euo pipefail
root=${1:?fresh extension root required}
original=${2:?original production root required}
preflight=${3:?original preflight root required}
previous=${4:?original successful day-100 dependency required}
[[ ! -e "$root" ]] || { echo 'Refusing duplicate extension root' >&2; exit 1; }
[[ "$previous" == 1435756 ]] || { echo 'Expected existing 35 N day-100 job' >&2; exit 1; }
mkdir -p "$root" logs
record="$root/submitted-jobs.tsv"
printf 'stage\tend_day\tjob\tafterok\n' > "$record"
for segment in $(seq 11 30); do
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" --job-name="ffr-beta35-d$((segment*10))" slurm/cascade35-extend.sbatch "$segment" "$root" "$original" "$preflight")
    printf 'production\t%s\t%s\t%s\n' "$((segment*10))" "$job" "$previous" >> "$record"
    previous=$job
done
job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/cascade35-300-gallery.sbatch "$root" "$original")
printf 'rendering\t300\t%s\t%s\n' "$job" "$previous" >> "$record"
cat "$record"
