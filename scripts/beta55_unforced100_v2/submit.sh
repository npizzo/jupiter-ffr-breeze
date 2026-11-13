#!/bin/bash
set -euo pipefail
root=${1:?new production root required}
preflight=${2:?new preflight root required}
[[ ! -e "$root" && ! -e "$preflight" ]] || { echo 'Refusing duplicate output roots' >&2; exit 1; }
mkdir -p "$root" "$preflight" logs
record="$root/submitted-jobs.tsv"
printf 'stage\tend_day\tjob\tafterok\n' > "$record"
previous=$(sbatch --parsable slurm/cuda122/beta55_unforced100_v2-preflight.sbatch "$preflight")
printf 'validation\t1\t%s\t\n' "$previous" >> "$record"
for segment in $(seq 1 10); do
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" --job-name="ffr-unforced55v2-d$((segment*10))" slurm/cuda122/beta55_unforced100_v2.sbatch "$segment" "$root" "$preflight")
    printf 'production\t%s\t%s\t%s\n' "$((segment*10))" "$job" "$previous" >> "$record"
    previous=$job
done
job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/beta55-unforced100-v2-gallery.sbatch "$root")
printf 'rendering\t100\t%s\t%s\n' "$job" "$previous" >> "$record"
cat "$record"
