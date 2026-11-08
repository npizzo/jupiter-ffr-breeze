#!/bin/bash
# Run from the Unity project root. The two output roots must be new.
set -euo pipefail
root=${1:?Provide a new production output root}
preflight=${2:?Provide a new preflight output root}
[[ ! -e "$root" && ! -e "$preflight" ]] || { echo 'Output roots must not already exist' >&2; exit 1; }
mkdir -p logs
thermal=$(sbatch --parsable slurm/large40-thermal.sbatch "$preflight")
thermal=${thermal%%;*}
printf 'thermal %s\n' "$thermal"
probe=$(sbatch --parsable --dependency="afterok:$thermal" --kill-on-invalid-dep=yes slurm/large40-preflight.sbatch "$preflight")
previous=${probe%%;*}
printf 'probe %s\n' "$previous"
for index in {1..8}; do
    job=$(sbatch --parsable --dependency="afterok:$previous" --kill-on-invalid-dep=yes slurm/large40.sbatch "$index" "$root" "$preflight")
    previous=${job%%;*}
    printf 'segment %s job %s\n' "$index" "$previous"
done
