# Isolated CUDA deployment repair for the unchanged v2 experiments

The v2 CPU validation passed all17 checks. Its first GPU gate then failed during
CUDA initialization because the original environment had been precompiled on a
node without an NVIDIA driver. No production simulation started in that chain.

[CUDA's official installation documentation](https://cuda.juliagpu.org/stable/installation/overview/#Precompiling-CUDA.jl-without-CUDA)
explains that runtime/compiler artifacts are selected during precompilation and
that an explicit toolkit preference is needed when precompiling without CUDA.
This repair sets both CUDA_Runtime_jll and CUDA_Compiler_jll to12.2 in a new
environment before imports. The corresponding UUIDs are listed in Project.toml
extras so Preferences can find them. Setup runs on an allocated L40S and must
verify CUDA functionality and the selected runtime before the chain proceeds.

The new environments/gamma2026_gpu environment copies the existing resolved
Manifest on Unity. Setup verifies package names, versions, tree hashes and
source revisions remain unchanged. The original environment and depot are not
modified. A new depot-gamma2026-gpu cache is first in JULIA_DEPOT_PATH; the old
depot supplies existing packages/artifacts as a fallback, with new writes going
to the new depot head. Runtime and compiler preferences invalidate the unusable
CUDA precompilation selection.

All numerical code, parameters, initial conditions, stopping conditions and
output settings remain the frozen scripts/gamma2026_v2 and configs/gamma2026_v2.
New Slurm wrappers select only the repaired environment and cache. The chain
reruns CPU validation and both GPU-grid gates before their dependent simulations.
A runtime report records the GPU, driver/runtime versions, Manifest/preference
hashes and unchanged-package check in the new output root's preflight directory.

Submit once on Unity to a new output root:

```sh
bash scripts/gamma2026_gpu/submit.sh data/gamma2026-v2gpu1
```

The script first submits a one-hour L40S setup job, then queues the existing
validation and eight-case sequence with successful-predecessor dependencies.
An optional second argument can refer to an already-submitted repair setup job
for that output root. It never modifies or cancels an earlier chain. Successful
runtime setup is deployment validation, not a scientific result.
