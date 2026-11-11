# Isolated CUDA runtime deployment

Ferrel revision 3 passed its 104 initialization checks and Unity CPU native-closure/restart gate. Its GPU quick job `1405390` then stopped before time stepping because CUDA's JLLs had been precompiled without an NVIDIA driver and had no appropriate runtime selection. This was a runtime initialization failure, not a failed physical-model or numerical-stability test.

Runtime setup **1409245 passed** on an NVIDIA L40S (host driver 580.178.04, CUDA runtime 12.2.0). The original and isolated manifests have the identical SHA-256 `0ed2ee7d2e5e72a5558acdda0a87adc921adc7609f408ffafe81111eb9fb567d`; the Float64 device execution probe passed. Initialization **1409338** then passed all 104 checks in the repaired environment. CPU model validation **1409339** and GPU quick **1409340** also passed native momentum/mass closure, restart, symmetric controls, tracer and cap checks. The GPU restart comparison passed all 17 checks; sampled GPU memory peak was 1.54% of capacity on this small diagnostic grid. This is not the production memory estimate. Full-grid pilot **1409341** is running; only its successful completion releases 48-hour cases **1409343–1409345** and the independent refinement **1409342**.

The repair leaves `scripts/ferrel2026_v3/`, `configs/ferrel2026/`, the root Julia environment and other running campaigns unchanged. All queue wrappers use `--project=environments/ferrel2026_gpu` and an isolated first depot, `depot-ferrel2026-gpu`, with the existing root depot as a read-through fallback.

On an allocated Unity GPU node, setup preserves byte-exact copies of the root Project and Manifest. The isolated effective Project adds only the runtime and compiler JLL UUIDs under `[extras]`, allowing explicit preferences to resolve on Julia 1.12. Both CUDA runtime and compiler are selected as **12.2 before importing or precompiling CUDA**, including on subsequent CPU nodes.

The root environment resolves CUDA/CUDACore 6.4.2, CUDA_Runtime_jll 0.25.0+3 and CUDA_Compiler_jll 0.6.4+0. The working gamma environment established runtime 12.2 compatibility with the same CUDA family on Unity's L40S nodes, but its compiler JLL patch differs. Therefore this repair retains every root package version and tree, and independently verifies the root stack's runtime selection and actual Float64 device allocation, kernel execution and reduction. No package update is performed. Setup fails if package signatures, original files or the effective Project change unexpectedly.

Every subsequent launch checks the isolated Project, Manifest, preferences and runtime-guard fingerprints. The unchanged v3 scripts continue checking their own complete scientific-source and configuration hashes. The queue remains setup → initialization/negative checks → CPU native closure/restart → GPU quick gate → 512² screening → bounded 48-hour cases and the independent 1024² check. No failed older simulation variant is resumed or modified.

To attach the new queue to an already submitted setup job, use:

```bash
bash scripts/ferrel2026_gpu/submit.sh NEW_PRODUCTION_ROOT RUNTIME_PREFLIGHT_ROOT SETUP_JOB_ID
```

Omitting the third argument submits setup once, requiring both roots to be fresh. Scientific interpretation remains the short transient mechanism screening described in the Ferrel campaign documents; a working runtime does not establish Jovian fidelity.
