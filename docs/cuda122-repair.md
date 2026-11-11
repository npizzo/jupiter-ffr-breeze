# CUDA repair and latitude campaign restart — November 11, 2026

Both previous latitude validation jobs failed before integration: CUDA runtime/compiler artifacts had been selected during precompilation without a GPU driver. The dependent production and rendering chains were cancelled. This was not a model instability.

The repair explicitly selects CUDA runtime and compiler 12.2 through a scoped Julia preference environment. A separate writable depot isolates compiled caches. The original Project.toml, Manifest.toml, model source, configurations and frozen extension source remain unchanged. This follows [CUDA.jl runtime selection](https://cuda.juliagpu.org/stable/installation/overview/) and [Julia preference inheritance](https://juliapackaging.github.io/Preferences.jl/stable/).

**GPU verification passed November 11 at 3:21:53 p.m. Eastern (20:21:53 UTC).** Job **1447116** completed successfully in 4 minutes 8 seconds on `uri-gpu010`, using an NVIDIA L40S, CUDA runtime 12.2.0, and compiler 12.2.140. The Float64 GPU kernel/reduction and three-step Breeze anelastic integration passed. The PASS record includes environment and probe hashes, all rechecked successfully after completion. Every restarted simulation checks that record, the hashes, and the actual GPU runtime before entering the unchanged driver. The full campaign restart, initial-condition matching, thermal/cap/tracer and throughput checks are still required.

| Campaign | Validation | Day 100 | Rendering |
|---|---:|---:|---:|
| 55°N heated beta plane | 1447338 | 1447348 | 1447349 |
| 35°N heated beta plane | 1447350 | 1447360 | 1447361 |

The 35°N extension continues with jobs **1447362–1447381** through day 300; its rendering job is **1447382**. The day-100 latitude comparison waits for both latitude chains. Rendering retains synchronized three-depth vorticity and passive tracer output. All computation and rendering stay on Unity.

The restart contains 45 dependent jobs plus the shared CUDA probe. The CUDA dependency is now satisfied. As checked at approximately 3:36 p.m. Eastern, both latitude validation jobs (1447338 and 1447350) are pending scheduler priority with no remaining dependency. Production remains gated on each latitude's full validation; successful CUDA verification does not establish full scientific validation. Full job/dependency records are in [cuda122-restart.json](cuda122-restart.json).

Failed, field-free output directories and original submission records were archived under `data/cuda122-failed-attempt-20261111` on Unity. Original output paths are reused only after confirming that no saved simulation frames exist, preserving the existing rendering/reference paths. No simulation fields were downloaded.
