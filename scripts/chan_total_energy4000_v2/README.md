# chan_total_energy4000_v2: native continuation, t=2000–4000

User authorized on 24 November 2026. Parent: `chan-total-energy-long-v1-run-2237362`. Equations: compressible; thermodynamic formulation: total energy at PR1165 pin8273be8.

Physics, grid256×256×80, Ω3 full axis67.5°N, periodic4×4×1 box, beta0, no jet, uniform constant bottom heat, local cap relaxation, WENO9 and no closure are unchanged. Native checkpoint continuation; no reseeding. The original potential-temperature extension is explicitly a flagged diagnostic despite its previously measured29.6% heat-ledger deficit and partly unstable cap. Its endpoint is4000, not40000.

CPU preflight runs now using the actual driver, including all output and restart paths, and binds every executable/configuration/dependency. A dependent short CPU job runs only after both donor and preflight succeed. It checks completion, finite positive state, mass/constraint/Mach bounds, configuration agreement, checkpoint iteration and SHA256, and storage. Total-energy release additionally requires the existing0.1%-of-bottom-heat energy check and reconstruction agreement. The anelastic compressible-energy proxy is not used as a conservation gate. Cap stability and thermal drift are recorded, not claims of equilibrium. A donor failure, early wall limit, source mismatch or missing checkpoint blocks release.

Release is locked and idempotent. The GPU job verifies source and donor fingerprints again, exercises serialized restart on GPU, then checks the native full-grid restart endpoint against donor KE, energy, minimum density and temperature before advancing. Existing runtime guards remain active. No automatic further extension. Hard allocation24GPUhours per case. Scientific work stays on Unity.

Output intervals remain diagnostics1, movies5, fullvolume250, rolling checkpoint50. Continuation heat ledgers and profile-drift baselines reset at the segment start; donor endpoint and full prior ledger remain available for cumulative analysis. No alteration of active parent source, environment or outputs.

V2 corrects the native-check timing: it refreshes derived thermodynamics and verifies the checkpoint immediately after loading, before run! and before any timestep. V1 placed this check in a scheduled callback that first fired after one step and therefore skipped the exact-time condition. The production native-check branch now runs in both CPU and GPU preflight. This changes validation and diagnostic baselines only, not physical equations or parameters.
