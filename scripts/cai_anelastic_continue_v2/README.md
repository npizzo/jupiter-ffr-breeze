# Cai B7 anelastic continuation to t=5000

User authorized correction on27November2026: resume the interrupted target rather than treating an allocation boundary as completion.

The initial restart is checkpoint_iteration126736.jld2 at t=248.44410780816244 from job2361149. Model, initialization, B7 physical configuration and dependencies are copied unchanged from cai_anelastic_long_v1. Resolution512x512x101, Omega2.56, beta0, bottom flux0.001, same transport and numerical guards. Native momentum, potential-temperature, tracer, clock, timestep and integrator state are restored. No thermal reset or new initial eddies.

Each segment requests one A100/H10080GB,4–12hours of scheduler backfill. On a normal wall-limit stop it checkpoints, verifies positive progress and source/receipt integrity, and submits exactly one dependent successor. The chain ends at t=5000. A failure, guard stop, missing checkpoint or ambiguous submission halts the chain without automatic retries. No extension beyond5000 is authorized. Other experiments are untouched.

Real-driver CPU preflight advances the production continuation path, reads its outputs, compares a serialized restart against uninterrupted evolution, checks transport/reference matching and exercises warning output. A controller self-test rejects failures, missing checkpoints, changed targets and no-progress outcomes. The receipt binds all code, batch files, configuration, pinned libraries, environment and Julia. Each GPU segment repeats a tiny restart check and verifies the actual checkpoint checksum, restored time and iteration before continuing the real state.

Diagnostics and movie intervals stay10, full volumes1000, rolling checkpoint interval100. As before, the reconstructed anelastic energy is a diagnostic proxy, not a compressible conservation identity. Integrated heat/proxy changes and profile drift restart their accounting at each allocation; these scopes are explicit in metadata. Maximum potential-temperature departure continues to use the original fixed reference. Approximation warnings remain visible and do not silently change physical parameters or stopping thresholds.

All scientific work and raw fields remain on Unity. Machine-readable state lives in runs/cai-anelastic-continue-v2/chain-state.json. Submission goes through chain.py start after LAUNCH_PASS.json; repeated start calls cannot create duplicate chains.

Revision2 uses the live scheduler deadline rather than potentially stale SLURM_JOB_END_TIME. It records both values. It also refreshes derived fields immediately after checkpoint loading, before segment reference/energy diagnostics. Previous continuation energy-proxy residual and segment profile-drift columns sampled stale initialization fields and must not be interpreted as cumulative physical budgets. Physics and prognostic checkpoint state are unchanged. The original pending successor is held during validation and replaced only after the new job is submitted.
