# Clean anelastic continuation, t=2000–4000

User-authorized continuation of job 2214011_1 from its native final checkpoint. Physics, 256×256×80 grid, Ω=3, beta=0, no jet, WENO9, no closure, heating/cooling and all output intervals are unchanged. This is not a dry-RCE redesign.

The model file is byte-identical to chan_clean_pair_v1. Driver selects the anelastic branch and checkpoint pickup. Only target time and restart settings change. Metadata identifies the checkpoint; the diagnostic energy-proxy ledger and profile-drift reference reset at t=2000. That energy proxy is not a conserved anelastic invariant.

CPU preflight runs the actual driver and all output/restart paths. The receipt binds sources, environment and native donor checkpoint. GPU validation precedes production in one 12-hour allocation. Release is idempotent; no automatic extension. All scientific computation stays on Unity.
