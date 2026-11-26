# Native continuation: chan-total-energy-seeded4000-v1

Authorized November 26: continue the seeded total-energy compressible run from t=2000 to t=4000. Donor job 2263614. No reseeding or parameter changes.

256×256×80; 4×4×1 periodic horizontal box; impermeable stress-free top/bottom; Ω=3 at 67.5°N; beta=0; no imposed jet; WENO9 with no closure; original bottom heat and local cap cooling retained. Total-energy environment is unchanged and pinned. Diagnostics every 1, planes every 5, volumes every 250, rolling checkpoints every 50 time units. Times are nondimensional.

Existing guards retained: positive density/temperature, mass error below 1e-7, Mach below 2, energy-ledger residual below 0.001 of accumulated bottom heating after startup. This does not assert low Mach or thermal equilibrium. The endpoint cap and mean-temperature drift are reported, not silently repaired.

CPU preflight exercises the actual driver and restart diagnostics, binds the source/environment fingerprint, validates the completed donor, then idempotently releases the one-GPU 24-hour job. GPU validation precedes production. The native endpoint is checked before advancing. A wall-limit stop preserves a checkpoint and is not reported as completion; no automatic extension.
