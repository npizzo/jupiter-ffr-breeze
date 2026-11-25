# Resume the original rest-start total-energy case to t=4000

Authorized November25 after cancellation of job2261539. Resume its last complete native checkpoint at t=2050 (iteration295540), replaying t2050-2096 in a new output directory. Preserve existing fields and checkpoint. No physical, numerical, resolution, forcing or formulation changes. Existing seeded job2263614 remains authorized separately through2000. No other simulations are launched.

The completed-checkpoint health gate verifies cancelled scheduler status, unchanged physical metadata, finite-state and energy checks, saved diagnostic values and checkpoint iteration/hash. The exact CPU/GPU restart validation and full-resolution native-state check remain mandatory. New source version preserves stopped driver. Existing native energy ledger thresholds remain, with the new segment ledger explicitly beginning afterrestart. Source revision and environment unchanged.
