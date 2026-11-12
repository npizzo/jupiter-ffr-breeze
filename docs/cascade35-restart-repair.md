# 35°N restart-validation repair — November 11, 2026

The CUDA repair passed, but latitude validation job 1447350 failed one of nineteen restart comparisons: the depth-integrated meridional mass streamfunction. This was a diagnostic tolerance error, not evidence of a blown-up flow. Native mass fluxes and all eighteen other checked fields passed; saved Float32 vorticity and tracer slices matched exactly.

A Unity audit (1459099) measured a streamfunction L2 difference of 3.9208144e-8 kg/(m s), relative difference 2.1202003e-8. Its native meridional mass-flux L2 difference was 1.6563240e-12 kg/(m² s). The original test incorrectly reused the same numerical absolute tolerance, 1e-10, for both quantities despite the extra length dimension in the integral. Horizontal velocity relative L2 differences were about 1.6e-14; this does not imply bitwise-identical restarts.

## Correction

The eighteen other checks retain their original tolerances. The streamfunction check now compares Psi/H using the same native-flux tolerance. In addition, it requires the actual streamfunction difference to be bounded by the actual flux difference propagated through the integral:

`norm(delta Psi) <= H * norm(delta native_rhov) + summation_roundoff_allowance`.

For uniform dz, the cumulative-sum operator's 2-norm is at most nz, so dz*nz = H gives this conservative bound. The summation allowance is `8*eps(Float64)*nz*H*max(norm(flux_a),norm(flux_b))`. This is an error-propagation bound rather than a fitted relaxation of all tolerances. The measured allowed bound was 1.9876027e-7 kg/(m s).

The Unity regression job **1459435 passed all four checks**, including rejection of an intentionally inconsistent streamfunction and a changed native mass flux. An earlier regression attempt, 1459406, passed the checks but failed while writing its mixed-type status record; that reporting error was corrected before resubmission.

Model equations, configuration, timestep algorithm, checkpoint mechanism and the 55°N source/chain are unchanged. Only the inactive 35°N validation and a new comparison helper were changed; the new helper is included in the campaign's source hashes. Original failed validation fields and submission records remain archived on Unity under `data/cascade35-restart-failed-20261111`.

## Resubmitted chain

Full GPU validation **1459466** must pass the restart checks, one-day production-resolution integration, actual initial-field matching, cap/tracer/heat-budget checks and throughput limit. **It completed successfully in 12 minutes 31 seconds.** All 21 restart checks, all 11 initial-field comparisons, the one-day production-resolution integration, stable cap, tracer, thermal-budget and throughput checks passed. The first 0–10-day production job **1459467 is RUNNING** on an NVIDIA L40S as checked November 11 at approximately 7:20 p.m. Eastern. The verification record and source/environment hashes are in [cascade35-restart-verified.json](cascade35-restart-verified.json).

The validation job alone used a scheduler time-limit override of 30 minutes (reduced from two hours to allow earlier backfill); production allocations are unchanged. The measured estimate is 1.515 wall hours per ten simulated days; this is an estimate, not a guaranteed completion time.

- Days 0–100: **1459467–1459476**.
- Day-100 comparison rendering: **1459477**, also waiting for the unchanged 55°N day-100 job **1447348**.
- Days 100–300: **1459478–1459497**.
- Day-300 rendering: **1459498**.

All jobs retain successful-completion dependencies. See [the complete submission record](cascade35-restart-v2.json). All scientific computation and rendering remain on Unity; no model fields were downloaded.
