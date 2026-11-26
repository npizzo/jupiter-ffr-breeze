# Fixed Omega=5 rest-start total-energy experiment

Authorized November 26, 2026: new unseeded compressible total-energy integration from t=0 to t=2000. Omega=5 remains fixed throughout; this is not a continuation or a rotation ramp.

Physical control is chan_total_energy_long_v1 (Omega=3). The only changed physical input is Omega. Same initial weak temperature perturbation and seed 314159; zero initial velocity; 256×256×80, 4H×4H×H, latitude 67.5N/full tilted rotation, beta=0, no jet, dry ideal gas, bottom flux0.03125, local upper cooling time20, WENO9, no explicit closure. Same pinned total-energy formulation, thermal switches and stopping tolerances. The split-explicit solver retains its rotation-dependent time-step bound.

An output-only addition saves cell-centered vertical velocity at the same three actual video heights every5 nondimensional time units. Vorticity, tracer and temperature slices remain every5; diagnostics1, checkpoints50, full volumes250. CPU/GPU readback validates the new planes against native face velocities.

Run the exact-driver CPU preflight before idempotent GPU submission. The fingerprint covers source/config/scripts/environment and pinned dependencies. GPU restart equivalence is checked inside the allocation before the production run. Existing mass/positivity/Mach/native-energy guards remain active; t5 full-grid energy check continues enforcing its threshold thereafter. Maximum allocation24hours, no automatic extension. All computation and output stay on Unity.
