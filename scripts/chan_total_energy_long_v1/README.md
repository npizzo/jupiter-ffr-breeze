# Authorized total-energy integration to t=2000

User authorization: 24 November 2026, run the long case if the implementation validation passed. CPU2233881 and GPU2234115 passed on Breeze8273be89d9357ca2ac572347c45db8c647e3c45c. Moving closed-box energy residual was zero; constant-input residual2.03e-13; both formulation restart errors zero.

Physical parameters and original initial state match chan_clean_pair_v1. The thermodynamic variable is total energy. WENO9/no closure,256×256×80, Omega3, latitude67.5N/full rotation, beta0, no jet, deterministic bottom heating and original cap/top cooling. No dry-RCE redesign and no reseeding. This uses the isolated new Breeze/Oceananigans environment; old runs remain unchanged.

The new dispatch gets an exact-driver CPU preflight and a GPU restart check inside one24-hour maximum allocation. The actual production model starts at t0 and advances directly to2000. At t5, record a full-grid native-energy check and continue without restarting. Every later diagnostic keeps enforcing residual<0.1% of gross bottom heat and native/reconstructed energy agreement1e-10 relative. External heat is sampled by diagnostic-cadence trapezoidal integration, so this is not a claim of exact RK-stage forced closure. Failure records a guard and stops; no hidden thermal compensation or automatic extension.

Output cadence unchanged: diagnostic1, three-depth vorticity/tracer5, checkpoint50, full volume250. Time units are nondimensional. Production GPU validation and first real advancement must be verified from logs, not inferred from submission. CPU preflight and source/environment fingerprints bind this exact version. Simulation, analysis and rendering stay on Unity.
