# Matched compressible total-energy preparation

Pinned Breeze PR1165 revision: `8273be89d9357ca2ac572347c45db8c647e3c45c`, updated 24 November 2026. Isolated Julia project and depot; Oceananigans0.118.1 required by the PR. Existing runs retain their original software.

The physical target remains the original clean Chan-style configuration: 256×256×80, box4×4×1, Omega3 at67.5N, full constant rotation, beta0, no jet, periodic horizontal boundaries, impermeable stress-free lids, bottom energy flux0.03125, top conductive temperature1 and local cap energy source -rho*(T-1)*S(z)/20. WENO9, no closure, original weak thermal-noise seed314159, no initial velocity. No stochastic forcing. No dry-RCE changes. Target comparison t0–2000 remains proposed, not submitted.

## What has been implemented

- Select total energy or potential temperature on the same new Breeze revision.
- Keep the same physical hydrostatic initialization and energy forcing; select the matching advected thermodynamic field.
- Diagnose theta from dry rho/T rather than accessing a theta prognostic field in the energy model.
- Integrate native rhoE for total energy; retain an independent staggered kinetic+internal+gravitational reconstruction and report their difference.
- Preserve original output, positivity, native checkpoint and restart paths. The external-source ledger remains trapezoidal at the diagnostic cadence; it is not advertised as exact RK-stage quadrature.

## Implementation validation only

The real driver has only preflight and gpu_validate modes; there is deliberately no long-run submission. It checks matched physical initialization, direct native energy forcing, a closed moving case, constant known energy input, original heating/cooling, output readback and serialized restart for both formulations. CPU uses the actual driver. A passing source/environment-bound receipt releases exactly one GPU validation allocation capped at20minutes; it cannot release production.

Closed/constant-source checks require relative energy error below1e-11; original-forcing check allows source quadrature error below0.1% of gross bottom heat. Reconstructed/native energy must agree within1e-11 relative. Actual measured errors are recorded in DRIVER_PASS.toml. These are small implementation checks, not a test of LSV formation or production-grid resolution. A developed-flow/full-grid performance and closure check remains necessary before committing to a long run.

## Environment history

Initial environment setup2233707 stopped before any simulation because the reused registry catalog ended at Oceananigans0.117.2. Setup_v2 copies the registry into the private depot and refreshes that copy. No shared registry, active environment or active model was changed. Validation CPU job2233881 depends on successful setup2233877. Results must be read from Unity; submission is not evidence that validation passed.
