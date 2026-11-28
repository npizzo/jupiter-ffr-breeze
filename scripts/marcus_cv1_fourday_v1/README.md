# Four-day CV1-like vortex persistence experiment

Authorized27November2026. This continues the successful reconstructed-atmosphere CPU2331883/GPU2331929 checks. It is a documented reconstruction of Zhang–Marcus CV1, not the author's original inputs and not spontaneous anticyclone formation.

## Fixed physical choices

Dry anelastic dynamics with Bannon-style linear buoyancy correction; full tilted rotation at23S, beta0, no convection or bottom heating. Prescribed anticyclone and background shear from the previously tested initializer, using Rv4600km, outer scale2000km, ellipticity1.5, vertical scales50/120km, centre8.2km. U_y=-9.1e-6/s in the central region. Positive relative vertical vorticity is anticyclonic here. Periodic horizontal domain105000x70000km, vertical domain-220..155km, stress-free impermeable lids. Same reconstructed reference.csv and profile_audit.json: Nmax0.02449/s; raster-temperature discrepancy9K RMS,31K maximum. Exact author constants/profile are not claimed.

384x384x256 Float64, WENO9 without closure. The vertical absorbing depth changes from the fixture37.5km to the documented31km; meridional absorption remains7000km. Both use the declared trial20000s timescale, damping toward the balanced zonal background; it is not claimed to be a recovered published coefficient. Record sponge overlap, compensation circulation and instantaneous sponge power. Passive dye starts in eight smooth zonal bands to make deformation visible; it does not affect dynamics.

## Tests and release

Real-driver CPU preflight exercises rest, shear and vortex branches, native buoyancy correction, all outputs, checkpoint readback and restart. Full source, profile, environment, Julia and pinned Breeze/Oceananigans sources are fingerprinted. GPU repeats the small checks, then runs384x384x256 through1800physical seconds. Adaptive advective CFL0.25; timestep at most4s and0.1/Nmax. Timing excludes first40steps; metadata separates initialization. No full-volume output at each diagnostic.

Pilot release requires finite positive state, mass-flux constraint<1e-8, Mach<0.8, final vertical/horizontal KE<0.3, horizontal KE within0.25..2 of its initial value, relative potential-temperature departure<0.25, and measured projected four-day cost<12GPUhours. These are exploratory startup/cost checks, not proof of vortex fidelity. The full-grid pilot must complete; a queue dependency alone cannot release production without GPU_PASS.toml and unchanged source receipt.

The dependent four-day run starts independently from the same initialized state. It is queued only after CPU preflight and depends on pilot success; it cannot run if the pilot fails. A12hour allocation checkpoints on orderly wall stop if the realized rate is slower than projected. No40/120day extension is authorized here.

## Outputs and interpretation

Target345600seconds (four Earth days). Save synchronized vorticity, tracer, temperature and vertical velocity nearz=-40,8.2,50km, every2160seconds, with actual grid heights relative to the one-bar level. Rolling checkpoints every21600seconds and at segment end. Record core-ellipse circulation, velocity extrema, Mach, thermal departure, horizontal/vertical KE and anelastic constraint. Cell-centred energy/sponge power are diagnostic proxies, not claims of exact discrete energy closure. All computation, rendering and raw data stay on Unity.

Assess initialized-vortex persistence, deformation and adjustment. Coherent survival does not establish formation or sign selection. The specified reference, shear, sponge and numerical departures remain part of interpretation.
