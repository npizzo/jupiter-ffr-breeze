# Polar convection resolution and heating pilot

This new, isolated dry Breeze campaign tests resolution, timing and bottom-heating response in a regional patch at **85°N**. It is not the full polar gamma-plane vortex-assembly experiment or a Ferrel-cell simulation. The local f-plane uses Jupiter's rotation rate and beta=0; its horizontal variation in Coriolis frequency is deliberately omitted at this stage.

The 6,000 × 6,000 × 120 km channel has periodic zonal boundaries, free-slip meridional and vertical walls, no imposed jet/shear, and identical full-depth random eddies of central wavelength 1200 km and Rossby number 0.1. Changing latitude from the earlier 55° runs changes f and the speed implied by fixed Ro. This is not an identical-initial-energy comparison to that campaign.

| Stage | Grid | Heat input | Duration |
|---|---|---|---|
| Initialization and spectral checks | two-level 512²/1024² initializer tests | all cases | no integration |
| Restart verification | 64² × 96 |162 W/m² |1200 s, continuous versus resumed |
| GPU construction and throughput |512² × 96 |0,54,162 W/m² |6 hours each |
| Selected fine pilot |1024² × 96 |162 W/m² |6 hours |
| First conditional science batch |512² × 96 |0,54,162 W/m² |10 days each, five 2-day segments |

The short tests establish execution and throughput, not convergence of mature convection. Production requires passing source-hash, initialization, restart, thermal, tracer, cap and output checks, plus estimated 2-day runtime below 4 wall hours. Allocations allow6 hours. The 512 grid uses one L40S and 96 GB host memory; the 1024 pilot requests one H100 and 128 GB host memory. Per-process GPU memory is sampled every 15 seconds, matched to its GPU UUID, and must remain below 95% of the device capacity. This is sampled process residency, not a guarantee of 5% transient headroom or a measurement of allocated array bytes; missing samples fail the resource gate. An allocation or out-of-memory failure blocks that job's dependent chain; no automatic CPU fallback occurs.

The first ten-day matched batch is intentionally bounded. Longer runs, ensembles, mechanical forcing matched to measured energy input, alternative heat-removal profiles, and moist extensions follow analysis of these results. The 1024 pilot runs independently after the 512 gate and does not release fine-grid production.

The thermal treatment is copied from the validated cascade campaign: neutral lower layer, stable cap above 80 km, one-day cap relaxation fully active above 90 km, vertical-velocity sponge above 108 km, and globally uniform fractional potential-temperature feedback to balance net imposed enthalpy input. These controls can influence circulation and must be varied before a transport mechanism is attributed to Jupiter. Bottom pressure remains 20 bar. No moist physics or calibrated atmospheric radiation is included.

Saved fields retain vorticity and passive tracer at three actual levels nearest 20, 60, 100 km, vertical velocity and potential temperature at those levels, depth-mean velocities/stresses, thermal budgets and spectral-transfer diagnostics. Additional cell-centered zonal means and density-weighted covariance fluxes versus latitude and height are saved for horizontal/vertical momentum, heat and tracer transport. They are not native discrete face fluxes or a closed Ferrel-cell momentum budget; residual/TEM analysis, full tendency accounting, jet forcing and boundary sensitivity are deferred.

All integration, validation, diagnostics and rendering run on Unity. Existing campaigns and their source files are unchanged. The new driver and dependencies are frozen under `scripts/convection2026`; completed gates reject changed code or configurations.

Submit on Unity with fresh roots:

```sh
bash scripts/convection2026/submit.sh data/convection2026-v1 data/convection2026-preflight-v1
```

Every submission is recorded in the new production root's `submitted-jobs.tsv`. Submission is not evidence that validation passed or simulations completed. Three-depth movies will accompany any scientific results; two-dimensional QG outputs elsewhere cannot supply fictitious depth panels.
