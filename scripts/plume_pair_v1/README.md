# Isolated plume / interface pair

Authorized 19 November 2026 after the existing-run plume audit. This is a bounded local mechanism test, not an LSV-formation simulation.

- Dry anelastic Breeze at the existing pinned version, Float64, WENO5 and explicit SSPRK3. Constant 55°N f; no beta, shear, random eddies, or moisture.
- Local 200 x 200 x 40 km domain, representing original heights 70–110 km. Grid 96 x 96 x 400: horizontal spacing 2.083 km, vertical spacing 100 m. Horizontal periodicity replaces the large channel walls; top/bottom remain impermeable and free-slip.
- Initial velocity is zero. Initial potential-temperature anomaly is `2 K * cos²(pi*r/80km) * cos²(pi*(z-75km)/10km)` inside r<40km and 70<z<80km, zero outside. Passive mass fraction equals the same shape without the Kelvin amplitude. No random seed is used.
- The initial draft proposed 0.02 K; this was increased to 2 K and disclosed before submission because the weak, broad pulse risked sub-grid penetration of the strong interface. This amplitude is deliberately not calibrated to the earlier convective field.
- Stable case: N²=2.5e-5 s^-2 above 80 km. Neutral control: N²=0. Hydrostatic reference rebuilt in both cases at the same bottom pressure, calculated from the original neutral adiabat at 70 km. Gravity, dry H–He gas constant and heat capacity unchanged.
- Unlike the original steady-heating runs, only finite initial heat is supplied. Bottom flux, cap relaxation, global heat controller and sponge are absent. Ordinary viscosity 1000 m²/s, original thermal diffusivity taper (1000 to 1 m²/s over original 60–80 km), tracer diffusivity 100 m²/s retained. These are important changes, not a smaller copy of the turbulent campaign.
- Initial velocity Rossby number is zero. Reported evolving Ro is domain horizontal RMS speed/(f*40 km), not the earlier 1200-km spectral-eddy definition. It is descriptive, not held fixed.
- Two simulated hours each; maximum timestep 2 s, advective CFL target 0.2. Hard Slurm budget two A100-hours per case, four combined; model wall limit 6600 s. No automatic extension/retry or additional parameter sweep.
- Save synchronized vorticity/tracer at three heights every minute, vertical profiles every minute, full states every ten minutes, checkpoints every thirty minutes. Render on Unity.

## Budget and validity

The diagnostic follows density-weighted layer-mean vertical vorticity in 70–80, 80–84 and 84–110 km. It applies the same linear centered curl to momentum increments and native momentum tendencies, using the exact three SSPRK3 weights every timestep. Native vertical advective flux is retained explicitly; horizontal advection is the difference between total and vertical advection. Coriolis and viscosity are separately retained. Horizontal forcing is absent; pressure-gradient curl is explicitly checked for each layer. Partial-layer horizontal divergence is **not** assumed zero. Cumulative budgets are saved every minute.

Exact closure of this discrete observable does not validate a continuum stretching/tilting decomposition or a convection closure. The reconstructed conditional stretching estimate is labeled separately.

The small CPU validation checks the curl, reference state, nonzero-plume native partition/step closure for both controls, and an actual serialized velocity restart. The same stage checks run on the production GPU. Native ledger history restarts from the segment origin; it is not silently carried across a model restart. Production is not released if the implementation check fails.

Fidelity screening reports vertical-velocity power at wavelengths shorter than four horizontal cells, characteristic w/N in vertical cells, tracer penetration, and tracer/wave proximity to boundaries. Fewer than eight w/N cells or over 10% near-grid w power are warning thresholds, not convergence criteria. A raised-lid follow-up is conditional on evidence of boundary contamination; it is not automatically queued. Pressure and waves can communicate boundary effects before tracer reaches the lid. Interpret the stable/neutral comparison over its useful pre-contamination interval.

Do not infer LSV formation, an inverse cascade, or a null physical mechanism from this short pair. If it resolves anticyclonic production at the interface, the next question is whether that vorticity persists and organizes—not whether the initial pulse makes a GRS.
