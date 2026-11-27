# Matched Cai B7 anelastic convection

Authorized27November2026 after measured anelastic/compressible/Boussinesq cost comparisons. This explicit new authorization supersedes the older requirement to wait for a compressible anticyclone before starting any anelastic run. It does not change the compressible B0 calibration or B7 release gate.

## Fixed setup

512x512x101 Float64, box4x4x1, periodic x/y, impermeable stress-free lids, vertical fixed rotation Omega2.56. Bottom heat flux0.001, top temperature1, conductivity0.00025. g8,R1,cp2.5. Rest-start with the same mean-free thermal perturbation amplitude1e-5 and seed314159 as compressible B7. No seeded velocity, beta, jet, cap, volumetric cooling or drag.

The actual discrete compressible reference density and pressure are copied into AnelasticDynamics. WENO9 and Cai grid-scaled Smagorinsky match the benchmark. Momentum uses rho0-weighted full deviatoric stress; no SGS scalar diffusion. The potential-temperature equation receives conduction and viscous dissipation as heat. These are corresponding approximations, not identical equations. Compressible evolves total energy; anelastic freezes reference density and eliminates acoustic modes. The small initial thermal perturbation maps to temperature differently in the two formulations.

Native anelastic potential temperature is used for stratification diagnostics. Reconstructed internal+kinetic+gravitational energy is a diagnostic proxy only: its residual is not labeled a compressible total-energy closure error. Every-step boundary heat quadrature is retained. Record maximum Mach, mean thermal drift, maximum relative potential-temperature departure and inferred density departure at fixed reference pressure. Mach>0.3 or potential-temperature departure>0.3 writes an approximation warning; these exploratory thresholds do not constitute validated asymptotic error bounds. The existing finite-state, positivity, mass-constraint and Mach1.2 stopping checks remain. Significant reference departure limits physical comparability even if numerically stable.

## Duration, cost and output

Scientific targett5000. First stable A10080GB allocation at most24hours, checkpoint before its wall limit. No automatic multi-day extension; assess attained time, thermal evolution, approximation warnings and cost before continuation. Benchmark rates imply roughly4–12A100-days to5000 depending on flow strength; this is not a promise of an LSV by that time. Stop time remains5000 even when the first segment ends earlier.

Diagnostics and synchronized three-depth fields nearz0.2,0.5,0.9 every10time units, including vorticity, passive tracer, temperature and vertical velocity. Full fields every1000; rolling native checkpoints every100 and at orderly segment end. Initial fields saved. All outputs and processing stay on Unity, in the existing30-day scratch workspace; no old data deleted.

## Launch checks

Full-source parse and real-driver CPU preflight test conductive equilibrium, matching reference/thermal perturbation, exact viscosity formula, viscous heat return, native restart, output readback, rotating/nonrotating and moving branches, and the approximation-warning path. Fingerprints bind source, configs, batch/environment, Julia and pinned solver code. The GPU allocation repeats the small restart/output test, then a target-gridt2check before production. The scientific run independently starts from rest after this check; it does not resume its pilot. Existing compressible jobs and shared dependencies remain untouched.
