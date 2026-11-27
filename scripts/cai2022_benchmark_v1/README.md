# Cai 2022 B0/B7 benchmark, v1

Authorized 2026-11-27. DOI 10.3847/1538-4357/ac3b4b, Table 1 and Appendix A.
B0 is calibration before B7 production. No anelastic run until a compressible anticyclone is established.

Dry total-energy ideal gas: R=1, cp=2.5, cv=1.5, gamma=5/3, H=1, g=8.
Initial T=5-4z, rho=T, p=T² (m=1, eta=4); conductive equilibrium unstable to convection.
Box 4×4×1; production 512×512×101 finite-volume cells. Paper reports grid points, so vertical discretization is not bitwise identical.
Periodic x/y; impermeable stress-free lids. Ttop=1; Fbottom=.001; constant conductivity .00025.
Adiabatic conductive fraction=K*g/cp/F=.8. No cap, volumetric cooling, imposed shear, beta, imposed drag or inserted vortex.
Rotation vertical, fixed Omega=0 (B0) or 2.56 (B7). Same zero initial velocity and seed314159 mean-free thermal perturbation amplitude1e-5.

Cai Eq A7: nu=.28² sqrt(dx dy) dz sqrt(2 sigma:sigma). Native constant Smagorinsky (NOT SmagorinskyLilly) with adjusted coefficient exactly matches this grid factor. No stability limiter, SGS scalar diffusion, entropy diffusion, or scalar Prandtl prescription; heat uses constant conductivity.
Momentum stress is 2 rho nu [sym grad u - div(u)I/3]. Total-energy forcing includes divergence(u dot stress), not an additional viscous-heating source. Kinetic dissipation returns to internal energy through diagnostic reconstruction from total energy.

Declared numerical departures: full nonlinear conservative stress retained; Cai reports dropping unspecified higher-order fluctuating-viscosity terms. WENO9 finite-volume advection instead of horizontal Fourier spectral method; Breeze split-explicit stepping rather than Cai semi-implicit method. The paper does not fully identify those omitted viscous terms or a final stopping time/perturbation recipe. These are a cross-code reproduction attempt; a negative result alone does not refute Cai.

Targets: B0 t=3000 (3 nominal thermal times), B7 t=5000. First production allocation limited to24hours, checkpoint before wall limit. No automatic long extension. B0 must show late-time top/bottom heat-flux balance and stable velocity statistics; compare reported rms .066 and convective flux ~.2F. B7 must form a persistent warm-core anticyclone with clockwise circulation, coherent across depth; a cyclone companion is expected and anti-enstrophy>50% is NOT required. Thermal relaxation and vortex saturation are separately assessed.

Release: complete-driver CPU equations/output/restart checks, then GPU equations/restart and 512²×101 memory/output/cost check to t=2. Validation receipt binds all code/configs/environment and solver source. No production release from stability alone. B7 additionally requires B0_CALIBRATION_PASS.json. B0 release after measured target-grid cost reviewed. No changes to frozen source after submission.

Outputs remain on Unity scratch workspace jupiter-cai-benchmark-20261127 (30day lifetime; renew/preserve results before expiry). Source and receipts stay in /work. Three actual heights near .2,.5,.9; save vorticity, passive tracer, T,w plus depth mean velocity every10time units. Full fields every1000; rolling checkpoints100. Energy input sampled every timestep. Full source/paper provenance in metadata. No old data deleted.
