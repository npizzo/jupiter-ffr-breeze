# Matched internal-deformation pair: completed results

Both Unity runs completed 25 initial turnover times (nondimensional t=2500). Scheduled diagnostics and movies end at t=2490, or 24.9 turnovers; the final scheduled output callback was not invoked at the floating-point stopping time. No endpoint is interpolated.

**There is no convection in these experiments.** The fluid is stably stratified, with no heating, cooling, beta or imposed jet. The aim is to test internal-deformation effects on pre-existing mixed-sign eddies, not convection-driven vortex formation or the compressible external Lamb mode.

| Case | L/L_D | Upper anticyclonic enstrophy share | Upper skewness | Final / initial KE |
|---|---:|---:|---:|---:|
| Weak deformation | 0.5 | 48.6% | +0.095 | 0.91 |
| Strong deformation | 2 | 43.3% | +0.560 | 4.95 |

Both began with 50% anticyclonic enstrophy. Last-five-turnover averages are 48.8% and 44.9%, respectively. The measure is the negative-vorticity contribution to squared vorticity at z/H=0.8984375, not vortex area or count. Positive vorticity is cyclonic for f>0. At the last saved time, the strong case also favors cyclonic squared vorticity at the other two sampled depths.

**This realization does not establish an anticyclonic preference.** Matching the velocity fields while reducing stratification changes their initial potential-energy reservoir: the strong case starts with approximately 16 times the small-amplitude APE proxy. Its KE grows while that proxy declines. Conversion of initial baroclinic energy is a plausible interpretation; no closed conversion budget was calculated. This is not a matched-total-energy comparison, and the outcome does not refute nonlinear shallow-water sign asymmetry.

## Configuration and checks

- Dry anelastic Breeze, 128×128×64, domain 8×8×1; periodic horizontally and impermeable stress-free lids.
- Vertical f=0.1, initial Ro=0.1, Re=2000; WENO9 and ordinary constant viscosity, no explicit scalar diffusion.
- Same initial random-phase first-baroclinic velocity field; seed 261130 selected using initial static stability only, before evolution. Balanced buoyancy and pressure fields differ with stratification.
- L=0.474079 is a velocity-gradient length, not a Fourier wavelength. U=0.00474079; initial turnover time L/U=100.
- L_D=c_1/f from the first internal mode. N²=0.0935561 and 0.00585015.
- Minimum recorded local N² stays positive: 0.08417 and 0.0009994. Peak Mach <0.0073, anelastic mass-flux constraint <1.9e-14, scalar conservation errors <6.6e-14.
- GPU startup and serialized restart checks passed, with zero restart field difference.
- KE plus small-amplitude APE proxy declines about 5.8% in both cases. This proxy is not an exact anelastic invariant; with viscosity present, its change is not an exact numerical energy residual.
- No resolution-convergence or inverse-transfer claim is made.

## Execution and publication

Unity jobs 2308454 and 2308455 completed successfully on A100 80 GB GPUs. Allocations lasted 9m29s and 8m50s including initialization and validation; integration took about 3.1m and 2.6m.

Unity rendering job 2310082 generated and fully decoded two 25-second movies. Each combines three-depth raw vorticity above absolute tracer concentration, with fixed scales shared between cases. The sampled heights are z/H=0.2421875, 0.6015625, 0.8984375.

**Movie publication remains pending.** The GitHub connector's encoded movie upload exceeded the automatic approval review budget; browser upload could not proceed because computer-use permissions were not granted. The finished assets and prepared gallery page are retained for publication. No broken movie links have been added to the public gallery.

See [summary.json](summary.json) for diagnostics and media checksums, [simulation source](../../scripts/deformation_pair_v1/), and [analysis/rendering source](../../scripts/deformation_pair_results_v1/). All field analysis, rendering and media decoding ran on Unity. Only finished publication assets were transferred.
