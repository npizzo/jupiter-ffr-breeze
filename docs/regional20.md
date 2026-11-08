# Longer f-plane and physical beta-plane campaign

The five-day runs establish sustained convection, but their small visible structures do not establish the absence of an inverse cascade. This campaign extends time, adds a physically parameterized beta plane, isolates the new sidewall geometry, and collects direct diagnostics of transfer within the depth-averaged mode.

## Physical scales and interpretation

The completed fine-vertical run has a final-day volume horizontal RMS speed near **1.56 m/s**. Dividing the 1,200-km domain width by this speed gives **8.9 Earth days**. This is a reference domain-crossing time, not a predicted inverse-cascade time; five days is short for domain-scale reorganization. Twenty days provide a useful extension but do not guarantee equilibration or a large vortex.

Using retained Ω = 1.7585 × 10⁻⁴ s⁻¹ and central latitude 55°N gives f₀ = **2.88095774 × 10⁻⁴ s⁻¹**. With the nominal equatorial-radius approximation a = 71,492 km, beta = 2Ω cos(55°)/a = **2.82167001 × 10⁻¹² m⁻¹ s⁻¹**. The radius convention follows the [NASA/JPL planetary parameter table](https://ssd.jpl.nasa.gov/planets/phys_par.html); this is a spherical local approximation, not an oblate-coordinate atmospheric model. The total Coriolis variation across the 1,200-km channel is only **1.175%** of f₀. Beta is not artificially amplified.

With the explicitly chosen convention L_beta = sqrt(U/beta), the reference Rhines scale is about **744 km**. Numerical factors and definitions differ in the literature; this value is a scale estimate, not a predicted jet wavelength. It is comparable to the box, so this first experiment has little room for multiple Rhines-scale structures. A larger domain may eventually be needed. [Rhines (1975)](https://doi.org/10.1017/S0022112075001504) describes how beta introduces waves that can limit upscale growth and favor anisotropic flow; beta does not simply accelerate vortex merger.

Our evolved volume vorticity/f₀ is near 0.6 and small-scale vorticity remains horizontally resolution-sensitive. Longer runtime alone cannot establish the strongly rotation-constrained regime needed for some inverse-cascade mechanisms. [Maffei et al. (2021)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/on-the-inverse-cascade-and-flow-speed-scaling-behaviour-in-rapidly-rotating-rayleighbenard-convection/65EA7F62ED610264B8588035016CFBF3) examine inverse transfer and large vortices under strong rotational constraint, with sensitivity to forcing and domain geometry. Their asymptotic regime is not presumed here.

## Experiments

All cases use bottom heating 162 W/m², prescribed upper cooling −162 W/m², neutral initial lower stratification, the stable cap and existing diffusivity/viscosity. There is **no imposed cyclonic shear**. These enhanced heat fluxes remain idealized sensitivities. Runs start from their prescribed initial states, not old checkpoints.

| Task | Configuration | Grid | Earth days | Purpose |
|---|---|---|---:|---|
| 1 | periodic f-plane | 256 × 256 × 64 | 20 | extend the existing geometry |
| 2 | f-plane channel | 256 × 256 × 64 | 20 | isolate meridional wall effects |
| 3 | physical beta-plane channel | 256 × 256 × 64 | 20 | isolate beta against task 2 |
| 4 | finer f-plane channel | 512 × 512 × 64 | 2 | horizontal-resolution sensitivity |
| 5 | finer beta-plane channel | 512 × 512 × 64 | 2 | horizontal-resolution sensitivity and matched beta comparison |

The short fine-grid cases compare the same early time windows; they do not establish the resolution of late-time organization. Vertical resolution is retained at 64 levels because the completed vertical comparison changed bulk statistics by 2–4%, while horizontal refinement changed vorticity by roughly 52%. The scientific emphasis remains on resolving that horizontal sensitivity.

A linear f(y) cannot join smoothly across a meridionally periodic seam. Channel cases therefore use periodic x, bounded y in [−600,600] km, and the existing bounded vertical interval [0,120] km. North/south walls are impermeable and free slip; scalar boundaries have no normal flux. The f-plane channel has exactly the same geometry and seed as its beta partner. Wall-compatible toroidal/poloidal mass fluxes use odd/even reflection of their potentials, giving discrete anelastic mass conservation, zero wall-normal flow, and initial eddy Rossby 0.1. Channel and periodic initial fields necessarily differ, so tasks 1 and 2 are a geometry sensitivity, not an identical-field beta control.

## Diagnostics for inverse transfer

Hourly output now includes full-depth density-weighted horizontal velocities U and V, plus depth-dependent velocity covariances Rxx, Rxy and Ryy. These fields permit a barotropic-mode energy analysis without treating a single 20-km slice as representative of the full depth.

`python scripts/barotropic_transfer.py DATA/CASE --output ANALYSIS/CASE_transfer.npz` computes:

- spectra of the density-weighted full-depth horizontal velocity after horizontal solenoidal projection and two-thirds spectral truncation;
- nonlinear self-transfer from −(U · grad)U and cumulative flux Pi(K) = −sum[k≤K] T(k), with **negative flux indicating upscale transfer within this projected barotropic mode**;
- work of the depth-dependent Reynolds-stress divergence, reported separately from barotropic self-transfer;
- the fraction of energy in zonal means, projection losses, high-wavenumber truncation losses and the self-transfer conservation residual.

For the channel, even/odd reflection implements free-slip parity before the spectral calculation. Physical wavenumbers account for the doubled extension in y. Compare physical wavenumbers and statistical time windows, not raw shell indices across geometries.

This diagnostic is **not a full 3-D cascade measurement or exact discrete WENO energy budget**. It omits pressure, buoyancy, viscous and numerical contributions to a closed total-energy balance. The stress-work estimate also uses centered/interpolated fields. Projection and truncation losses must be inspected before interpretation. Consistent time-averaged upscale flux across a resolved scale interval would support an inverse transfer claim for this mode; growth in low modes or a large visible vortex alone would not. [Lemm, de Wit & Kunnen (2025)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/scalebyscale-kinetic-energy-flux-calculations-in-simulations-of-rotating-convection/7EDC0DA2FC127E7CCFC37D6CE825F621) show why transfer diagnostics and spatial location matter in rotating convection. Our first diagnostic is deliberately narrower than their full-flow analysis.

Five-minute bulk histories, hourly sections and vorticity movies continue. Movies normalize by central f₀ and show actual north–south coordinates; paired beta/f-plane movies use one fixed color scale and explicit labels. Monitor the evolving cap over the longer integration as well as energy, heat transport and spectral transfer.

## Validation and execution

The new driver lives in `scripts/run_regional.jl`; the earlier `scripts/run.jl` remains unchanged for reproducibility. GPU preflight checks compare the new periodic implementation to the legacy driver, compare BetaPlane(beta=0) to FPlane in the same channel, verify the imposed Coriolis gradient, check resting hydrostatic flow, and verify initialization amplitudes, wall-normal velocity, finite output, prescribed source powers and thermal inventory. Spectral checks cover analytic divergence-free modes, channel reflection, Parseval energy and conservation of self-transfer. The diagnostic and movie renderer are also exercised on actual preflight outputs.

Submit `sbatch slurm/beta-check.sbatch NEW_CHECK_ROOT` first. The production array `slurm/regional20.sbatch` must depend on successful completion of the preflight job. It runs at most one L40S at a time, with four CPUs, 24 GB host memory and a 12-hour limit per case. Previous throughput suggests order ten GPU-hours overall, with uncertainty from channel pressure transforms and finer grids; this is an estimate, not a measured campaign cost. Output is isolated in `data/regional20-v1`. No restart facility has been added, so a failed or timed-out case requires a new output root. A compact launch record records actual job IDs and verification status.

## Submitted campaign

Final preflight **66937526** completed successfully in 2 minutes 47 seconds on an L40S. Production array **66937528** was submitted with an `afterok` dependency and is now running, with the remaining cases queued at concurrency one. No production outcomes are claimed yet. The preceding preflight 66937496 also passed; the final run repeats validation against the finalized driver. See [exact parameter and launch record](../configs/regional20-parameters.json).

The offline transfer diagnostic passed analytic single-mode and nonzero-triad checks, with conservative nonlinear work and correct Parseval normalization. On actual short channel output, solenoidal projection removed about 8.4 × 10⁻⁶ of retained-mode energy; the high-wavenumber truncation fraction remained below 7 × 10⁻⁷, and total self-transfer residuals were of order 10⁻²³ m²/s³. These validate the preflight calculation, not the accuracy of a future turbulent cascade measurement. Production loss diagnostics must be inspected separately. The channel comparison video was rendered and visually checked for correct labels, signed y coordinates and fixed f₀ normalization.

## Interim finding after the three long runs

All three 20-day runs completed, but part of the stable cap becomes nearly neutral or weakly unstable. The late energy increase therefore cannot be interpreted as an inverse cascade alone. See the [interim comparison](../gallery/regional20-interim/README.md). The finer-grid comparisons remain underway.
