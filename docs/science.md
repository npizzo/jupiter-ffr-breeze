# Scientific rationale and comparison plan

The target is sustained morphology: folded filaments, interacting lobes, embedded eddies, and intermittent convective outbreaks in cyclonic regions. GOFLOW training is outside this project's scope.

Rapidly rotating Rayleigh–Bénard convection is a substantive baseline: it already supports upscale energy transfer and organized large vortices. The alternative hypothesis is that density stratification, upper-layer stability, background shear, and eventually Jovian moist convection improve FFR morphology at practical cost.

## Sequence

1. Validate dry stratified convection on CPU and GPU; verify hydrostatic rest, heating sign, finite states, reference density and pressure, and meaningful vertical motion.
2. Run sufficient convective/eddy turnover times to assess organization. Four-hour pilots are not an equilibrated morphological benchmark.
3. Vary stable-layer structure and domain width; establish resolution and numerical-diffusion sensitivity.
4. Add cyclonic shear and a shear-only control. Track external energy input from the prescribed forcing.
5. Introduce water thermodynamics only after testing heavy-vapor buoyancy and condensation in a hydrogen–helium atmosphere. Earth microphysics is not assumed portable unchanged.
6. Compare with the existing RBC implementation when available. Match useful nondimensional measures and display the same diagnostic or passive tracer at the same resolution; do not treat different Rayleigh-number definitions as equivalent.

## Morphology benchmark

Visual comparison is the first screen. Across many snapshots and seeds measure filament width, curvature and folding; lobe and eddy sizes; cyclonic organization; and persistence. Repeat at refined resolution. Compare GPU-hours needed for spin-up and a sustained sample at a specified morphology quality, not just seconds per step. Preserve a distinction between vorticity/temperature structures and cloud imagery.

## Numerical decisions

Start with anelastic dynamics to retain density stratification while filtering sound. The first setup is f-plane; beta and spherical curvature are future controlled additions. Current Breeze PR #1103 provides compressible Yin–Yang support but excludes anelastic multi-region dynamics. Full-sphere coverage and overlap interpolation are unnecessary for the initial regional test.

Explicit diffusivities are numerical model parameters, not Jupiter's molecular diffusivities. The pilot uses Float64 as a validation reference. Benchmark Float32 only after comparing relevant fields and budgets. The pilot saves horizontal and vertical slices rather than large 3D archives. No automatic restart or campaign escalation is enabled.

## Reading

- [Breeze Yin–Yang PR #1103](https://github.com/NumericalEarth/Breeze.jl/pull/1103)
- [Rubio et al., upscale energy transfer](https://doi.org/10.1103/PhysRevLett.112.144501)
- [Guervilly et al., large-scale vortices in rotating RBC](https://arxiv.org/abs/1403.7442)
- [Fletcher et al., observed FFR structure](https://doi.org/10.1029/2025JE009315)
- [Hyder et al., FFR morphology in shallow water](https://arxiv.org/abs/2206.08475)
- [Sankar et al., shear and vapor transport](https://doi.org/10.3847/PSJ/adcc21)
- [Leconte et al., heavy-condensable inhibition](https://doi.org/10.1051/0004-6361/201629140)
- [Julien et al., efficient rapidly rotating convection](https://doi.org/10.1016/j.jcp.2025.114274)
