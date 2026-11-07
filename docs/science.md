# Scientific rationale and comparison plan

The target is sustained morphology: folded filaments, interacting lobes, embedded eddies, and intermittent convective outbreaks in cyclonic regions.

The objective is to establish whether Breeze can sustain convectively driven flow and reproduce FFR-like morphology. Density stratification, upper-layer stability, cyclonic shear, and Jovian moist convection are candidate ingredients to test individually and together.

## Sequence

1. Validate dry stratified convection on CPU and GPU; verify hydrostatic rest, heating sign, finite states, reference density and pressure, and meaningful vertical motion.
2. Run sufficient convective/eddy turnover times to assess organization. Four-hour pilots are not an equilibrated morphological benchmark.
3. Vary stable-layer structure and domain width; establish resolution and numerical-diffusion sensitivity.
4. Add cyclonic shear and a shear-only control. Track external energy input from the prescribed forcing.
5. Introduce water thermodynamics only after testing heavy-vapor buoyancy and condensation in a hydrogen–helium atmosphere. Earth microphysics is not assumed portable unchanged.
6. Compare against observed FFR morphology across multiple seeds and resolutions. Record runtime and memory to guide feasible parameter studies.

## Morphology benchmark

Visual comparison is the first screen. Across many snapshots and seeds measure filament width, curvature and folding; lobe and eddy sizes; cyclonic organization; and persistence. Repeat at refined resolution. Compare GPU-hours needed for spin-up and a sustained sample at a specified morphology quality, not just seconds per step. Preserve a distinction between vorticity/temperature structures and cloud imagery.

## Numerical decisions

Start with anelastic dynamics to retain density stratification while filtering sound. The first setup uses the traditional f-plane approximation (vertical component of rotation only); beta, the horizontal rotation component, and spherical curvature are future controlled additions. At 55°N, the omitted horizontal rotation component can matter for deep convective plumes, so this approximation needs a sensitivity test. Current Breeze PR #1103 provides compressible Yin–Yang support but excludes anelastic multi-region dynamics. Full-sphere coverage and overlap interpolation are unnecessary for the initial regional test.

Explicit diffusivities are numerical model parameters, not Jupiter's molecular diffusivities. The pilot uses Float64 as a validation reference. Benchmark Float32 only after comparing relevant fields and budgets. The pilot saves horizontal and vertical slices rather than large 3D archives. No automatic restart or campaign escalation is enabled.

## Reading

- [Breeze Yin–Yang PR #1103](https://github.com/NumericalEarth/Breeze.jl/pull/1103)
- [Fletcher et al., observed FFR structure](https://doi.org/10.1029/2025JE009315)
- [Hyder et al., FFR morphology in shallow water](https://arxiv.org/abs/2206.08475)
- [Sankar et al., shear and vapor transport](https://doi.org/10.3847/PSJ/adcc21)
- [Leconte et al., heavy-condensable inhibition](https://doi.org/10.1051/0004-6361/201629140)

## Initial nondimensional regime

Using the 80-km convecting depth, prescribed lower-layer buoyancy frequency squared of −10⁻⁷ s⁻², and explicit viscosity/diffusivity of 1,000 m²/s gives nominal Ra ≈ 4.1 × 10⁶, Ek = ν/(fH²) ≈ 5.4 × 10⁻⁴, Pr = 1, and convective Rossby number √(−N²)/f ≈ 1.1. These are local dry-layer estimates in a stratified atmosphere. The first pilot is therefore not in the asymptotic rapidly rotating limit. Varying superadiabaticity and diffusivity will test stronger rotational constraint. Numerical advection also contributes effective diffusion.
