# Small instability-forced 2D Oceananigans experiment

This exploratory experiment runs on Unity, independently of the Breeze and
Ferrel-cell campaigns. It tests spontaneous sign selection, not a preferred
sign relative to rotation. There is no beta, heating, cap, or vertical direction.

The forcing is the band-limited negative-viscosity operator from the supplied
`beta_channel.f90`: Aζ = a P_band (-∇²)ζ. The previously validated channel
projection matches the original Fourier–sine sums to roundoff. Here the box is
doubly periodic to remove channel-wall effects. Modes with kx=0 remain excluded
from forcing, but are allowed in the evolving flow.

## Bounded pilot

| Parameter | Value |
|---|---|
| Grid and box | 128 × 128; 2π × 2π |
| Forcing band | 6 ≤ \|k\| ≤ 8 |
| Feedback amplitude a | 0.009375 |
| Largest injection growth rate | 0.6 |
| Constant Laplacian viscosity | 0.001 |
| Cubic velocity drag | −\|u\|²u, with staggered interpolation |
| Initial RMS speed | 0.02; random seed 1701, modes 4–10 |
| Coriolis parameter | f=1; no beta |
| Time step | 0.01, dimensionless |
| Target | 120 injection growth times = 200 time units |
| Compute ceiling | One CPU; 50-minute job, 40-minute pilot limit; no extensions |

The band is moved down from 18–20 to give 16–21 grid cells per forced wavelength
at modest cost. Adjusting a preserves the original maximum injection growth
rate. The fixed viscosity and cubic drag differ from the supplied Fortran
RMS-dependent viscosity and from the paper's hyperviscosity. This is not a
quantitative paper reproduction. No mapping to Earth/Jupiter days is assigned.

The spatial mean of the cubic-drag acceleration is removed, holding the
zero-mean-velocity subspace. Without this projection, nonuniform drag can
generate a uniform current that undergoes f-dependent inertial oscillations;
that is not vortex-sign selection. This projection does not alter drag work
when velocity has zero mean. The damping-enabled f-flip validation detected
this issue before any turbulence run, and the check was retained unchanged.

`check.jl` must pass before `run.jl` starts. Checks cover analytic cubic drag
and viscous decay, native RK vorticity increments, third-moment accounting,
reflection symmetry, and f-independence with both damping terms enabled.
Existing operator validation also checks band growth and the original Fortran
projection. All scientific checks run on Unity.

The budget uses all three native Wray RK stages. Each stage's exact identity
q_new³−q_old³=(q_new−q_old)(q_new²+q_new q_old+q_old²) apportions the increment
among the actual operators. The pressure-projected vorticity increment is
checked independently against their sum. Advection/Coriolis remains an explicit
term: incompressible 2D advection should conserve the domain third moment in
the continuum, but its discretization need not. A numerical contribution above
10% of the summed absolute physical step contributions blocks a mechanism
claim. Algebraic closure alone does not certify physical attribution.

Outputs include vorticity frames every 1.2 growth times, skewness, energy,
enstrophy, low-wavenumber energy fraction, operator budgets, a restart state,
configuration/source hashes, and a final review. A partial wall-time-limited
run is explicitly labeled. One realization is not evidence of ensemble sign
preference, and low-k energy growth is not a measured spectral flux.

## A possible later rotation-sensitive extension

For a parallel inviscid shear, f(f+ζ)<0 is an inertial-instability criterion.
For an axisymmetric barotropic vortex the relevant discriminant instead is
Φ=(f+2V/r)(f+ζ), with signed azimuthal speed V. Stratification, viscosity,
non-axisymmetry and finite geometry can modify applicability and growth.
Negative vorticity by itself does not establish instability.

These instabilities involve motions excluded by strict horizontal 2D dynamics.
A reduced closure could represent energy lost to those motions using a
nonnegative additional viscosity ν_inst=C ℓ² sqrt(max(0,−Φ)), applied as the
divergence of a viscous stress. That form dissipates kinetic energy in a
periodic incompressible domain and preserves net circulation. The expression
is a proposed closure, not a derived or validated model. Identifying V/r in
non-axisymmetric turbulence and calibrating C and ℓ need independent evidence.
The parallel-shear proxy f(f+ζ) must not be substituted silently for Φ.

No such correction is enabled in this pilot. First establish the unbiased
baseline. Then test a proposed correction against isolated-vortex 3D growth
and saturation, reflection/f-reversal covariance, and energy loss before
interpreting any selected sign. An empirically imposed anticyclone sink would
not explain the Breeze result by itself.

References: van Kan et al. (2024), https://arxiv.org/html/2308.08789v2;
Billant & Gallaire (2005), https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/generalized-rayleigh-criterion-for-nonaxisymmetric-centrifugal-instabilities/EBD5A1C0273A1D8616C5F1FC75085DCF.
