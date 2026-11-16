# A minimal 2D instability-forcing experiment

This is an isolated Oceananigans CPU prototype, not a change to the Breeze
convection or beta-plane campaigns. Scientific validation runs on Unity.
It ports the **linear instability operator**, not the entire Fortran solver.

## Question this model can answer

Can sign-symmetric, flow-dependent injection produce persistent one-signed
vortex cores through spontaneous symmetry breaking? That is distinct from
asking why rotation prefers cyclones over anticyclones in convection.

For a horizontal, incompressible, constant-density 2D flow with constant f,

    ∂t ζ + u·∇ζ = A ζ + ν∇²ζ + Dζ,
    A = a P_band (-∇²).

The curl of the Coriolis acceleration is -f ∇·u = 0. There is no f-dependent
vorticity production and no vortex stretching. Changing f may change pressure;
it must not change the vorticity dynamics of the zero-mean periodic test.
The mean-flow caveat matters: on a periodic domain a spatially uniform velocity
can undergo inertial oscillations. All current tests use zero mean velocity.

A spatial reflection flips vorticity while preserving these equations and
the band forcing. A plain sign flip of velocity at fixed coordinates is not
the correct forward-time nonlinear symmetry test. A reflection-symmetric
ensemble can therefore have E[S] = 0 but E[|S|] > 0, where S is vorticity
skewness. One realization of either sign is not evidence of a preferred
cyclonic sign. The opposite-sign shield or background must also be measured:
a periodic velocity field has zero domain-integrated vorticity.

## What was read from beta_channel.f90

The supplied file's `add_instability` (lines 409–462) projects ζ onto
cos(kx)sin(ly) and sin(kx)sin(ly), multiplying each admitted coefficient by
`gamnu * (k*k + l*l)`. Its Euler kick is dt times this tendency.
The domain is Lx=2π, Ly=π, with periodic x and zero ψ and ζ at the y walls.
The current defaults are gamnu=0.0015, instability shell 18–20, and exclusion
of kx=0 modes. The maximum linear growth rate is 0.6 in the code's time units.
There is no conversion from these units to Earth days without a scale choice.

The separate stochastic forcing has amplitude 0.5 and shell 4–6; initial
noise is in 14–16 with initial RMS speed 0.15. There is a wall sponge.
The scalar viscosity is recomputed as ζ_rms min(dx,dy)^2. Thus grid refinement
with an unchanged `viscosity_factor` changes the governing regularization.
These ingredients are NOT silently treated as the 2024 paper's setup.

Source SHA256: `38fc879875e73a8929381832bc1b534c849a50c569d042997f47d32d71fd0468`.
The user's original Fortran file is neither modified nor copied into this repo.

## Implementation

`core.jl` supports a periodic horizontal box or a free-slip channel with a flat
vertical direction. Channel spectral projection uses FFTs of an odd extension
in y, equivalent to the Fortran Fourier–sine projection. A separate direct
trigonometric sum in `validate.jl` checks this equivalence, including zonal
mode inclusion/exclusion.

The multiplier uses physical continuous k², exactly as the Fortran operator.
The conversion from curl forcing to face-centered momentum forcing uses the
**discrete** Laplacian and compatible C-grid differences, so its discrete curl
recovers the requested forcing and its discrete divergence vanishes. Feedback
is evaluated at every RK tendency evaluation, not only once per timestep.
This is the same continuous-time forcing with improved temporal integration;
it is not a bitwise replication of the Fortran's split Euler kick or Arakawa
advection scheme.

The initial prototype is CPU/Float64 only. No GPU performance claim is made.
There is no stochastic forcing or wall sponge in this prototype. Both are
deliberately absent in the minimal pure-instability experiment. Options expose
the Fortran's global RMS-dependent viscosity and a velocity-cubic drag
`-alpha |u|² u`. Cubic drag must be applied to velocity, not replaced with a
different law `-alpha ζ³`. These damping options require their own validation
and resolved-regime checks before scientific production.

Example construction (dimensionless; not a production prescription):

```julia
include("core.jl")
using .Instability2D, Oceananigans
op = BandOperator(256, 256; kmin=18., kmax=20., allow_zonal=true)
model = build_model(op; f=1., nu=1e-4)
# Seed a zero-mean, divergence-free field with velocity_from_curl and seed!.
sim = Simulation(model; Δt=0.001, stop_time=1., verbose=false)
feedback = attach!(sim, op; amplitude=0.0015, cubic_drag=1e-4)
# run!(sim) only after adding the required diagnostics and cost limits.
```

## Why the third-moment experiment is cleaner

In a periodic continuum model, incompressible advection gives exactly zero
contribution to d<ζ³>/dt. For ordinary viscosity,

    d<ζ³>/dt = 3<ζ² Aζ> - 6ν<ζ |∇ζ|²> + 3<ζ² Dζ>.

There is no depth-average closure, cap, sponge, or reconstructed baroclinic
stress in this equation. The injection and damping are prescribed, accessible
operators. Nevertheless, finite-volume advection need not conserve the third
moment. Its numerical contribution must be measured, not hidden in a residual.
Required production outputs are the native moment increment, the forcing and
damping contributions, the advection contribution, and the remaining closure
error, sampled every timestep with RK-consistent accounting.

Even linear instability forcing does not automatically amplify normalized
skewness: if all occupied modes have the same growth rate, amplitudes grow
uniformly and normalized skewness is unchanged. A sign-selection mechanism
must involve mode coupling, damping, and/or competition between coherent
vortices. Correlating injection with ζ² alone does not establish that mechanism.

## Relation to the literature

van Kan, Favier, Julien & Knobloch (2024), *From a vortex gas to a vortex crystal
in instability-driven two-dimensional turbulence*, JFM 984:
https://arxiv.org/html/2308.08789v2

Their standard configuration is doubly periodic, with shell 33–40, ν*=0.002,
eighth-order hyperviscosity ν4=1e-14, cubic velocity drag coefficient 1e-4,
and a mixing parameter between instability and random forcing. Their drag
coefficient is called β in the paper; it is NOT a planetary beta-plane term.
The pure-instability case has no random forcing. The paper finds spontaneously
chosen same-sign shielded cores, often a vortex gas/crystal rather than one
domain-scale vortex. Its long transient matters: Figure 1 shows an already
symmetry-broken state at σt=378 and a much later dense state at σt=9312.
These are reported examples, not predicted transition times for this prototype.

The 2022 precursor is https://doi.org/10.1017/jfm.2022.935.
For an exact paper reproduction, a dealiased spectral method with appropriate
stiff-linear integration is a better fit than blindly porting eighth-order
diffusion to an explicit finite-volume RK scheme. This implementation is a
minimal test of the forcing mechanism, not that reproduction.

## Bounded next scientific experiment

1. Pass operator, modal growth, divergence, reflection, and f-independence
   checks on small Unity CPU grids. No physical inference from these checks.
2. Add RK-consistent third-moment accounting and validate the damping options.
   Reject the configuration if numerical advection contributes more than 10%
   of the mechanism-sized third-moment tendency, using interval RMS or a
   non-cancelling norm rather than a nearly zero net six-hour integral.
3. Use a periodic 2π square, forcing shell 18–20, no stochastic forcing, zero
   mean flow, and a small unbiased broadband initial field. Start with fixed
   viscosity plus cubic drag, then compare the Fortran RMS-viscosity choice
   only as a separately labeled closure-sensitivity test. No beta.
4. Estimate cost with one short 256² pilot before scheduling long runs.
   Set a fixed compute ceiling and checkpoint. Use instability growth times,
   not arbitrary Earth days. Lack of selection in a short pilot is inconclusive.
5. Only if the budget and cost checks pass, use several independent seeds and
   reflected partners. Reflected partners check symmetry, not independent
   statistical samples. Confirm a candidate regime at 512² with identical
   physical damping and band parameters, and at half timestep.

Measure signed/absolute skewness, extreme-core number and circulation by sign,
shield profiles, core lifetimes, low-wavenumber energy fraction, injection into
each sign population, and spectra. An inverse-cascade claim additionally needs
spectral energy flux. A one-signed vortex gas with little large-scale energy is
not the same outcome as the mature Breeze single-vortex state.

If sign selection is robust here, it demonstrates that sign-symmetric feedback
is sufficient for one-signed populations in this model. It does not establish
that Breeze selects sign by the same mechanism. A reproducible preference
relative to f in Breeze still needs a rotation-sensitive 3D/convection closure.
No f-sensitive term should be inserted merely to manufacture the desired bias.
