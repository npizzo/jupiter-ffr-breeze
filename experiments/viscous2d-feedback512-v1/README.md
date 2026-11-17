# Polynomial feedback plus linear instability: controlled LSV hypothesis

Three new horizontal 2D Oceananigans experiments test whether a rotation-dependent quadratic vorticity source changes large-scale organization when combined with the existing linear instability forcing. These are phenomenological experiments, not a derived convection closure or an established resonance mechanism. The user authorized this quadratic-feedback/cubic-saturation model on November 17; the cubic term is specifically part of this hypothesis, not the velocity drag rejected in earlier campaigns.

The equation is ∂tζ + u·∇ζ = ν∇²ζ + Lζ + P₀(bfζ² − cζ³). Here P₀ removes the spatial mean; L has Fourier multiplier a|k|² in band 12–16 and is zero elsewhere, including kx=0. All new cases use 512² cells, a doubly periodic 2π square, ν=0.0001, a=0.0002, c=0.0002, beta=0, and the same original initial velocity file (initial band12–16).

| Case | b | f | Initial velocities |
|---|---:|---:|---|
| Cubic-only control | 0 | +1 | Original |
| Biased plus | 0.002 | +1 | Original |
| Biased minus | 0.002 | −1 | x-reflection of original |
| Existing linear-only reference | 0 | +1 | Original; completed output reused; c=0 |

The minus-f case is a symmetry test, not an independent random realization. No selective viscosity, linear drag, velocity-cubic drag, or hyperviscosity is added. The polynomial term can supply kinetic energy; it is not assumed to be pure dissipation. The cubic contribution damps enstrophy for zero mean ζ. Maximum positive local fractional growth from bfζ−cζ² is 0.005 for |f|=1, compared with approximately 0.02 net linear growth near the middle of the forcing band. This is an exploratory coefficient choice, not a fit to Breeze.

Production targets dimensionless t=400, with saved vorticity/scalars every0.5 and checkpoints every10. A timestep ceiling0.0025 is reduced to target advective CFL0.3 and reaction step0.1. Stage guards stop on CFL≥0.5 or reaction step≥0.2; a proposed timestep below1e−5 stops for review. Each segment has a bounded runtime and retains a checkpoint; failures hold subsequent jobs. No repeated automatic resubmission is configured.

## Numerical implementation and checks

Ordinary viscosity uses the native ScalarDiffusivity. The polynomial is evaluated pointwise on native corner vorticity; its zero Fourier mode is removed. The unchanged linear spectral operator and polynomial curl are added, inverted with the native discrete Laplacian, and converted to momentum tendencies by the native rotated gradient. Three cached FFTs per stage suffice. This is a discrete pointwise closure, not a dealiased spectral-Galerkin polynomial model; high-wavenumber enstrophy is recorded as a resolution indicator.

Before release, Unity checks the independently evaluated source curl, mean-circulation preservation, actual momentum energy work versus curl inversion, the quadratic third-moment identity, cubic enstrophy work, reflected-f symmetry, recovery of the existing linear-only implementation when b=c=0, serialized restart, timestep refinement, exact output timing under adaptive stepping, finite-state/CFL and measured speed/allocation. A separate one-time production-driver smoke check verifies actual saved frames, scalar output and checkpoint creation before full production.

## What would count as progress?

The quadratic source has the exact semi-discrete moment identity

(d⟨ζ³⟩/dt)quadratic = 3bf[⟨ζ⁴⟩ − ⟨ζ²⟩²].

This identity is tested independently and recorded at every saved state. It is a source-term benchmark, not a claim that the model's full third-moment budget has been integrated or that discrete advection preserves all continuum Casimirs.

Scalar records include absolute kinetic energy in 0<|k|≤2, total energy/enstrophy, signed third moment/skewness, mean vorticity, high-wavenumber enstrophy, and direct energy input from the quadratic/cubic terms, both globally and into |k|≤2. An increase in large-scale energy must therefore be assessed alongside direct low-mode forcing; it cannot automatically be called an inverse cascade. The cubic-only case separates effects of saturation from the quadratic bias. The existing linear-only case supplies the original reference. A single seed and resolution do not establish robustness or convergence.

All computation, diagnostics and rendering remain on Unity, one CPU per job. On completion, the gallery will show a synchronized four-panel comparison (existing linear-only, cubic-only, biased +f, reflected −f) and an individual +f movie, with fixed shared symmetric asinh vorticity colors. These models have no depth layers or passive tracer.
