# Anticyclone-selective dissipation in decaying 2D turbulence

This experiment tests whether explicitly suppressing anticyclonic vorticity promotes a dominant cyclonic vortex. It is a phenomenological mechanism test, not a derivation of a convection or centrifugal-instability closure.

The model is horizontal, incompressible Oceananigans flow in a doubly periodic square of side 2π, with f = +1, no beta, no forcing and no drag. Positive vertical vorticity is counterclockwise and cyclonic. It uses the same saved initial velocities as the existing 512² unforced control, whose completed output is reused. The initial turbulence occupies wavenumbers 12–16. Integration runs from dimensionless t = 0 to 400, with dt = 0.0025 and ordinary viscosity ν = 0.0001.

## Added operator

Let s = sign(f) and q_A = s min(s ζ, 0). The vorticity equation is

∂t ζ + u·∇ζ = ν ∇²ζ + ν_A ∇²q_A,

with ν_A = 0.0004. Away from sign interfaces, the total viscosity is 5ν in anticyclonic regions and ν in cyclonic regions. The sharp selector is continuous and monotone; its derivative changes at zero. No smoothing, instability forcing, random stirring, drag or hyperviscosity is added.

A periodic velocity field must have zero domain-mean curl. Simply damping negative ζ pointwise would violate this condition. The flux-form operator above instead preserves circulation and dissipates kinetic energy and enstrophy:

- dE/dt from the selective term = −ν_A ⟨q_A²⟩ ≤ 0.
- dZ/dt from the selective term = −ν_A ⟨H(−sζ)|∇ζ|²⟩ ≤ 0, interpreted weakly at sign interfaces.
- ⟨∇²q_A⟩ = 0.

These identities follow by periodic integration by parts. The actual implementation uses native C-grid differences: F_u = −ν_A D_y⁺q_A and F_v = ν_A D_x⁺q_A. Its discrete curl is the five-point Laplacian of ν_A q_A. Periodic summation by parts gives the corresponding discrete energy and enstrophy inequalities. Fluxes at sign interfaces can also change neighboring cyclonic cells: this does not remove anticyclonic circulation while leaving all cyclonic circulation untouched. A dominant cyclone still requires compensating anticyclonic circulation elsewhere.

## Release checks and interpretation

A short Unity validation must pass before production: identical initial field, independent curl stencil, zero mean vorticity tendency, energy and enstrophy work identities, reflection combined with f reversal, serialized restart, halved-timestep comparison, finite-state/CFL checks, and measured speed/allocation limits. The selective callback uses local loops and no extra FFTs. One CPU is requested per job; stages stop on failure rather than releasing dependent jobs.

Saved scalar diagnostics record energy, enstrophy, signed enstrophy, skewness, mean vorticity and the selective dissipation rates every 0.5 model time. These are not a claim of a closed total discrete vorticity budget or proof of inverse energy transfer. Movies will compare this case with the existing unforced decay using a fixed shared symmetric asinh vorticity scale, width 0.1, with no clipping. This genuinely 2D model has one plane and no tracer.

A change in vortex sign or morphology would show what this imposed asymmetric dissipation can do. It would not establish that convection produces this operator. Stronger total dissipation is also introduced; a symmetric dissipation-matched control would be a further discriminator if the first test shows an interesting response. No such additional production run is included in this first test.
