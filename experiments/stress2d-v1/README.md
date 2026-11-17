# Projected nonlinear stress: a minimal sign-selection hypothesis

Authorized November 17, 2026. This is a new isolated Oceananigans campaign, separate from the existing Breeze and beta-plane experiments. All computation, initial-field preparation, diagnostics, rendering and media validation run on Unity.

## Equation and assumptions

For horizontal incompressible flow U, define Z = curl U and S = (grad U + grad Uᵀ)/2. P is the orthogonal Fourier projection onto |k|≤16 in a periodic 2π square, including zonal modes. The model applies

    Ubar = P U
    Rdev = P[-2 νT(Zbar) Sbar]
    momentum forcing = -div Rdev
    νT(Zbar) = ν0 + ν1 Zbar/f + ν2 (Zbar/f)^2.

Ordinary molecular viscosity ν=0.0001 acts on every resolved mode. There is no additional van Kan instability forcing, linear drag, velocity-cubic drag, hyperviscosity, beta or temperature equation. The constant Coriolis term is retained in Oceananigans and separately checked to have no effect on the periodic horizontal divergence-free dynamics. Sign dependence comes from the stress coefficients, not from direct 2D Coriolis acceleration.

This is the strain-aligned member of a larger constitutive family. The scale separation, short correlation time and smooth weak response of the unresolved convection are assumptions, not demonstrated facts. No convective eddies are resolved in this 2D experiment. Coefficients are exploratory, not inferred from Breeze. The test can show what this closure does, not derive it from rotating convection or prove a resonance mechanism.

## Cases and initial conditions

| Case | ν0 | ν1 | ν2 | f | Initial state |
|---|---:|---:|---:|---:|---|
| Ordinary decay | 0 | 0 | 0 | +100 | Original pattern |
| Symmetric stress | -0.0002 | 0 | 0.01 | +100 | Original pattern |
| Biased +f | -0.0002 | -0.001 | 0.01 | +100 | Original pattern |
| Biased -f | -0.0002 | -0.001 | 0.01 | -100 | x-reflected original |

All cases use 256² cells and target dimensionless t=400. The existing original 512² initial condition has its vorticity in a low Fourier band; Unity verifies the absence of power above the new Nyquist limit, samples the same corner-vorticity pattern on 256², and reconstructs divergence-free velocities with the native discrete Laplacian. All positive-f cases have identical velocities. The negative-f case is a reflection check, not an independent seed. Time and f are nondimensional, not Earth days or Jupiter's dimensional rotation rate.

The initial vorticity RMS is approximately 2.758, giving RMS Z/f≈0.0276. The coefficients put the net local viscosity zeros at Z/f≈+0.162 and −0.062 in the biased case, and ±0.1 in the symmetric case. These zeros describe a stress response, not predicted vortex equilibria. Maximum linear growth is approximately 0.0256 near |k|=16. Low modes also receive direct stress forcing: growth of a large vortex is not, by itself, evidence of an inverse cascade.

## Conservative implementation

The packed complex FFTs apply P independently to the real and imaginary velocity components. Deviatoric normal stress is stored at cell centers and shear stress at corners. Normal strain is (∂xU−∂yV)/2, shear strain is (∂yU+∂xV)/2, and center vorticity is the arithmetic mean of the surrounding four corner values. Ryy=-Rxx. The three polynomial stress contributions are then projected back to the same band before taking their native discrete divergence. No division by a small vorticity and no reactive mean-removal source are used.

The input support is |k|≤16 and the highest polynomial degree is three, so stress products have support no larger than |k|≤48, below the 256² Nyquist limit. Thus the nonlinear stress evaluation is not aliased. The projection is self-adjoint and the staggered divergence/strain operators are adjoints. At each saved state the independently computed momentum work is checked against

    work = -4 mean[νT(Zcenter) Sxx² + νT(Zcorner) Sxy²].

The cubic contribution is nonpositive. The quadratic contribution has no universal sign on an arbitrary field. Mean momentum and circulation are conserved by the stress. Spatial derivatives are retained: for a one-dimensional shear the source reduces to Dyy[νT(Z) Z], rather than a local polynomial reaction.

## Validation and bounded cost

Before production, a single Unity CPU checks independent force evaluation, circulation, momentum work, cubic damping, the constant-coefficient discrete Laplacian limit, the exact shear reduction, spectral cutoff, reflected-f symmetry, cancellation of constant-f Coriolis, serialized restart, timestep refinement, output timing and performance. A separate production-driver smoke check verifies actual files and checkpoints. Production is released only after both pass and the source hashes match the gate.

The timestep ceiling is 0.005 with advective CFL target 0.3 and an explicit stress-rate guard. A timestep below 0.00025 or filtered max |Z/f|≥0.5 stops for review. There is no automatic retry or parameter search. Each production case has one CPU, 4 GB memory, a one-hour Slurm ceiling and a 55-minute internal wall limit. Partial checkpoints are retained. All four cases together have a four CPU-hour production ceiling.

## Output and interpretation

Vorticity frames and scalar diagnostics are saved every model time unit; checkpoints every ten units. Records include total and |k|≤2 energy, enstrophy, mean vorticity, signed third moment/skewness, high-wavenumber content, filtered RMS/max |Z/f|, area fraction above 0.2, strain/f and the energy and third-moment contributions of each stress term. The latter are source contributions, not an integrated closure of the full time-discrete budget. An algebraic energy identity does not validate the physical stress assumption.

The gallery will show synchronized vorticity evolution for all four cases, and an individual biased +f movie, with a shared fixed asinh color scale. These are genuinely horizontal 2D simulations with one plane and no passive tracer. The weak-response interpretation must be judged against the recorded Rossby and strain measures; a finite-state run alone does not establish validity or an LSV mechanism. One seed and one resolution do not establish convergence or statistical robustness.
