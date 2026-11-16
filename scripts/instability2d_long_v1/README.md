# Longer matched 2D forcing-on/off experiment

This is a continuation of the completed 128² Oceananigans instability-forcing
pilot, not a new Breeze convection run. The user requested a longer forced and
unforced comparison to determine whether inverse energy transfer is operative.
All integration, analysis, plots, movie rendering and validation run on Unity.

## Controlled comparison

Both branches load exactly the same completed pilot checkpoint at sigma*t=120
(t=200, sigma=0.6). The forced branch retains band-limited negative viscosity;
the unforced branch sets only its amplitude to zero. Both continue to
sigma*t=1200 (t=2000), adding 1080 reference growth times. The reference sigma
also labels the unforced clock; it is not an instability rate in that branch.
This is a mature-state forcing switch, not two runs independently initialized
from small random noise. The initial RMS speed is about 0.361, versus 0.02 at
the start of the original pilot. Checkpoint SHA256 is recorded in both cases.

Everything else is unchanged: 128², 2pi periodic box, f=1, no beta, native
centered advection, dt=0.01, forcing band 6–8 excluding kx=0, coefficient
0.009375 when enabled, Laplacian viscosity 0.001, cubic velocity drag coefficient
1, and removal of the mean drag acceleration. There is no rotation-sensitive
closure, tracer or vertical dimension. A paired single-plane vorticity movie
is the honest visualization; three artificial depths must not be created.

## What establishes inverse transfer

The Fourier diagnostic uses the native discrete Laplacian eigenvalue lambda:

    E_k = |qhat_k|² / (2 N² lambda_k)
    T_j,k = Re(conj(qhat_k) * tendencyhat_j,k) / (N² lambda_k)
    Pi(K) = - sum_{|k| <= K} T_advection,k

Here N=nx*ny is the number of horizontal cells, q is native corner vorticity,
and Fourier transforms are unnormalised. The zero mode is omitted; the
zero-mean, divergence-free velocity has kinetic energy exactly equal to this
sum. Native advection includes Coriolis; its curl vanishes in the incompressible
periodic constant-f setting to solver precision. The energy conservation and
Parseval checks test that assumption directly.

Shell K means K-1 < |k| <= K, with explicit shell zero. Negative Pi means net
transfer toward larger scales across K. We report its time history and window
means for K=1–5, below the forcing band, alongside energy spectra, integral
length and absolute energy at |k|<=2. A rising low-k fraction alone does not
establish inverse transfer: preferential small-scale decay can produce it.
The unforced case has a finite decaying energy reservoir, not steady injection.
With forcing at 6–8 there is only a short interval below injection; a broad
inertial range or a -5/3 law is not demonstrated by this small box.

Separate injection, viscous and drag spectral work is saved every 100 steps
(0.6 reference growth times). Direct forcing below k=6 must be zero. Native
Runge–Kutta energy increments are also apportioned at every stage using
Delta E = -mean[(psi_new+psi_old) Delta q]/2, and independently checked against
kinetic energy storage. The previous exact vorticity/third-moment accounting
and its 10% numerical-advection gate remain in place; that gate is about sign
attribution, not a substitute for spectral-transfer evidence. Algebraic budget
closure alone does not establish resolution convergence or a physical mechanism.

## Validation and compute bounds

Before production, rerun the inherited ten analytic/symmetry checks, then test
both forcing settings on the actual checkpoint. Verify native energy closure,
Parseval energy/work equivalence, no low-k forcing leakage, and an independent
manufactured upscale exchange whose flux sign and normalization are known.
Both source checkpoint hashes must match. Production refuses to overwrite data.

Two single-CPU array tasks, each capped at four wall hours with an internal
3.5-hour stop. Partial runs retain checkpoints and return failure so a dependent
complete-run renderer cannot publish them as complete. No automatic extension,
parameter sweep or resolution change. A failure of the original CFL<0.5 guard
halts the affected case for inspection rather than silently altering dt.

Vorticity is saved every 200 steps (1.2 reference growth times), so each complete
branch has 901 frames. Checkpoints are written every 2000 steps and at exit.
The dependent analysis produces a synchronized two-panel movie with fixed
shared physical color limits, full decoding, spectral plots and numerical
summaries. Completed assets still require review before gallery publication.

## Interpretation

Longer integration tests whether the pilot was transient; it does not guarantee
condensation. Van Kan et al. (2024), especially section 6.2, show that shielded
vortices and nonlinear damping can suppress inverse transfer:
https://arxiv.org/html/2308.08789v2. The present Laplacian-viscosity pilot is not
a quantitative reproduction of their hyperviscous experiments. Damping is kept
fixed in this first pair so the effect of switching forcing is interpretable.
No conclusion about a sign tied to f can follow from this strictly horizontal
2D model without an independently justified three-dimensional closure.
