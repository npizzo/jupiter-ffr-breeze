# Finite-window convective response to a prescribed barotropic vortex

Derived November 18, 2026. This document specifies a proposed Boussinesq response calculation. It is a reduction of the perturbation equations, not a replacement of the original anelastic Breeze run and not a saturated convection model. The reference cap and boundary conditions must be recorded with the chosen input profile.

## Parent equations and what the response means

Let U(x,y) be a prescribed, horizontally divergence-free, depth-independent velocity. Let B(z) be a frozen horizontally averaged reference buoyancy profile and N²(z)=dB/dz. Linear perturbations obey

    ∂t u + U·∇h u + u·∇h U + f zhat×u = −∇h p + ν Δ u,
    ∂t w + U·∇h w                    = −∂z p + b + ν Δ w,
    ∂t b + U·∇h b + N²(z) w          = ∇·(κ(z)∇b) − λ(z)b,
    ∇h·u + ∂z w = 0.

Here u is the perturbation horizontal velocity. The prescribed U has no vertical velocity or vertical shear. Viscosity and thermal diffusion are ordinary second-order operators. The optional λ is the linearized thermal cap relaxation, not momentum drag; include it only if explicitly part of the declared configuration. No beta, selective vortex damping, polynomial reaction or stochastic stirring is implied.

Impermeable, stress-free momentum boundaries are w=0 and ∂z u=0. The old specified heat-flux boundary implies ∂z b=0 for perturbations when the imposed flux is spatially uniform and unchanged. Setting b=0 instead would change the thermal boundary condition.

Because U and B are held fixed, this computes finite-time growth and a *virtual* tendency of the barotropic vortex. Actual barotropic energy does not increase in this calculation. Holding the backgrounds fixed corresponds to external maintenance that must not be confused with self-sustained turbulence. Stop before perturbations cease to be small; no saturation or stationary-covariance assumption is justified.

Removing the imposed heat flux while retaining the same frozen B leaves these perturbation equations unchanged. An informative control therefore needs a separately declared reference gradient, such as replacing only negative lower-layer N² by zero while retaining the stable cap, or a measured unheated-run profile. Call the former a **convective-gradient-removed control**, not a literal switch-off-heating experiment.

## A vertical Galerkin basis that accommodates the cap

Full-height sine modes can spuriously stabilize the strongly capped reference by forcing vertical motions through the stable upper third. A localized basis can be selected using an independent one-dimensional stability calculation. This is a representation choice; it does not alter the input N² profile.

Choose linearly independent real vertical-velocity functions g_n(z), n=1,…,M, satisfying g_n=0 and g_n''=0 at both physical boundaries. Choose a horizontal basis h_a(z), a=1,…,Nh, that spans **both** their derivatives g_n' and the toroidal horizontal-velocity profiles of the selected full-reference eigenmodes. All h_a have zero depth mean and satisfy h_a'=0 at the walls. Retaining the toroidal shapes prevents an avoidable loss of the selected eigenmode when its toroidal component differs from g_n'. One retained convective mode can therefore require two horizontal shapes.

Write

    u(x,y,z)=Σ_a u_a(x,y) h_a(z),
    w(x,y,z)=Σ_n w_n(x,y) g_n(z),
    b(x,y,z)=Σ_l b_l(x,y) β_l(z).

The β_l are orthonormal buoyancy basis functions satisfying homogeneous Neumann conditions; cosine functions including the constant mode are a simple choice. Their number need not equal the number of velocity modes. Define the exact derivative map J by

    g_n' = Σ_a h_a J_an,
    H J = <h_a,g_n'>_z.

If the derivative reconstruction does not close within numerical tolerance, reject the basis. The coefficient constraint is

    ∇h·u + J w = 0.

Use the depth-average inner product <a,b>_z=Hdomain⁻¹∫a b dz and define

    H_ab=<h_a,h_b>_z,       W_mn=<g_m,g_n>_z,
    C_ml=<g_m,β_l>_z,       B_ln=<β_l,N² g_n>_z,
    Vh_ab=<h_a',h_b'>_z,    Vw_mn=<g_m',g_n'>_z,
    Kx_lr=<β_l,κ β_r>_z,    Kz_lr=<β_l',κ β_r'>_z,
    Γ_lr=<β_l,λ β_r>_z.

H and W are positive-definite mass matrices when their respective bases are independent. The special choice h_n=g_n'/q_n is sufficient when toroidal profiles lie in this same span. For normalized full-height sine/cosine modes, H=W=I and J=diag(nπ/Hdomain). For unnormalized sine/cosine modes, the familiar factor one-half appears in H and W; do not mix the conventions.

The unprojected coefficient tendencies are

    F_u = −U·∇h u − u·∇h U − f zhat×u
          + ν(Δh u − H⁻¹ Vh u),
    F_w = −U·∇h w + W⁻¹ C b
          + ν(Δh w − W⁻¹ Vw w),
    F_b = −U·∇h b − B w + Kx Δh b − Kz b − Γ b.

These formulas presume a depth-independent ν. A variable ν requires the analogous weighted derivative matrices. No vertical pressure approximation is needed beyond the chosen velocity subspace: enforce continuity by an orthogonal projection in the kinetic-energy metric.

## Pressure projection

For each horizontal Fourier wavevector k=(kx,ky), stack the coefficient vectors v=(u_x,u_y,w). Let

    M = diag(H,H,W),
    D = (i kx I_Nh, i ky I_Nh, J).

The divergence-free subspace is Dv=0. Its energy-orthogonal projector is

    P = I − M⁻¹ D† (D M⁻¹ D†)⁻¹ D.

Advance ∂t v=P F_v. For nonzero horizontal k the small Nh×Nh solve can be precomputed for each horizontal |k|². At k=0, J has full column rank and continuity requires w=0 while horizontal coefficients remain unconstrained. Set the vertical tendency to zero and retain the horizontal tendency there; this avoids inverting the singular pressure-multiplier matrix when Nh>M. P satisfies DP=0, P²=P and P†M=MP. Thus pressure does no work on a constrained state, even when H and W are not diagonal. Project the initial perturbations and every Runge–Kutta stage. Ordinary Laplacian diffusion is projected along with the other tendencies for a general basis.

For orthonormal Fourier momentum modes this reduces to P=I−g g†/(|k|²+q_n²), g=(i kx,i ky,−q_n), separately for each vertical mode. This limiting case supplies an independent algebra check on the general implementation.

An additional exact limiting test uses one normalized full-height sine velocity mode, uniform N², zero viscosity/diffusion, no background vortex, and any finite Neumann thermal basis. Let c be its buoyancy-overlap row C. The active wave branch obeys

    ω² = [f² q² + N² ||c||² |k|²] / [q² + |k|²].

For ω²<0 the growing branch has rate sqrt(−ω²). The finite-rank factor ||c||² is essential: it approaches one as the thermal basis is enlarged, but a finite cosine basis cannot represent a sine exactly. This benchmark tests the implemented coupling against an independently derived scalar formula.

## Exact instantaneous exchange identity

The retained perturbation kinetic energy is

    K' = 1/2 <u_xᵀH u_x + u_yᵀH u_y + wᵀWw>_h.

Its depth-averaged horizontal stress is

    R_ij = Σ_ab H_ab u_i,a u_j,b.

The corresponding virtual barotropic power is

    P_BT = <U·(−∇h·R)>_h = <R_ij ∂j U_i>_h.

For periodic horizontal boundaries and a conservative/dealiased discretization,

    dK'/dt = −P_BT + <wᵀ C b>_h − Dν.

Dν is the positive kinetic dissipation computed using H, W, Vh and Vw. Mean advection, Coriolis acceleration and pressure contribute zero global work individually. The mean-strain work is exactly the opposite of P_BT. This identity remains valid with a prescribed mean; it is a prospective exchange, not an observed mean-flow energy increase. Verify both the stress-divergence and stress-strain evaluations independently.

The buoyancy-variance identity is

    d(1/2 <bᵀb>_h)/dt = −<bᵀ B w>_h − Dκ − <bᵀΓb>_h.

For uniform nonzero N², B=N² Cᵀ and K'+<bᵀb>/(2N²) is the inviscid, unforced perturbation energy in the absence of mean shear. For unstable N² that invariant is indefinite. For a sign-changing capped N² profile, do not divide by N² or call buoyancy variance positive available potential energy. Instead retain the two explicit work identities. All matrix quadrature must respect the same depth normalization.

## Reflection and matched controls

Use a horizontal reflection, for example (x,y,z)→(−x,y,z), with

    (U_x,U_y)→(−U_x,U_y),  (u_x,u_y,w,b)→(−u_x,u_y,w,b),
    f→−f,

and evaluate the right-hand fields at (−x,y,z). Vertical vorticity reverses sign, while cyclone/anticyclone identity relative to f is unchanged. K', P_BT and buoyancy work are invariant. This is an exact equation symmetry; simply changing f and negating ζ at identical coordinates is not the full symmetry transform for a generic disturbance.

At fixed f, equal-amplitude opposite-sign vortices are *not* symmetry-equivalent and are the physical comparison. Match the disturbance distribution and energy across the two vortices, using several predeclared seeds or orientation phases to avoid an accidental seed alignment. A compensated vortex satisfies periodic zero net circulation; its extended weak opposite-sign ring is part of the state, not another independent realization.

For an axisymmetric profile V(r), require the inertial discriminant (f+2V/r)(f+ζ)>0, including its regular r→0 limit, before labeling it centrifugally stable. This condition does not exclude all shear or three-dimensional instabilities. For a periodic approximation, quantify departure from axisymmetry and retain a margin from zero rather than interpreting rounding-level positivity as stability.

## Release checks and interpretation bounds

The following are proposed numerical gates, not universal physical thresholds:

- Matrix positive definiteness and quadrature refinement; all matrix entries converge relative to their relevant nonzero scale, with target 10⁻⁸.
- Projection: relative divergence and idempotence/energy-adjoint residuals below 10⁻¹¹ for Float64 synthetic states; scale the divergence by |k| and q times velocity.
- Instantaneous work: independently evaluated pressure/Coriolis/mean-advection zero work, equal/opposite strain exchange, and full kinetic/thermal identities below 10⁻¹⁰ relative to the sum of absolute contributing powers. Use an absolute floor for zero-power states.
- Uniform-background wave/convective-onset benchmark and temporal refinement. A mixed sine/cosine thermal truncation does not reproduce the single-sine analytic dispersion exactly at finite thermal rank; test convergence against a sufficiently resolved independent vertical reference rather than weakening a tolerance until it passes.
- Full sign-reflection tests should agree to roundoff before chaotic or nonlinear effects are introduced. For this linear finite-window problem, relative histories at 10⁻⁹ or better are a reasonable target.
- Short production-driver output and cost check before releasing the planned response cases.
- The one-/two-velocity-mode comparison is an *adequacy* test. Both must reproduce the independently resolved vertical instability that supplies the intended convection. If a full-height low-rank basis suppresses it, do not tune the thermal gradient to compensate.

Leading vertical eigenfunctions may be complex. Include independent real/imaginary components when needed, and declare the resulting count honestly. The horizontal basis augmentation above must retain the selected toroidal profiles as well as g derivatives. Test the full eigenvector residual and growth rate after projection, then test the basis across the declared horizontal disturbance band: optimizing at one horizontal wavenumber alone does not establish adequate response after a vortex scatters disturbances to other wavenumbers.

Report actual and fluctuation-energy-normalized P_BT, perturbation amplification, and a linear-validity measure such as perturbation speed relative to vortex speed. Positive growth alone is not positive vortex reinforcement. If both vortices lose energy, the model has not explained their maintenance on its own. A robust reinforcement difference while both remain centrifugally stable would establish a useful candidate feedback; it would still not demonstrate spontaneous LSV formation or saturated convection.

The capped reference originates in an anelastic model, while these equations use Boussinesq density weighting. That approximation must remain explicit. A later quantitative comparison can require a density-weighted projection; a successful Boussinesq test alone does not validate it for the full density contrast.
