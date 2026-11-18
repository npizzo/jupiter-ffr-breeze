# Nonlinear Galerkin convection: equations and interpretation

This experiment asks whether a small, fully coupled projection of rotating convection can develop a preference for one sign of depth-averaged vorticity. Its coefficients come from the fluid equations. There is no prescribed polynomial vorticity feedback or sign-dependent forcing.

This is a **fixed-temperature, free-slip Boussinesq reference problem**, not a quantitative reduction of the earlier capped Breeze configuration. The initial purpose is to establish a conservative, symmetry-correct nonlinear model before attempting that adaptation. The mode-selection reference is [Rademacher and Welter, *Energetic consistency and heat transport in Fourier-Galerkin truncations of free slip 3D rotating convection*](https://arxiv.org/abs/2505.07678), particularly Criteria 3.1 and 3.3 and equation 56. We use both horizontal phases and an anisotropic mode box; we do not impose their optional phase locks.

## 1. Physical problem and units

Take a right-handed coordinate system with upward vertical coordinate \(z\), and horizontal periodicity on \([0,L]^2\), \(0\leq z\leq H\). For positive \(f\), positive vertical vorticity is cyclonic. Write total buoyancy as

\[
B(x,y,z,t)=B_{\rm ref}(z)+\theta(x,y,z,t),\qquad
B_{\rm ref}(z)=-Gz,\quad G>0.
\]

The irrelevant additive constant in \(B_{\rm ref}\) may instead put zero buoyancy at the upper wall. The background buoyancy force is hydrostatic and is absorbed in the pressure. The equations are

\[
\partial_t\boldsymbol u+\boldsymbol u\cdot\nabla\boldsymbol u
+f\boldsymbol e_z\times\boldsymbol u
=-\nabla p+\theta\boldsymbol e_z+\nu\nabla^2\boldsymbol u,
\qquad \nabla\cdot\boldsymbol u=0,
\]

\[
\partial_t\theta+\boldsymbol u\cdot\nabla\theta
=Gw+\kappa\nabla^2\theta.
\]

Both walls satisfy

\[
w=0,\qquad \partial_z u=\partial_z v=0,\qquad \theta=0.
\]

Thus total buoyancy is fixed at both walls. The imposed unstable temperature difference supplies convective energy through boundary heat transport. It is **not** a prescribed uniform bottom heat flux, and there is no stable cap in this reference.

Lengths are scaled by \(H\), times by \(|f_0|^{-1}\), velocities by \(|f_0|H\), and buoyancy by \(|f_0|^2H\). Hence \(H=1\) and \(f=\pm1\). The dimensionless parameters are

\[
E=\frac{\nu_{\rm dimensional}}{|f_0|H^2}=\nu,
\quad Pr=\frac{\nu}{\kappa},
\quad Ra=\frac{G}{\nu\kappa},
\quad Ta=\frac{f^2}{\nu^2}.
\]

These times are inverse-Coriolis units, not Earth days or rotation periods. The first comparison uses \(Pr=1\); the executed configuration records the actual viscosity, supercriticality and box width. There is no beta term, momentum drag, hyperviscosity, random stirring, moisture or passive tracer.

## 2. Galerkin subspaces and nonlinear coefficients

Let \(k_0=2\pi/L\), \(\boldsymbol k_h=k_0(m_x,m_y)\), and use

\[
\begin{aligned}
(u,v)&=\sum_{\boldsymbol m,n}(u_{\boldsymbol m n},v_{\boldsymbol m n})
e^{i\boldsymbol k_h\cdot\boldsymbol x_h}\cos(n\pi z),\\
(w,\theta)&=\sum_{\boldsymbol m,n}(w_{\boldsymbol m n},\theta_{\boldsymbol m n})
e^{i\boldsymbol k_h\cdot\boldsymbol x_h}\sin(n\pi z).
\end{aligned}
\]

Reality requires conjugate horizontal coefficients. Incompressibility requires \(i\boldsymbol k_h\cdot\boldsymbol u_h+n\pi w=0\). There is no sine mode at \(n=0\). The horizontal-mean vertical velocity vanishes identically. The spatially constant horizontal velocity is initialized to zero and remains zero.

The retained subspaces are:

| Field | Nonzero horizontal wavevector | Horizontal mean |
|---|---|---|
| Horizontal velocity | \(|m_x|,|m_y|\leq K_h,\;0\leq n\leq N_v\) | Both components, \(1\leq n\leq2N_v\) |
| Vertical velocity | Same horizontal box, \(1\leq n\leq N_v\) | Zero by incompressibility |
| Buoyancy perturbation | Same horizontal box, \(1\leq n\leq N_v\) | \(1\leq n\leq2N_v\) |

The first model uses \((K_h,N_v)=(4,2)\), and the nested extension uses \((6,3)\), with the same physical box and coefficients. The mean thermal extension is essential: products of retained convective modes generate sum/difference vertical harmonics through \(2N_v\). Mean horizontal shear harmonics are also retained to preserve the integrated vorticity balance. Retaining both mean velocity components through \(2N_v\) slightly exceeds the sufficient upper index \(2N_v-1\) in the published full hierarchy. These are modest finite-dimensional systems, not a claim that a particular 13-mode truncation can describe an LSV.

Denote orthogonal projection onto the divergence-free velocity space by \(P_u\), and scalar projection by \(P_\theta\). The actual ODE is

\[
\dot{\boldsymbol u}=P_u[-\boldsymbol u\cdot\nabla\boldsymbol u
-f\boldsymbol e_z\times\boldsymbol u+\theta\boldsymbol e_z]
+\nu\nabla^2\boldsymbol u,
\]

\[
\dot\theta=P_\theta[-\boldsymbol u\cdot\nabla\theta+Gw]
+\kappa\nabla^2\theta.
\]

For an orthonormal basis \(\boldsymbol\phi_a,\chi_r\), the quadratic coefficients are integrals such as

\[
N_{abc}=-\langle\boldsymbol\phi_a\cdot
 (\boldsymbol\phi_b\cdot\nabla)\boldsymbol\phi_c\rangle,
\qquad
T_{rbs}=-\langle\chi_r(\boldsymbol\phi_b\cdot\nabla)\chi_s\rangle.
\]

They can be evaluated by exact Fourier convolution or sufficiently padded Fourier quadrature on a doubled vertical period of length 2. This doubled domain is an implementation device: horizontal velocity is even and \(w,\theta\) are odd in \(z\). The physical domain remains \(0\leq z\leq1\). Velocity, convective fluctuations and the horizontally averaged thermal profile all evolve. In particular, \(\overline\theta(z,t)\) is not reset to zero or artificially held fixed.

## 3. Exact instantaneous balances

Here \(\langle\cdot\rangle\) is a normalized average on the **physical** domain. Define

\[
K=\tfrac12\langle|\boldsymbol u|^2\rangle,\qquad
Q=\tfrac12\langle\theta^2\rangle,\qquad
C=\langle w\theta\rangle,
\]

\[
D_\nu=\nu\langle|\nabla\boldsymbol u|^2\rangle,\qquad
D_\kappa=\kappa\langle|\nabla\theta|^2\rangle.
\]

Orthogonality, incompressibility and the walls give

\[
\dot K=C-D_\nu,\qquad
\dot Q=GC-D_\kappa.
\]

Advection conserves kinetic energy and scalar variance separately; Coriolis does no work. Thus \(K-Q/G\) is an inviscid invariant when \(G>0\), but it is indefinite and must not be called a positive energy norm.

The additional potential-energy consistency check uses the linear moment

\[
M=\langle(1-z)\theta\rangle.
\]

Horizontal-mean \(w=0\) removes the \(Gw\) term from this moment. Integration by parts gives

\[
\dot M=-C+\kappa\langle(1-z)\nabla^2\theta\rangle,
\]

so

\[
\frac{d}{dt}(K+M)=-D_\nu+
\kappa\langle(1-z)\nabla^2\theta\rangle.
\]

The thermal diffusion term includes the wall heat exchange and must not be discarded. At finite truncation this identity needs the stated mean thermal modes; kinetic energy and scalar variance conservation alone do not establish it. For a mean sine expansion \(\overline\theta=\sum_n a_n\sin(n\pi z)\),

\[
M=\sum_n\frac{a_n}{n\pi},\qquad
\int_0^1(1-z)\sin(n\pi z)\,dz=\frac1{n\pi}.
\]

Using an unmodified linear weight over the entire doubled vertical domain gives the wrong physical moment. The implementation must use the analytic coefficient formula or independent accurate quadrature on the physical half-domain.

A useful boundedness check can also be derived. Put \(\ell=P_\theta(1-z)\), with the same retained mean thermal modes. Up to a constant, \(K+Q/G+2M\) equals the positive shifted quadratic expression

\[
\mathcal E=K+\frac{1}{2G}\langle(\theta+2G\ell)^2\rangle,
\]

whose derivative is

\[
\dot{\mathcal E}=-\nu\langle|\nabla\boldsymbol u|^2\rangle
-\frac\kappa G\langle|\nabla(\theta+G\ell)|^2\rangle
+\kappa G\langle|\nabla\ell|^2\rangle.
\]

This is a mathematical bound for the finite model, not a new physical energy reservoir. The finite-mode heating and diffusion remain able to change its value.

For completeness, the integrated vector-vorticity check is

\[
\langle\partial_t\boldsymbol\omega\rangle
=\langle(\boldsymbol\omega\cdot\nabla)\boldsymbol u
+f\partial_z\boldsymbol u+\nu\nabla^2\boldsymbol\omega\rangle.
\]

Its horizontal components can be nonzero because mean horizontal shear changes between the walls. The domain-average vertical vorticity is zero by horizontal periodicity. Mean shear modes in the selection above prevent the projection from dropping required boundary contributions.

## 4. Depth-averaged vorticity and reciprocal energy transfer

Write

\[
\boldsymbol U(x,y,t)=\int_0^1\boldsymbol u_h\,dz,\qquad
\boldsymbol u'_h=\boldsymbol u_h-\boldsymbol U,\qquad
R_{ij}=\int_0^1u'_iu'_j\,dz.
\]

The barotropic velocity is horizontally divergence-free. Since it has zero spatial mean, constant \(f\) is absorbed into its horizontal pressure after projection. Its equation is

\[
\partial_t\boldsymbol U+P_h(\boldsymbol U\cdot\nabla_h\boldsymbol U)
=-P_h\nabla_h\cdot R+\nu\nabla_h^2\boldsymbol U.
\]

Let \(\mathcal P_h\) denote retained horizontal Fourier-mode projection, \(Z=\partial_xU_y-\partial_yU_x\), and

\[
F_Z=-\mathcal P_h[\nabla_h\times(\nabla_h\cdot R)]_z.
\]

Then the finite model obeys

\[
\dot Z=-\mathcal P_h(\boldsymbol U\cdot\nabla_h Z)+F_Z+\nu\nabla_h^2Z.
\]

This is the proposed starting point for a minimal closure: diagnose \(F_Z\) from evolving convective modes rather than prescribe its sign dependence. Define

\[
E_{BT}=\tfrac12\langle|\boldsymbol U|^2\rangle_h,\qquad
P_{BT}=\langle R_{ij}\partial_jU_i\rangle_h.
\]

Then

\[
\dot E_{BT}=P_{BT}-\nu\langle|\nabla_h\boldsymbol U|^2\rangle_h,
\]

\[
\frac{d}{dt}(K-E_{BT})=C-P_{BT}
-\nu\langle|\nabla\boldsymbol u'|^2\rangle.
\]

Here \(\boldsymbol u'=(\boldsymbol u'_h,w)\). The same transfer must enter with opposite signs. A large energy fraction in \(E_{BT}\) is insufficient evidence of an inverse transfer: report its absolute energy and the signed power \(P_{BT}\).

## 5. Vorticity sign and a truncation caveat

A purely horizontal, incompressible constant-\(f\) model contains no dynamically active \(f\) in its vorticity equation. Any response to the sign of \(f\) here must therefore enter through the three-dimensional convective degrees of freedom and their stress \(R\).

The full equations have the exact reflection

\[
x\mapsto-x,\quad (u,v,w,\theta)\mapsto(-u,v,w,\theta),
\quad f\mapsto-f,
\]

with all fields evaluated at the reflected point. This maps \(Z\mapsto-Z\), preserves energies and reverses vorticity skewness. The negative-\(f\) control must use this reflected initial state, rather than merely changing \(f\) with arbitrary unmatched phases. Passing this test demonstrates the correct symmetry; it does not by itself establish a robust preferred sign.

Use \(S_Z=\langle Z^3\rangle_h/\langle Z^2\rangle_h^{3/2}\) only when barotropic enstrophy is appreciably above numerical noise. A single barotropic sine wave cannot have nonzero spatial skewness. Several noncollinear modes and zero-sum triads are needed; retaining conjugates and both horizontal phases avoids imposing a false sign preference by construction.

Crucially, finite Fourier-Galerkin advection generally does **not** conserve the cubic vorticity moment. Its exact budget is

\[
\frac{d}{dt}\langle Z^3\rangle_h
=-3\langle Z^2\mathcal P_h(\boldsymbol U\cdot\nabla_hZ)\rangle_h
+3\langle Z^2F_Z\rangle_h
+3\nu\langle Z^2\nabla_h^2Z\rangle_h.
\]

The first term cancels for the untruncated continuous equation but need not cancel after projection. It must be measured separately: attributing all third-moment production to convection would repeat a budget error. Quadratic energy and barotropic enstrophy conservation do not remove this issue. A candidate convective sign preference is credible only if it survives the nested truncation and cannot be explained solely by this projected-advection contribution.

## 6. Linear onset and bounded experiment selection

Before nonlinear integration, validate the retained linear modes independently. With \(q^2=k_h^2+(n\pi)^2\), the characteristic equation is

\[
(\sigma+\kappa q^2)\{q^2(\sigma+\nu q^2)^2+f^2(n\pi)^2\}
-Gk_h^2(\sigma+\nu q^2)=0.
\]

At \(\nu=\kappa\), its eigenvalues are

\[
\sigma_0=-\nu q^2,\qquad
\sigma_\pm=-\nu q^2\pm
\sqrt{\frac{Gk_h^2-f^2(n\pi)^2}{q^2}},
\]

and stationary onset for a specified mode is

\[
G_c(k_h,n)=\frac{\nu^2q^6+f^2(n\pi)^2}{k_h^2}.
\]

Select \(G\) from the retained-mode threshold and choose the horizontal fundamental so growing convection has larger wavenumber than the box mode. Record the eigenvalues for both truncations, not just a continuum critical estimate. The same physical parameters must be used in the extension; retuning its heating to force agreement would defeat the test. These low cutoffs provide only a few horizontal scales, so they cannot establish a broad inverse-cascade range.

The finite campaign comprises two positive-\(f\) seeds, the reflected negative-\(f\) counterpart of one seed, a matched zero-imposed-gradient control, and one nested mode extension of the first seed. The chosen parameters are \(k_0=20\), \(L/H=2\pi/20\), \(\nu=\kappa=10^{-5}\), and \(G=2\min G_c\) over the base model's retained nonzero horizontal convective modes and vertical indices. The minimum is evaluated on Unity, recorded, and held fixed in the extension. Integrations end at \(|f_0|t=200\), using \(\Delta t=0.02\), subject to the timestep and cost checks. The initial velocity and buoyancy RMS are \(10^{-4}\) and \(10^{-5}\), respectively, with no initially seeded barotropic or horizontal-mean fields. The two seed identifiers are 20261118 and 20261119. The extension embeds the first base initial state exactly, adding zero amplitudes for newly retained modes.

The \(G=0\) control retains the same initial buoyancy perturbation and its buoyancy force. It therefore has a finite initial thermal reservoir; it is not pure unforced Navier–Stokes decay. It tests the effect of maintaining an unstable background gradient.

Numerical validation and measured Unity cost precede release, with a projected total budget of 1800 CPU seconds for the nonlinear integrations. The first outcomes sought are finite nonlinear saturation, closed transfers, reproducible sign statistics and stability to added modes. Failure to generate a sign preference, failure of saturation, or loss of the preference under extension are useful negative outcomes. None authorizes an automatic longer or larger campaign. The small box accommodates only a few critical convective wavelengths and is expressly a finite-mode mechanism test, not a resolved Jovian layer or large-scale-vortex DNS.

## 7. What adapting to the capped Breeze setup would require

For a nonuniform reference profile, the thermal equation instead contains

\[
\partial_t\theta+\boldsymbol u\cdot\nabla\theta
=-w B'_{\rm ref}(z)+\nabla\cdot(\kappa\nabla\theta)
+\nabla\cdot(\kappa\nabla B_{\rm ref})+S,
\]

with the actual heating, cooling and cap terms in \(S\). Any part already balanced by the reference-state source must be accounted for once. A stable cap changes the coupling matrix through \(B'_{\rm ref}(z)\); it does not justify freezing the evolving mean buoyancy perturbation. Fixed flux walls need a compatible scalar basis and explicit boundary terms, rather than the homogeneous Dirichlet sine basis used here. Variable density additionally changes continuity, the pressure projection and the inner product used for conservation.

Those changes need fresh coefficient integrals and balance checks. The published free-slip, fixed-temperature consistency argument cannot simply be asserted for a capped or anelastic setup. Agreement in sign here would establish a mechanism worth testing in the original setup; disagreement would not by itself refute that setup's observed vortex.
