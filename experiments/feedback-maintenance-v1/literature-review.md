# Signed vortices: an adversarial review and a bounded alternative

November 18, 2026. This review accompanies the [maintenance experiment](README.md). It proposes a future model and response test; it does not authorize or report their execution. Scientific computation remains on Unity. The f-plane sign-selection question stays separate from beta-plane Ferrel-cell work.

## What the polynomial establishes

The present polynomial is a useful phenomenological hypothesis, but its preferred sign is built into its source. For periodic, zero-mean vorticity and

\[
Q_2=bf(\zeta^2-\langle\zeta^2\rangle),
\]

our direct calculation gives

\[
\left.\frac{d\langle\zeta^3\rangle}{dt}\right|_{Q_2}
=3bf\left[\langle\zeta^4\rangle-\langle\zeta^2\rangle^2\right].
\]

The bracket is the variance of squared vorticity and is nonnegative. Thus the quadratic term prescribes the sign of its third-moment production. Whether that bias organizes a coherent large vortex remains an emergent question. Skewness and its reversal under a reflected change of f verify implementation, but cannot establish the source's convective origin.

Mean subtraction preserves circulation; it does not derive a momentum transport. A solenoidal velocity force producing Q involves an inverse Laplacian, so a local vorticity reaction corresponds to a spatially nonlocal momentum force. This is not automatically unphysical, but its energy exchange and scale dependence require explanation.

A Ginzburg–Landau coefficient describes a specified amplitude near a bifurcation. It cannot be transferred directly to the full turbulent vorticity field without identifying the amplitude, expansion parameter, eliminated dynamics, and conservation constraints.

## Why the stress experiment is not a decisive rejection

The previous replacement changed three things together: 512² became 256², band-limited instability forcing became a projected stress acting throughout its retained band, and local polynomial reactions became derivatives of a nonlinear stress. The scale support and saturation mechanism therefore changed alongside the proposed physical closure. Less convincing morphology cannot isolate which change mattered.

The stress implementation's conservation and energy-work checks remain useful. They do not validate its assumptions of strain alignment, locality, short memory, or a weak dependence on local vorticity. Nor does one seed establish robust sign selection. The current maintenance pair is narrower: it changes only the quadratic output into the declared lowest modes of an already developed state.

## What convection studies support

[Guervilly, Hughes and Jones (2014)](https://arxiv.org/abs/1403.7442) discusses plume asymmetry and clustering alongside possible anticyclonic instability and weakening of rotational constraint. Its proposed mechanism is not a settled exclusion of alternatives. Particularly relevant is its suggestion to examine a few convective structures interacting with imposed large vorticity patches.

[Rubio et al. (2014)](https://doi.org/10.1103/PhysRevLett.112.144501) demonstrates feedback between large barotropic flow and convective fluctuations. The unresolved object in a depth-averaged equation is a correlated stress, not merely the distribution of independent random kicks. A [public author manuscript](https://atoc.colorado.edu/~jweiss/website/publications/RubioEtAl2014.pdf) is available.

[Favier, Guervilly and Knobloch (2019)](https://arxiv.org/html/1901.04269v1) finds subcritical condensates: a sufficiently strong initial vortex can reorganize convection into a state that maintains it. The study identifies nonlocal transfer from small three-dimensional motions to large two-dimensional motions, involving correlations absent from prescribed forcing. Consequently, dependence on direct large-scale feedback is not itself evidence of an artifact. Maintenance and spontaneous formation are separate tests.

[Van Kan et al. (2025)](https://arxiv.org/html/2409.08536v3) finds that cyclone–anticyclone asymmetry disappears toward the asymptotic rapid-rotation limit along its studied parameter path. The reported transition is not a universal threshold for our configuration. A finite-Rossby-number explanation must retain the terms that break the leading-order symmetry.

Constant f disappears from the strictly two-dimensional incompressible vorticity equation. Adding a height-independent third velocity component does not supply the missing feedback: a useful extension requires nonzero vertical structure and its associated divergent horizontal motion.

## The strongest alternative: retain an evolving vertical mode

[Thomas and Vishnu (2022)](https://doi.org/10.1175/JPO-D-21-0270.1) projects hydrostatic Boussinesq dynamics onto a barotropic mode and one high baroclinic mode. At increasing Rossby number, a weak baroclinic component modifies coherent structures with preferential cyclonic survival. Importantly, it also enhances forward transfer and erodes vortices. This is a derived, inexpensive example of sign-sensitive feedback, not a recipe for convective LSV formation.

The publisher's indexed abstract was accessible during this review although direct full-text access returned an error. The [author's publication page](https://jimthomasedu.github.io/publications/) supplies an [author PDF](https://github.com/jimthomasedu/jimthomasedu.github.io/raw/master/files/TV22.pdf). The [2024 primary follow-up](https://journals.ametsoc.org/view/journals/phoc/54/9/JPO-D-23-0132.1.xml) presents the same model equations. Its imposed forcing and hyperdissipation are not adopted here.

[Benavides and Alexakis (2017)](https://arxiv.org/html/1701.05162v1) provides a complementary precedent: a horizontal flow coupled to one nonzero vertical mode captures changes in cascade direction. It is neither rotating nor convective, so it supports the projection architecture rather than the proposed physical mechanism.

**Our proposed extension** is a finite-f, nonhydrostatic Galerkin model retaining barotropic velocity, one or two vertical velocity/temperature modes, and an evolving mean thermal structure. It is not an established convection closure. For an illustrative stress-free basis only,

\[
\mathbf u_h=\mathbf U+\mathbf u_1\cos(qz),\qquad
w=w_1\sin(qz),\qquad
\nabla_h\cdot\mathbf u_1=-q w_1,
\]

the depth-averaged horizontal stress is

\[
\mathbf R=\tfrac12\mathbf u_1\mathbf u_1,
\qquad
\partial_t\mathbf U+\mathbf U\cdot\nabla_h\mathbf U
=-\nabla_hP-\nabla_h\cdot\mathbf R+\nu\Delta_h\mathbf U.
\]

The fluctuation equations retain Coriolis acceleration, pressure, buoyancy, and vertical acceleration. Their dynamics determine the stress and its sign instead of prescribing a signed local source. Constant-f barotropic Coriolis acceleration can be absorbed into pressure, while the divergent fluctuation component still responds to rotation.

This does not justify replacing Breeze's cap by a sinusoid. The eventual vertical basis and reference state must represent the existing stratification and boundaries. Buoyancy flux also changes the thermal profile: at minimum retain its mean/second-harmonic response so convection can saturate. A permanently fixed unstable gradient and a lone unsaturated linear mode would be inadequate. One-versus-two nonzero vertical modes is an adequacy check, not a claim of convergence.

## A small response experiment before a long campaign

1. **Verify the projected equations.** Check continuity, inviscid nonlinear exchange between barotropic and fluctuation energy, and heat transport. Retain ordinary viscosity and thermal diffusion; introduce neither beta nor linear drag. Do not launch until the model's assumptions and bounded cost are stated.
2. **Prescribe matched compensated vortices.** Compare equal-amplitude cyclonic and anticyclonic profiles with matched perturbations and thermal conditions. First choose both centrifugally stable. For an axisymmetric horizontal vortex check the epicyclic discriminant
   \[
   (f+2V/r)(f+\zeta)>0,
   \]
   rather than applying the planar-shear criterion alone. This does not exclude every possible instability.
3. **Measure growth and feedback separately.** Calculate disturbance growth and barotropic power
   \[
   P_{\rm BT}=\langle\mathbf U\cdot[-\nabla_h\cdot\mathbf R]\rangle.
   \]
   Compare actual power and power normalized by fluctuation energy. A rapidly growing disturbance need not maintain its parent vortex. Use finite-window response initially; a linearly unstable subsystem does not have stationary statistics merely because averaging is performed.
4. **Test the proposed cause.** Include the reflected reversal of f, an unheated counterpart, and the second vertical mode. Distinguish convective growth from generic wave–vortex response. Only then consider a freely evolving barotropic flow.

Stop before long integrations if energy exchange fails, fluctuation growth is unsaturated, the result depends qualitatively on the vertical truncation, or no reproducible sign difference appears. If both signs lose energy, this candidate does not yet explain preferential maintenance. A robust preference while both vortices remain centrifugally stable would challenge centrifugal instability as the complete explanation. None of these short tests alone proves spontaneous LSV formation.

## Other ideas and their limits

[Gowthaman and Thomas (2024), S-2DV](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004438), uses a quadratic modification inside a Jacobian, rather than a reaction, and reproduces cyclone-preferred morphology. The authors explicitly describe the construction as ad hoc. Our calculation shows that without its external field, forcing, or dissipation, smooth periodic solutions conserve every integral of a function of vorticity. It is an interesting transport control, not a derivation of convective skewness production.

[Parfenyev et al. (2021)](https://arxiv.org/html/2104.09914v1) derives different cyclone and anticyclone profiles through fluctuation interaction with an established rotating vortex. This makes wave-mediated maintenance a concrete alternative to selective centrifugal destruction, but its mechanical forcing is not convection.

[Marston, Chini and Tobias (2016)](https://doi.org/10.1103/PhysRevLett.116.214501) supplies a generalized quasilinear hierarchy for selected mean–fluctuation interactions. Its utility here would require testing the omitted interactions, not assuming them negligible. [Bouchet and Simonnet (2009)](https://doi.org/10.1103/PhysRevLett.102.094504) explains sensitivity of condensate topology to geometry and state; morphology alone therefore cannot diagnose a cascade or an f-dependent mechanism.

The recommended sequence is to complete the present maintenance intervention, then test an explicitly derived fluctuation response before investing in another long phenomenological parameter campaign.
