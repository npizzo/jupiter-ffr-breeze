# Centrifugal-instability audit of the existing mature f-plane LSV

This is postprocessing only. No simulation is launched, restarted or modified.
All calculations and figure rendering run on Unity. The source is
`data/fplane-point4-v2/late_heated`, the completed days-90–100 heated continuation
at 256 × 256 × 64, with the original capped dry anelastic physics and constant
positive f. Eleven daily full 3D exports and 81 three-hourly movie states exist.

The analysis follows the questions raised by Favier, Silvers & Proctor (2014),
section IV.E and figure 11: https://arxiv.org/pdf/1408.6483. It tests whether
their proposed centrifugal-disruption explanation is supported by this run;
it does not assume that the model configurations are equivalent.

## Vorticity distributions

Native staggered velocities are interpolated to cell centers and differentiated
with centered differences, periodic in x and free-slip parity in y. This matches
the saved movie diagnostic; all eleven coincident fields are checked against
the separately stored movie slices. Native face-corner curl is also computed,
and threshold fractions are reported separately to expose averaging sensitivity.

Statistics use all x and the central 80% of the channel width. The outer wall
strips are stored separately. Each snapshot and horizontal cell receives equal
weight at a given height; the vertical vorticity PDFs are not density weighted.
Raw and central skewness are both reported. Fixed vertical bands are 0–10,
10–60, 60–80 and 80–120 km, with the cap beginning at 80 km. These geometric
bands are not fitted thermal-boundary-layer thicknesses. The lowest saved cell
is centered at 0.9375 km, with vertical spacing 1.875 km.

Histograms span −4 to 4 in zeta/f, with explicit outside-range fractions.
Plots zoom to −2 to 2. Exceedances of both −1 and +1, actual extrema, and
mean-profile N² are measured without smoothing. A PDF bend is visual evidence,
not an automatically established instability threshold or fitted breakpoint.

## Vortex profiles and discriminant

At each height and daily state, select up to three strongest local anticyclonic
extrema with zeta/f ≤ −0.25, separated by the nine-cell local-maximum stencil,
whose horizontal velocity gradient has complex eigenvalues. The latter is a
rotation-dominance screen, not an instability criterion. The search is confined
to the same central region so sampled circles do not cross sidewalls.

Sample 96 azimuths on each of twelve radial rings (1–12 cells) around each
center. Subtract the center velocity before assessing axisymmetry; full-circle
mean azimuthal velocity is unchanged by a uniform translation. Compute

    Phi = (f + 2V/r)(f + zeta_ring),

and independently replace zeta_ring by (1/r) d(rV)/dr. Ring means and linear
interpolation are diagnostics of the existing fields, not a new filtered model.

The declared screen requires both discriminants below −0.01 f² on at least two
adjacent rings between 4 and 10 cells, plus an axisymmetric azimuthal-energy
fraction at least 0.5 over those rings. Results below four cells are excluded;
no resolved-vortex claim is made there. This is a conditional inviscid,
axisymmetric, barotropic screen. Failure does not exclude smaller, non-axisymmetric,
stratification-dependent or otherwise untested instabilities. Passing does not
establish a growing eigenmode in the anelastic convective simulation.

Analytic checks include uniform flow, periodic differentiation, three solid-body
vortices (including a strong anticyclone), a localized Gaussian vortex with a
negative-discriminant annulus, sign/rotation covariance and independent native
versus centered curl equivalence. The Gaussian error tolerance applies at the
same 4–10-cell radii used in the screen. Smaller-radius interpolation errors are
recorded separately rather than treated as resolved physics.

## Temporal and energy evidence

Three-hourly horizontal velocities exist at the lowest movie level only
(19.6875 km). Up to eight rotation-dominated extrema per sign/frame are matched
one-to-one within four cells of a center advected by the local depth-mean
velocity, with a permitted strength ratio of 0.25–4. These are tentative pattern
tracks. Endpoints can reflect merging, weakening, rank truncation or a failed
match; they are not measured instability-induced deaths or unbiased lifetimes.

Daily full fields give a reference-density-weighted depth mean U and deviations
u'. The centered diagnostic

    P = -u'_i u'_j partial_j U_i,  i,j in {x,y},

is positive for transfer from horizontal depth-mean flow to deviations. The
companion thermal-production proxy is g w(theta−theta_bar)/theta_bar. These
are computed at cell centers, with mass-weighted regional averages and
vorticity-conditioned statistics. They are not a closed native perturbation-
energy budget; the central-region boundaries permit transport, and pressure,
diffusion, numerical and forcing terms are not fully reconstructed here.

The inverse-f time is about 0.96 hours and the inertial period about 6.06 hours.
Daily 3D sampling cannot resolve centrifugal growth on the former scale.
Three-hourly single-level velocities do not recover missing 3D motions. This
audit can identify susceptibility and correlations, but cannot turn those into
a causal growth/saturation demonstration without additional resolved evidence.

## Publication

The gallery uses the already published three-depth vorticity/tracer movie of
this same days-90–100 run. Raw 3D exports and local map arrays stay on Unity.
Figures, measured diagnostics, method, source hashes and input-file checksums
are published. No other campaign is changed.
