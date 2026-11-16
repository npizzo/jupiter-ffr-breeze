# Random stirring: random core sign versus core sign fixed to f

Two additional 2D Oceananigans controls, requested by the user, complement the
existing instability-forced and unforced long runs. Both start from the same
completed pilot checkpoint at sigma*t=120 and integrate to 1200 on Unity.
The existing two runs and their source files are unchanged.

## Meaning of the sign control

A periodic velocity field has zero domain-integrated relative vorticity. We
therefore cannot inject everywhere-positive vorticity while retaining this
periodic model. Each event injects a localized positive-core, negative-surround
vorticity template, or its negative. The template is the inverse transform of
unit coefficients on 6<=|k|<=8, excluding kx=0 as in the instability forcing.
It is normalized to unit RMS velocity using the compatible native curl/inverse.
Its zero Fourier mode vanishes; its core has positive vorticity and its third
moment is positive. The outer oscillatory lobes compensate its circulation.
The exclusion of kx=0 means the template is not perfectly isotropic.

At every event, translate the template to an independent uniformly random grid
location. The cases differ ONLY in the sign multiplying this template:

- `random_sign`: independent equiprobable positive/negative cores.
- `cyclonic_sign`: core sign is `sign(f)`, positive here because f=+1.

Both cases draw identical locations using seed 23017 and identical sign coins
using seed 23018, with the latter ignored by the fixed-sign case. Thus the
event locations stay paired. Multiplying a zero-mean Gaussian field by a fixed
global sign would not change its distribution; non-Gaussian pulses are required
for this polarity comparison. Random-sign pulses have a sign-symmetric law;
the f-aligned pulses have a prescribed odd-moment asymmetry but zero spatial
mean. This is an imposed-asymmetry response experiment, not a derivation of a
rotation-sensitive closure or an explanation of spontaneous Breeze selection.
In particular, for zero mean vorticity and an independent uniform translation,
the expected domain third-moment increment is zero for the random-sign law and
`sign(f) A^3 mean(template_q^3)` per event for the fixed-core law: the cross
terms average away. This prescribed third-moment input is recorded explicitly.
Whether its sign reaches a large-scale coherent vortex is a separate outcome.

## Energy and time control

Every 0.2 dimensionless time units (40 steps at dt=0.005), add a divergence-free
velocity pulse A*g. The template g has kinetic energy 1/2 and

    A = sqrt(2 epsilon * 0.2).

Epsilon is the instantaneous work of the original instability forcing on the
common branch checkpoint, measured on Unity with the compatible native operator.
It is fixed for both stochastic runs. The actual jump is

    Delta E = A <u dot g> + A²/2.

Uniform independent translations give zero EXPECTED cross work, even for fixed
core sign, because the template's spatial mean is zero. Both controls therefore
have identical expected injection, forcing spectrum, event cadence and amplitude.
Actual realized input fluctuates and is measured separately; it is not silently
rescaled to make each kick inject positive or equal energy. This also avoids
introducing state-dependent amplitude feedback. The reference rate matches the
instability run at branching, not necessarily its later mean work.

Physical parameters remain those of the long pilot: 128², periodic 2pi square,
constant f=1, no beta, viscosity 0.001, cubic drag coefficient 1, mean drag
removed and native centered advection. The pulse runs use dt=0.005, half the
original long-pair step, because a startup impulse exceeded CFL=0.5 at dt=0.01.
The forcing amplitude, event interval and physical damping are unchanged. Sigma=0.6 is the common plotting
clock, not a growth rate of stochastic forcing and not a physical day mapping.
The event is an exactly applied instantaneous jump between completed RK steps;
this defines a finite-cadence pulse process, not Gaussian white noise.

## Diagnostics and validation

Native RK budgets are closed before each pulse. Its exact vorticity, energy and
third-moment increments are then added as a separate forcing contribution.
Oceananigans caches next-step tendencies at the previous step's endpoint;
after a pulse the native tendencies and the initial-stage audit are refreshed
from the changed state. A startup test caught this ordering requirement before
production. The original budget tolerances were retained.
Independent total storage checks include both RK evolution and pulses. Each
event logs its time, location, sign, quadratic energy, flow/pulse cross work,
total work and third-moment increment. Exact per-event shell energy increments
are in `kick_spectra.csv`; their low-k values must vanish. Instantaneous
`spectra.csv` describes inter-pulse RK tendencies, so its injection column is
zero; it must not be mistaken for zero total forcing. Postprocessing combines
event-integrated injection with the inter-pulse advection/viscosity/drag data.

Checks cover zero circulation, discrete divergence, native curl, spectral
support, positive-core skewness, matched expected injection, paired random
locations, both sign laws, velocity-reflection/f-reversal covariance, exact
event and RK accounting, two 20-time-unit probes from the real checkpoint and
a timestep-halving comparison with unchanged physical pulse times.
No timestep or damping parameter is changed automatically. The CFL<0.5 guard
is checked after each pulse as well as at each RK stage.

As before, negative native advective energy flux below k=6 is evidence of
upscale transfer; a growing low-k fraction alone can reflect selective decay.
Spectra, absolute low-k energy, integral wavelength, skewness and signed budgets
are retained. With this short range below forcing, no broad inertial-range
claim follows. Finite-sample sign imbalance is recorded, not relabeled as a
preferred ensemble sign. One pair does not establish ensemble statistics.

## Outputs and bounded execution

Two single-CPU jobs with four-hour ceilings and internal 3.5-hour stops. Every
1.2 reference growth times, save a vorticity frame; 901 frames per completed
case. Checkpoints include both random-number-generator states and event counts.
The four-case renderer depends on BOTH old and new arrays completing and verifies
the common initial field, event schedule and color scales. It creates a 2x2
vorticity movie, spectra and flux plots on Unity. There are no invented depth
panels or tracers in this horizontal 2D model. Finished assets require review
before publication. No new 3D simulation, beta effect or anticyclone sink is added.

For context on why forcing statistics and nonlinear damping matter, see van
Kan et al. (2024): https://arxiv.org/html/2308.08789v2. The present localized
pulse prescription is a new control, not their exact stochastic forcing.
