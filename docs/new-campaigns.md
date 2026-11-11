# Polar vortex and convection pilots — November 10, 2026

The first staged batches are submitted on Unity. Submission does not establish
physical validation or completed science. All computation, diagnostics and
rendering remain on Unity; existing campaigns and their sources are preserved.

| Campaign | First batch | Compute | Required before production |
|---|---|---|---|
| [Polar γ-plane reference](gamma2026_v2.md) | Three polar/flat seed pairs at 1024²; one matched pair at 2048²; 75 Earth days per case | One L40S per job | CPU checks and separate GPU agreement, stability and throughput checks at each grid |
| [Dry Breeze convection](convection2026.md) | Matched 0, 54, 162 W/m² cases at 512² × 96; ten Earth days per case, in two-day segments | One L40S per job | Initialization, restart, thermal, cap, tracer, resource and throughput checks |
| Breeze refinement pilot | One 162 W/m² case at 1024² × 96; six simulated hours | One H100 | Successful 512² pilot |
| [Ferrel mechanism screening](ferrel2026.md) | Matched 0/162 W/m² cases and a combined boundary sensitivity at 512² × 96; 48 hours each | L40S; H100 for an independent 1024² check | Own initialization, restart, native momentum/mass closure, axisymmetric controls, resource and throughput gates |

The γ-plane experiment is an independent, supplementary-information-informed
reproduction attempt using GeophysicalFlows SingleLayerQG. It includes the
published quartic taper and ±5% initial spectral band. It has no convection.
The 75-day stage tests early segregation and numerical consistency; mature
vortex crystals require much longer integrations. Current filtering, new random
seeds and adaptive timesteps distinguish this implementation from the historical
calculation.

The Breeze experiment uses a 6000 × 6000 × 120 km regional channel at 85°N,
Jupiter's rotation, a local f-plane, 1200-km initial eddy wavelength and Ro=0.1.
It tests dry bottom-heating response and practical resolution. Its cap relaxation,
sponge and global thermal feedback remain idealizations. It is not yet a
γ-plane convection experiment or a Ferrel-cell simulation.

## Submitted jobs and validation

- Breeze CPU initialization job **1394136** completed: 52 initialization checks
  and the spectral-transfer self-tests passed. GPU pilot **1394137** also passed
  all three heating cases, restart/thermal/tracer/cap checks and resource checks.
  Fine pilot **1394138** is queued; first production segment **1394139** started.
  The remaining production chain ends at **1394153**.
- γ-plane CPU validation **1397848** passed all 17 checks, but GPU gate
  **1397849** failed because CPU-side precompilation had not selected a CUDA
  runtime. Its science jobs were cancelled automatically before starting.
  The [isolated CUDA repair](gamma2026_gpu-runtime.md) passed on an L40S in
  **1401300**, with unchanged numerical sources and resolved package versions.
  Replacement CPU validation is **1401301**, GPU gates are **1401302/1401309**,
  and science cases are **1401303–1401308, 1401310–1401311**.
- Ferrel initialization **1401859** passed. Quick numerical closure/restart
  validation is **1401860**, followed by four full-grid six-hour screens in
  **1401861**. The fine check is **1401862**; the three 48-hour cases are
  **1401863–1401865**, conditional on the successful full-grid gate.
- The preliminary γ-plane chain stopped at an overly strict bitwise FFT equality
  check; its dependent jobs were automatically cancelled before any science run.
  The isolated replacement uses a roundoff tolerance and the recovered reference
  initialization. A subsequent submission referencing the aged-out completed
  environment job was rejected before creating any job; the replacement uses
  the verified installed environment directly.

The [machine-readable launch record](new-campaigns-launch.json) gives every
dependency and the latest recorded validation status. Current output roots are
`data/gamma2026-v2gpu1`, `data/convection2026-v1`,
`data/convection2026-preflight-v1`, `data/ferrel2026-v1` and
`data/ferrel2026-preflight-v1` on Unity. Failed gates stop their dependent
chain; they do not alter older experiments.

## Decisions after this batch

Assess energy and enstrophy losses, segregation, seed sensitivity and resolution
before extending γ-plane integrations toward mature crystals. Assess thermal
balance, resolved convective scales, energy transfer and measured GPU cost before
extending the Breeze ensemble or adding matched mechanical forcing. Stronger
heating is a sensitivity experiment, not an inferred Jovian heat flux.

The Ferrel pilot now implements sampled native momentum closure, a mass-consistent
streamfunction, weak balanced jets, axisymmetric controls and a combined boundary
sensitivity. Its 48-hour tests assess early numerical and mechanism credibility;
steady circulation, statistical significance and Jovian fidelity remain unproven.
Direct tracer fluxes are used in the neutral layer, where standard TEM formulas
are ill-conditioned. Extended experiments require satisfactory budget closure,
eddy-to-mean transfer, control comparisons and refinement.
Multiyear γ-plane runs, 4096² refinement, moisture and long convection ensembles
have not been submitted in this first batch.

New Breeze results will include synchronized three-depth vorticity/tracer movies
on the [GitHub Pages gallery](https://npizzo.github.io/jupiter-ffr-breeze/).
γ-plane movies will be labeled two-dimensional. Morphology alone will not be
used as evidence of inverse energy transfer.
