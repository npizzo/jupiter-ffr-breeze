# Polar vortex and convection pilots — November 10, 2026

The first staged batches are submitted on Unity. Submission does not establish
physical validation or completed science. All computation, diagnostics and
rendering remain on Unity; existing campaigns and their sources are preserved.

| Campaign | First batch | Compute | Required before production |
|---|---|---|---|
| [Polar γ-plane reference](gamma2026_v2.md) | Three polar/flat seed pairs at 1024²; one matched pair at 2048²; 75 Earth days per case | One L40S per job | CPU checks and separate GPU agreement, stability and throughput checks at each grid |
| [Dry Breeze convection](convection2026.md) | Matched 0, 54, 162 W/m² cases at 512² × 96; ten Earth days per case, in two-day segments | One L40S per job | Initialization, restart, thermal, cap, tracer, resource and throughput checks |
| Breeze refinement pilot | One 162 W/m² case at 1024² × 96; six simulated hours | One H100 | Successful 512² pilot |

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
  and the spectral-transfer self-tests passed. GPU pilot **1394137** follows;
  fine pilot **1394138** and production jobs **1394139–1394153** depend on it.
- γ-plane replacement CPU validation is **1397848**. The 1024² GPU gate is
  **1397849**, followed by six cases **1397850–1397855**. The 2048² gate is
  **1397856**, followed by two cases **1397857–1397858**.
- The preliminary γ-plane chain stopped at an overly strict bitwise FFT equality
  check; its dependent jobs were automatically cancelled before any science run.
  The isolated replacement uses a roundoff tolerance and the recovered reference
  initialization. A subsequent submission referencing the aged-out completed
  environment job was rejected before creating any job; the replacement uses
  the verified installed environment directly.

The [machine-readable launch record](new-campaigns-launch.json) gives every
dependency and the latest recorded validation status. Output roots are
`data/gamma2026-v2a`, `data/convection2026-v1` and
`data/convection2026-preflight-v1` on Unity. Failed gates stop their dependent
chain; they do not alter older experiments.

## Decisions after this batch

Assess energy and enstrophy losses, segregation, seed sensitivity and resolution
before extending γ-plane integrations toward mature crystals. Assess thermal
balance, resolved convective scales, energy transfer and measured GPU cost before
extending the Breeze ensemble or adding matched mechanical forcing. Stronger
heating is a sensitivity experiment, not an inferred Jovian heat flux.

Ferrel-cell production remains a later stage: the new Breeze outputs add zonal
means and covariance transports, but a closed momentum budget, residual
circulation diagnostics, jet controls and boundary sensitivity are still required.
Multiyear γ-plane runs, 4096² refinement, moisture and long convection ensembles
have not been submitted in this first batch.

New Breeze results will include synchronized three-depth vorticity/tracer movies
on the [GitHub Pages gallery](https://npizzo.github.io/jupiter-ffr-breeze/).
γ-plane movies will be labeled two-dimensional. Morphology alone will not be
used as evidence of inverse energy transfer.
