# 35°N continuation through 300 Earth days

Authorized November 11, 2026. Continue the heated 35°N beta-plane case from its full day-100 checkpoint through day 300. Twenty additional ten-day jobs depend on successful completion of the original day-100 job 1435756. They retain the same grid, initial realization's trajectory, latitude, beta, heating, cap, controller, diffusivities and saved output cadence. No initial state is regenerated as a replacement for checkpoint pickup.

The extension uses the unchanged, hash-checked `scripts/cascade35` model source. Its separate continuation driver verifies the original GPU gate, source/configuration/environment hashes, previous segment status and time, and full checkpoint before each segment. It saves into `data/cascade35-300-v1`; original day 0–100 output and checkpoints are preserved. Stable-cap, thermal-accounting, finite-state, tracer, zonal-output and transfer checks are retained. A failure blocks dependent segments.

Each segment requests one L40S GPU, four CPUs and 32 GB RAM, with a six-hour time limit. The existing day-100 comparison with 55°N remains queued independently. The continuation does not wait for that comparison movie to render.

After day 300, a dependent Unity CPU job joins the original and continued outputs into a read-only view, verifies all 2401 saved times, and renders the full 0–300-day movie with three depth columns, vorticity above passive tracer and fixed documented scales. It also produces zonal-wind time-versus-y diagrams, zonal kinetic-energy fractions and 50-day block profiles for assessing persistent jets and late-time drift. Full media decoding is required before publication to the [GitHub gallery](https://npizzo.github.io/jupiter-ffr-breeze/#beta35).

[Submission record](cascade35-300-launch.json). All computation, diagnostics, rendering and media validation remain on Unity. Submission is not a completed validation or a science result.

Day 500 is a conditional follow-up, not part of this submitted chain. Review day-300 numerical/physical checks and jet evolution before extending. Block statistics describe drift; they do not provide independent-sample confidence intervals or prove equilibrium. Longer duration does not address inadequate domain size, resolution, or deep spherical circulation absent from this 120-km channel.
