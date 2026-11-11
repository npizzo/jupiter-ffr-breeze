# Completed Ferrel publication observer — layout revision 2

Run only in a Unity allocation. This tool reads completed output without modifying
scientific sources, fields or checkpoints. No data need to leave Unity except the
finished publication assets. No script here submits or publishes by itself.

The default input roots are `data/ferrel2026-gpu-v1` and
`data/ferrel2026-preflight-gpu-v1`. `--stage pilot` reads the four 512² six-hour
pilot cases. `--stage production` reads the three 48-hour cases and the two
six-hour axisymmetric controls. `--include-fine` additionally requires the
successful `PASS1024.toml` gate and the independent six-hour refinement.

Example full production submission, after the existing simulation chain:

```bash
sbatch --kill-on-invalid-dep=yes --dependency=afterok:1409345:1409342 \
  slurm/ferrel2026-publication-v2.sbatch \
  --stage production --include-fine --output rendered/ferrel2026-production-v2 \
  --completed-job 1409341 --completed-job 1409342 \
  --completed-job 1409343 --completed-job 1409344 --completed-job 1409345
```

For a separate, clearly labeled pilot report, use the pilot's successful job:

```bash
sbatch --kill-on-invalid-dep=yes --dependency=afterok:1409341 \
  slurm/ferrel2026-publication-v2.sbatch \
  --stage pilot --output rendered/ferrel2026-pilot-v2 --completed-job 1409341
```

Use a fresh output directory for every rendering attempt. Successful prerequisite
Slurm exit states, complete metadata, complete chronological frame coverage,
source/configuration fingerprints, native mass reconstruction, finite fields,
mean/eddy tracer decomposition and symmetric-control properties are required.
The caller must list all relevant successful simulation jobs; metadata COMPLETE
alone is not proof that post-run validation passed.

Each case gets an actual saved-plane movie with three columns and two rows:
vorticity above passive tracer. Heating/control, boundary-sensitivity and
axisymmetric pairs additionally receive synchronized four-row comparisons.
All cases share fixed-in-time scales: symmetric asinh vorticity in ζ/f0 with
width 0.1 and full per-depth extrema; asinh tracer with width 0.01 and physical
concentration ticks. The tracer range covers [0,1] and any observed overshoots.
Every saved frame appears once. The default six frames/s gives 8.2-second movies
for 49 hourly 48-hour frames. Seven-frame pilot/control movies use one frame/s
so their early evolution remains watchable. H.264 yuv420p, faststart, decode/frame-count/duration
validation, and SHA-256 manifests are mandatory. Assets exceeding 24 MB fail
publication preparation rather than silently dropping data. Preview the encoded
movies before publishing; numerical/decode checks do not establish legibility.

Diagnostics include native mass streamfunction (kg/(m s), per unit zonal length),
sampled exact SSPRK3 momentum closure, approximate centered eddy→mean energy
transfer with its native-advection remainder, and direct mean/eddy tracer fluxes.
Time means use trapezoidal integration of saved snapshots. The first 12 hours are
excluded from 48-hour summaries as an explicit startup convention; equilibration
is not demonstrated. Six-hour pilots/controls/refinement are reported separately
within this early window, with coarse/fine and axisymmetric comparisons restricted
to coincident times. Variability/coherence metrics are descriptive, not confidence
intervals. No inverse-cascade or steady-Jovian-cell claim is generated.

Output includes movies, posters, diagnostic figures, compact CSV tables,
`report.json`, a gallery `README.md`, and `publication-manifest.json`. Only a
READY_FOR_REVIEW manifest means rendering finished; the parent task must still
inspect visuals, write a scientifically qualified report, publish via GitHub,
update the existing Pages landing page without replacing other content, and
verify public playback.

## Layout-only change from revision 1

Revision 1 passed allocated reader/encoding QA (job 1411081), but its automatic
movie x-axis labels overlapped in the three-column frame. Revision 2 changes only
movie tick placement: four evenly spaced x ticks (0, 4000, 8000, 12000 km for this
domain), and five evenly spaced y ticks. Fields, normalization, frames, movie size,
encoding, diagnostics, validation and scientific sources are unchanged.

Run `slurm/ferrel2026-publication-probe-v2.sbatch` to check the same completed
unheated six-hour case. It guards this revision's renderer hash and writes only
to `rendered/ferrel2026-schema-probe-v2`. The parent must inspect the encoded
movie/poster; only after that visual check should the published revision 1 layout
be marked superseded. Existing revision 1 source, observer jobs and outputs remain
untouched. Replacement observer jobs use fresh revision 2 output directories.
