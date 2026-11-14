# 35°N horizontal-resolution audit, days 300–301

[Watch the synchronized comparison on the gallery](https://npizzo.github.io/jupiter-ffr-breeze/#ferrel35-resolution301).

Both Unity simulations completed, as did the independent review and full movie decode. Three synchronized depth columns show coarse vorticity/tracer above fine vorticity/tracer. Actual heights are 19.69, 59.06 and 100.31 km; all nine three-hourly states are shown at one frame per second.

The 256²×64 and 512²×64 cases retain the 6000×6000×120 km domain, original forcing and physical coefficients, with the same 15-second timestep ceiling. The original coarse 30-second diagnostic is a timestep control. Fine momentum is a conservative divergence-preserving refinement of the actual day-300 checkpoint; its initial kinetic energy is 1.7075% lower.

On common coarse meridional cells, the full-day reconstruction/eddy RMS ratio decreases from 0.452227 to 0.204897. Absolute reconstruction acceleration falls from 2.31425e-6 to 7.44897e-7 m/s², while the resolved eddy contribution changes from 5.11744e-6 to 3.63547e-6 m/s². The coarse 15- and 30-second results agree at these displayed digits. The final-18-hour ratios are 0.459533 and 0.209114, so the reduction persists beyond initial adjustment.

This is encouraging resolution sensitivity, not proof of convergence or eddy-driven Ferrel cells. The fine-grid correction remains substantial, the resolved contribution also changes, and this is only one day after remapping. Column enthalpy storage remains negative (−22.44 versus −22.24 W/m² over the full day), despite closed numerical thermal accounting. The vertical grid and boundary treatments were not varied.

Unity jobs: preflight 1718914, coarse 1718940, fine 1718941, review/render 1719193. The coarse production run used an A100 and the fine run an H100, so their elapsed times are not a controlled hardware cost benchmark. No raw fields were downloaded. Only verified finished publication assets were transferred for the GitHub upload.

See [full review](resolution_review.json), [media verification](PASS.json), [fixed movie scales](coarse_over_fine_300_301days.json), and [source/methods](../../scripts/ferrel35_resolution_v1/README.md).
