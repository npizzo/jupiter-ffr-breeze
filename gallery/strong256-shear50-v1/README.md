# Higher-resolution shear cases: days 0–50

[Watch all three movies](https://npizzo.github.io/jupiter-ffr-breeze/gallery/strong256-shear50-v1/).

Completed 256 x 256 x 64 simulations of no shear, cyclonic initial shear and anticyclonic initial shear. Same physical conditions as the preceding 128-squared cases: 6000-km box, 1200-km random eddies, 5.5022 m/s initial random RMS, Ro=0.1, dry anelastic f-plane at 55 N, heating 162 W/m2 and stable cap. Each movie combines vorticity over tracer at three depths, all 401 three-hourly states, with fixed color scales shared among these fine-grid cases. Coarse and fine movies have independently calculated vorticity limits.

All three model integrations reached day 50. Their wrapper exited afterward because it assumed the continuation would re-emit day 1. That frame was retained in the original segment; timeline and checkpoint-provenance checks passed without rerunning the simulations. Duplicate restart-frame equality was unavailable, not passed. The replacement renderer verifies all timestamps and full movie decoding. See completion-audit.json and publication.json.

This is a resolution comparison, not a claim of convergence or preferred vortex sign. All simulation, scientific checking, rendering and media validation ran on Unity. Only finished publication assets were transferred for upload.
