# Completed long-run publication — November 14, 2026

[Watch the completed runs](https://npizzo.github.io/jupiter-ffr-breeze/#completed-nov14).

| Experiment | Saved interval (Earth days) | Saved frames | Simulation / rendering |
|---|---:|---:|---|
| 35°N heated beta-plane | 0–300 | 2,401 | Final segment 1459497; original renderer 1459498; contrast renderer 1707676 |
| 55°N unforced above heated beta-plane | 0–100 | 801 per case | Unforced final segment 1650252; renderer 1650253 |
| 35°N diagnostic restart | 200–210 | 81 | Simulation 1694155; renderer 1694156 |

All science computation, analysis, rendering and full-video decoding were performed on Unity. Only finished publication assets were transferred for GitHub browser upload. No model fields were downloaded. The running wider-channel experiment was not modified.

All three primary movies use three simultaneous height columns (19.6875, 59.0625 and 100.3125 km above the bottom), vorticity followed by passive-tracer concentration. The 55°N comparison has four rows: unforced vorticity/tracer above heated vorticity/tracer. Each movie contains all saved three-hourly states, played at 12 frames/s, with no temporal interpolation.

Vorticity uses fixed symmetric asinh normalization at each height (linear width 0.1 in ζ/f₀), spanning each movie's full extrema. The 55°N pair shares limits between cases. Positive vorticity is cyclonic in the northern hemisphere. Primary tracer movies use the established BrBG contrast, centered at C=0.0793082 with two-sided asinh width 0.001 and physical concentration range 0–1. The mapping stays fixed across time, depths and matched cases. The original 300-day tracer-scale movie is retained as a linked alternative.

The 55°N unforced case has zero boundary heat flux, no cooling and no sponge. Its saved initial fields match the heated case. This compares the full forced setup with free decay, rather than isolating only bottom heating.

The 35°N restart passed the native momentum-budget numerical gate over days 200–210. Its native-minus-centered advection remainder is retained explicitly: the centered three-hourly mean/eddy diagnostic is not an exact native operator partition. Numerical closure does not establish Ferrel cells. The 300-day zonal-flow history and its diagnostic record are also included, without claiming equilibrium or an inverse cascade from the movie.

Public text metadata replaces the absolute Unity repository prefix with relative paths. Original source hashes remain in source records; the publication manifest gives checksums of the actual published files. Existing gallery movies and concurrent repository changes are preserved.
