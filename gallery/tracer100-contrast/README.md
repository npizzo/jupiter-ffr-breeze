# Full 100-day movies with improved tracer contrast

[Watch on the gallery](https://npizzo.github.io/jupiter-ffr-breeze/#tracer100-contrast).

The user requested a clearer view of late-time passive-tracer gradients. These variants use unchanged physical concentration C with a fixed brown–teal color scale centered at C = 0.0793082 and asinh width 0.001. The mapping is shared across both cases, all depths and all times; no per-frame normalization or clipping is used. Vorticity scaling is unchanged. Original movies are retained.

[Rendering method and audit](../../docs/tracer100-contrast.md). All three movies were rendered and fully decoded on Unity by job array 1437852. Each contains all 801 saved times, from day 0 through day 100, at 12 frames/s. Vorticity normalization was checked against the original movies. The audit ran as job 1437416. [Publication provenance and hashes](provenance.json) and [contrast audit](audit.json) accompany the movies. The gallery provides a button to start at day 75.
