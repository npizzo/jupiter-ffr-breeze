# Complete 100-day heating/control comparison

[Watch all three movies on the gallery](https://npizzo.github.io/jupiter-ffr-breeze/#cascade100).

Both cases completed all ten ten-day segments. Heated final job 66992567 and unheated final job 66992568 exited successfully. Movies include all 801 saved times, from day 0 through day 100 inclusive at three-hour intervals. Each has synchronized columns at 19.6875, 59.0625 and 100.3125 km, with vorticity above passive tracer. The comparison and separate movies share the same fixed scales across both cases and all times: symmetric asinh vorticity (width 0.1 in ζ/f₀), and asinh concentration (width 0.01, physical concentration ticks).

The dry f-plane experiment uses a 6000 × 6000 × 120 km channel, 256 × 256 × 64 cells, Jupiter rotation at 55°N, no beta or imposed shear, and identical 1200-km random eddies. Bottom heating is 162 versus 0 W/m². The stable cap and thermal-feedback model remain unchanged. [Full experiment design](../../docs/cascade100.md).

Scientific computation, rendering and complete media decoding ran on Unity. Job 1431877 rendered the complete package after correcting numeric segment ordering (the initial rendering attempt rejected nonchronological frames). Job 1432183 checked every frame and full decoding. Smaller web copies were encoded by job 1432273 after GitHub rejected the larger upload; dimensions, frames, times and scales were retained. Original movies remain on Unity.

[Provenance and movie hashes](provenance.json), [decoding checks](media-validation.json), and individual movie JSON files record source/configuration hashes, physical times and encoding. Segment metadata are included. This publication presents output, not a new claim of an inverse cascade or numerical convergence. Tracer concentration is not cloud brightness.
