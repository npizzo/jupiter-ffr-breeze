# Native-face diagnostic: Earth days 300–301

[Watch the three-depth vorticity/tracer movie and read the budget review](https://npizzo.github.io/jupiter-ffr-breeze/#ferrel35-face301).

Completed on Unity with unchanged 35°N heated beta-plane physics. The original-driver comparison passed 23 checks (19 physical-field comparisons); all 2,880 production timesteps were audited. The native momentum and potential-temperature budgets close. The WENO-minus-linear-product reconstruction contribution is 45.2% of the resolved eddy acceleration RMS in the stated interior, so numerical transport remains significant. This is a diagnostic validation, not a Ferrel-cell identification or resolution-convergence result.

The movie contains all nine three-hourly states, including the original day-300 checkpoint state, at one frame per second. It preserves three synchronized depths, vorticity above passive tracer, with fixed scales and no interpolated frames. It was rendered and fully decoded on Unity.

- `review.json`: measured budgets, region, checks and limitations.
- `publication.json`: source revision, Unity jobs, GPU and restart provenance.
- `PASS.json`: full media-validation receipt and asset hashes.
- `diagnostic_300_301days.json`: exact movie times and color scales.
- [Method and source](../../scripts/ferrel35_face_v3/README.md).
- [Configuration](../../configs/ferrel35-face-v3/diagnostic.toml).

The column enthalpy budget uses cell-center reference Exner weights and is not a full total-energy budget. Its one-day storage rate is −22.44 W/m², despite the imposed external-heating feedback. Wall, cap and resolution sensitivity remain to be established.
