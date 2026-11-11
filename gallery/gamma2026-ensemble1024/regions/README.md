# Fixed-region gamma-plane diagnostics

Completed two-dimensional, unforced barotropic QG cases; all76 daily fields are included. These fixed regions describe actual signed vorticity and relative enstrophy, not tracked vortices, convection or inverse energy transfer.

Masks: r/Lγ <1, 1–3, 3–4, 4–5 and ≥5; cumulative interiors r/Lγ <3, <4 and <5. Each annulus includes its lower boundary and excludes its upper boundary. The exterior is restricted to the periodic computational square. Fractions use whole-domain denominators of the same sign/quantity. Undefined moments are blank/null.

Native energy and relative-enstrophy checks retain the frozen publication thresholds. **Unresolved total-PV diagnostic failures remain**; see regions.json and native_audit. No tight PV-conservation claim is made.

The figure shows every selected case; CSVs contain all regions/times, JSON records initial/final values, matched changes and provenance. Positive paired change means polar minus flat change, not necessarily positive absolute separation. Changing contents of fixed masks does not measure individual vortex migration.

[Watch the already verified gamma comparison movies](https://npizzo.github.io/jupiter-ffr-breeze/#gamma75). This observer does not generate new movies.
