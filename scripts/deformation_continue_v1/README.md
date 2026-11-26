# Matched deformation pair: continuation to 100 turnover times

Authorized November 26, 2026: continue the two inexpensive unforced anelastic cases to 100 initial turnover times, then to 200 if needed. These are turnover counts, not the native model time: tau=L/U=100, so targets correspond to native t=10000 and 20000.

The model, stratification, rotation, viscosity, resolution and random realization are unchanged. There is no heating, cooling, convection, beta or jet. Both start from their own completed 25-turnover checkpoint.

## Release and cost controls

- The CPU preflight executes the actual continuation driver, including output readback and serialized restart equivalence for both stratifications. The full Julia sources are parsed before execution.
- A receipt binds all executable source, batch/configuration files, Julia environment, Breeze source and the two immutable input checkpoints.
- The production GPU job repeats restart equivalence before loading its full-resolution checkpoint.
- A 25-to-30-turnover segment measures continuation cost. The remaining 30-to-100 segment is released only if its projected runtime is under 40 minutes; each entire GPU allocation is limited to one hour. This is a cost ceiling, not an expected duration. The first campaigns took about 3 minutes of integration for 25 turnovers.
- Existing finite-state, positive-temperature, Mach, anelastic constraint and scalar-conservation guards remain. Scalar conservation is referenced separately to the start of each restart segment and explicitly labeled as such.
- Same output cadence: diagnostic/three-depth frames every 0.1 turnover, full volumes every 5, rolling checkpoints every 2.5. A final explicit output call repairs the old omitted-endpoint issue without changing stepping or dynamics.
- Submission records prevent duplicate releases. No changes to other running experiments.

## Assessment before 200

No automatic physical parameter changes. At 100, assess both completed cases together: vortex mergers and persistence, upper-plane sign share/skewness and three-depth structure, kinetic and potential-energy histories, numerical checks and measured runtime. If the flow is still reorganizing and the extension remains inexpensive, continue the matched pair to 200 using their own checkpoints. Do not treat a cyclonic result as a failure that justifies indefinite extension. No extension past 200 is authorized here.

The v1 driver executes through 100 only. The 200-turnover decision and any subsequent submission must be recorded after the 100-turnover output is examined; this directory does not imply an unattended extension has been scheduled.
