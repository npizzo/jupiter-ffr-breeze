# Preliminary Unity validation

These tests establish basic numerical operation; they do not establish FFR morphology or superiority to rotating RBC.

## CPU results

- Hydrostatic rest: maximum vertical velocity in the saved vertical section was 1.6 × 10⁻²⁹ m/s after 60 s.
- Heating-only column budget: the density-weighted potential-temperature increase matched the imposed bottom flux with relative error 3.0 × 10⁻⁸ after 60 s.
- Seeded convection smoke test: completed 60 s with finite fields; volume RMS vertical velocity reached 0.0186 m/s.
- The smoke-test reference density changes by a factor of approximately 9.24 across the cell centers. Reference pressure spans approximately 0.93–19.13 bar; the lower boundary is at 20 bar.

[Machine-readable diagnostic values](validation.json). The velocity maximum above is a saved vertical section, not a whole-domain maximum.

Preparation plus the first CPU smoke test took 9 min 23 s and peaked at about 4.2 GiB host memory, including installation and compilation. These numbers describe environment preparation, not a solver benchmark. The smoke-test measured run interval was 12.7 s, including output and first-step compilation but excluding process startup and model construction.

## GPU pilot

A Float64 128 × 128 × 64 run completed four simulated hours in 733 steps on an NVIDIA L40S. The short seeded CPU/GPU comparison passed; relative L2 differences were below 7 × 10⁻¹² for the saved velocity and vorticity fields, and about 1.3 × 10⁻¹⁵ for potential temperature.

The pilot developed convective columns beneath the stable cap. Peak whole-domain vertical speed among logged snapshots was 20.79 m/s; volume RMS vertical velocity peaked near 3.12 m/s and ended at 2.80 m/s. The final horizontal vorticity slice reached |ζ/f| ≈ 1.25. These snapshots do not show convincing folded-filament morphology. Initial thermal perturbations and weak superadiabaticity provide energy in this transient; sustained convection driven by the imposed bottom flux has not yet been demonstrated.

### Cost and limits

- Entire GPU batch, including first-time CUDA installation, compilation, smoke test, agreement checks, and pilot: **9 min 53 s**.
- Pilot `wall_seconds`: **36.01 s**, including simulation initialization, first-step compilation, and output; excluding process startup and model construction.
- Interval after the first nonzero-time output: **13.21 s** for the remaining 13,800 simulated seconds. This is a throughput observation, not a fair benchmark against another code.
- Batch maximum host RSS: approximately **5.64 GiB**.
- Slurm's sampled GPU-memory maximum: **776 MiB**. Sparse sampling may miss short peaks; this is not a guaranteed VRAM requirement. Direct high-frequency measurement is needed before aggressive GPU downsizing.

The initial GPU request was one GPU, 16 GiB host memory, and one hour. Future pilot scripts request 8 GiB and 30 minutes based on these measurements. Runtime is short enough to make controlled parameter studies practical. Longer integration, larger domains, imposed-shear controls, stronger-rotation regimes, and resolution tests are needed before assessing FFR realism or comparing cost against RBC.

[Gallery](../gallery/pilot/README.md) · [Run parameters](../gallery/pilot/metadata.json) · [Diagnostic time series and resource measurements](pilot-results.json).
