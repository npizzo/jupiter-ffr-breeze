# GPU callback compatibility fix: unchanged SI physics

The repaired CUDA environment successfully passed runtime validation and CPU
checks. The subsequent GPU test failed while constructing a QG model: the
planetary-PV callback captured the string-valued geometry setting. A GPU kernel
cannot receive that non-isbits captured string.

This isolated v3 resolves the geometry name into Boolean flags on the CPU,
then constructs the callback using only those flags and the floating-point trap
radius. The expressions for the polar, flat and zero PV geometries are unchanged.
A fail-fast isbits guard and explicit none/polar/flat construction tests prevent
this deployment failure from being mistaken for a scientific failure.

All SI parameters, initial conditions, time integration, energy limits, output
settings and the eight-case grid/seed comparison remain those of the frozen v2
configuration files. The new scripts live in scripts/gamma2026_v3; all v2 files
remain unchanged. The already-verified gamma2026_gpu environment and isolated
runtime depot are reused without modification.

The entire CPU gate and both GPU gates rerun against the new source hash before
any dependent production run. No old validation hash is reused. The first CPU
job needs no dependency on the completed runtime-setup job.

```sh
bash scripts/gamma2026_v3/submit.sh data/gamma2026-v3
```

All execution and diagnostics stay on Unity. This is a callback compatibility
repair, not a change to the physical model or an assertion of reproduction.
