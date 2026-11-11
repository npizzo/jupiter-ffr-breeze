# Pressure-diagnostic initialization revision

The [physical design and controls](ferrel2026.md) are unchanged. Current Unity
jobs are unit 1403142, CPU closure/restart 1403143, GPU quick 1403144,
full-grid pilot 1403145, fine pilot 1403154 and 48-hour cases 1403155–1403157.
The [launch record](new-campaigns-launch.json) records validation status.

The first GPU diagnostic gate stopped before time stepping because its native zonal pressure-gradient mean was `0.0005639795803662246`, above the unchanged `1e-10` tolerance. The periodic zonal pressure gradient must telescope to zero when the reference density depends only on height.

The pinned Breeze source fills pressure-anomaly halos after each Poisson solve. Its ordinary state update refreshes prognostic and velocity halos; the separately stored pressure anomaly is not a prognostic field. Our direct initial pressure assignment preceded the first Poisson solve. An unrefreshed periodic pressure halo therefore contaminates the initial diagnostic stencil. Revision 2 records the actual native pressure-gradient mean before and after the explicit halo fill so Unity can verify this diagnosis rather than assume it.

Only `scripts/ferrel2026_v2/` and `slurm/ferrel2026-v2*.sbatch` implement this revision. The original source and independent original CPU probe remain unchanged. Revision 2 reuses the immutable `configs/ferrel2026/` physics settings in fresh output roots. It fills pressure halos after initial pressure assignment and before pressure diagnostics, retaining the original zero-mean assertion. Physical interior values, governing operators and scientific parameters are unchanged.

The new chain requires a CPU native-closure/restart smoke check before the GPU quick gate:

1. Initialization and source-hash checks.
2. Unity CPU smoke: 64 × 64 × 24, 6,000-km random seed, 120 seconds at steps no larger than 5 seconds; continuous versus restart comparison, every-step native SSPRK3 closure and both axisymmetric controls.
3. The existing 96² GPU diagnostic/restart gate, 512² screening, 1024² early comparison and bounded 48-hour experiments.

The CPU smoke is an API and diagnostic check, not a physical-resolution validation. Production remains blocked unless every prerequisite succeeds. All numerical execution stays on Unity.

Pressure anomaly is not part of Breeze's prognostic checkpoint. Immediately
after restart, the pressure diagnostic therefore uses initialized pressure until
the first pressure solve. It tests zonal cancellation only; it is not a claim
that the previous physical pressure field was restored. The independent measured
momentum-increment check still tests the effects of the actual pressure projection.

Unity audit 1403720 passed: all original physical-field restart comparisons passed;
the maximum streamfunction difference was 3.252e-10 kg/(m s), below the measured
mass-flux integration bound of 2.024e-9 kg/(m s). The difference from the integrated
native mass-flux change was 4.337e-19 kg/(m s), within the 2.414e-15 roundoff bound.
This supports a dimensionally appropriate derived-field comparison while retaining
all original physical-field tolerances. The failed chain remains preserved.
