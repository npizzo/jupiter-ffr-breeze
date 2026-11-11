# Dimension-aware streamfunction restart validation

**Unity validation passed:** initialization job 1405388 passed all 104 checks,
including deliberate corruption rejection. CPU job 1405389 passed native
momentum/mass closure, continuous/restarted agreement, both axisymmetric controls
and tracer checks. GPU quick check 1405390 and full-grid pilot 1405391 precede
fine pilot 1405392 and the three 48-hour cases 1405393–1405395.
See the [physical experiment design](ferrel2026.md) and
[launch record](new-campaigns-launch.json). Numerical readiness is not a claim
that steady or deep Jovian Ferrel cells have been reproduced.

Revision 3 changes only the validation of the derived mass streamfunction. Governing equations, parameters, native operators and all twelve original physical/native-field restart criteria remain unchanged (`rtol=1e-9`, `atol=1e-10`). The pressure-halo initialization fix from revision 2 is retained. Original revisions and output remain untouched.

Unity audit job **1403720** passed after reading the revision-2 continuous and restarted output. All twelve physical/native criteria passed. The measured maximum velocity difference was `8.53e-14 m/s`; the maximum native meridional mass-flux difference was `1.68667e-14 kg/(m² s)`.

The streamfunction has units `kg/(m s)` and integrates those mass fluxes over the 120-km depth. Its maximum difference was `3.25163e-10 kg/(m s)` and its array-norm difference was `2.58990e-9 kg/(m s)` (relative norm `2.43018e-8`). Applying the same absolute norm tolerance as the physical velocity/mass-flux fields was dimensionally inappropriate.

The audit verified the actual relationship

\[
\delta\Psi_{j,k}=-\Delta z\sum_{\ell<k}\delta\overline{\rho v}_{j,\ell},
\qquad
\|\delta\Psi\|_\infty\le H\|\delta\overline{\rho v}\|_\infty+\epsilon_{\rm sum}.
\]

Here the measured integral bound was `2.02401e-9 kg/(m s)`. The explicit Float64 summation-roundoff envelope was `2.41407e-15 kg/(m s)`, and the discrepancy between the observed streamfunction difference and the integrated mass-flux difference was only `4.33681e-19 kg/(m s)`. Both individually saved streamfunctions also reconstructed from their own native mass fluxes within this roundoff envelope.

The shared revision-3 helper therefore checks:

- The observed streamfunction difference equals the integral of the measured native mass-flux difference within the Float64 roundoff envelope.
- Both the whole-depth maximum bound and the sharper, pointwise accumulated absolute-mass-flux bound hold.
- Each saved streamfunction reconstructs individually from its own native mass flux, with a zero bottom boundary value.

The roundoff envelope is `8*(Nz+2)*eps(Float64)*S`, where `S` is the largest absolute streamfunction magnitude or vertically integrated absolute native-mass-flux throughput across both outputs. Every factor has the appropriate streamfunction units. It is not an arbitrary new absolute tolerance. The physical-field tests remain independent, so this criterion cannot accept an incorrect restart merely because its streamfunction is internally consistent.

Both CPU and GPU restart comparisons call the same helper. The Unity unit stage first confirms acceptance of an identical synthetic streamfunction, then deliberately changes one interior streamfunction value while keeping its native mass flux identical. The exact acceptance predicate used by the production checker must reject that corruption on its integral identity, whole-depth bound and individual reconstruction checks. The isolated chain remains initialization → CPU native-closure/restart smoke → GPU quick gate → 512² screening → bounded 48-hour experiments and the independent 1024² check. It reuses unchanged `configs/ferrel2026/` settings and requires fresh output roots. All numerical validation runs on Unity.
