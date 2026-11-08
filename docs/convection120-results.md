# Five-day 30× convection: completed results

All four 120-hour runs completed and passed their output, forcing-budget, initial-amplitude and thermal-inventory checks. Bottom heating sustains vigorous convection through day five while the matching cooling-only flows decay. Increasing vertical resolution changes final-day bulk measures by only 2–4%. Horizontal vorticity is still sensitive to horizontal resolution, and a large folded-filament system is not established.

[Videos and snapshots](../gallery/convection120/README.md) · [Experiment](convection120.md) · [Metrics](../gallery/convection120/summary.json)

## Final-day comparison

Means below cover **96–120 Earth hours**. Vertical velocity and vertical-vorticity RMS are volume averages over the entire domain; column kinetic energy is density-weighted. The final resolved heat-flux maximum is a maximum over height at 120 hours, not an exact boundary-flux or total-energy closure.

| Case | Vertical RMS [m/s] | Column KE [J/m²] | ζ RMS / f | Final maximum upward resolved heat flux [W/m²] |
|---|---:|---:|---:|---:|
| 64 levels, heated | 1.228 | 194089 | 0.590 | 127.682 |
| 64 levels, cooling only | 0.041 | 10579 | 0.033 | 0.005 |
| 128 levels, heated | 1.253 | 201563 | 0.602 | 125.429 |
| 128 levels, cooling only | 0.039 | 10355 | 0.033 | 0.000 |

Heated runs have about 18–19 times the final-day kinetic energy and 30–32 times the vertical RMS speed of their controls. Their diagnosed buoyancy work remains about 21–22 W/m², compared with approximately zero in the controls. Relative to hours 24–48, heated-run kinetic energy rises only 1.55% (64 levels) and 1.81% (128 levels); vertical RMS rises 2.79% and 2.03%. These nearly flat bulk histories support sustained convection over this window, rather than a decaying initial burst. They do not establish statistical stationarity or seed independence.

## Resolution and rotation

At fixed 256² horizontal resolution, doubling vertical levels from 64 to 128 changes final-day vertical RMS by **+2.08%**, kinetic energy by **+3.85%**, and vertical-vorticity RMS by **+2.09%**. The last snapshot's maximum upward resolved heat flux changes from 127.68 to 125.43 W/m². Bulk convection is relatively insensitive to this vertical refinement.

The horizontal test remains less reassuring. Comparing hours 24–48 at 128² and 256² with 64 vertical levels, horizontal refinement changes vertical RMS by **+3.70%**, energy by **+6.36%**, and vorticity RMS by **+51.65%**. Thus vertical robustness does not imply horizontal convergence of morphology. Initial mode sampling and thermal noise differ with grid resolution; this is a sensitivity study with one seed and two grids, not a convergence proof.

Jovian f, beta = 0 and zero imposed shear remain unchanged. Initial eddy Rossby is 0.1; evolved volume ζ RMS / f is about 0.59–0.60 on day five, with larger local values. Rotation remains present, but an asymptotically small-Rossby or geostrophic-turbulence regime is not demonstrated by these runs.

## Spatial organization, heat transport and cap

The movies and stills show many small vortices, short elongated features, and vertical convective structures beneath the stable cap. They do not show a clear large-scale folded-filament system. At 120 hours, horizontal modes 1–2 contain only about 0.9–1.4% of the heated near-20-km slice's non-mean horizontal kinetic energy. These are single-height spectra, not full-depth statistics or transfer budgets; they cannot establish the presence or absence of an inverse cascade.

Resolved upward transport peaks near the bottom, then declines through the layer; a large maximum alone does not imply an equilibrated through-layer heat flux. The mean potential-temperature contrast between 82.5 and 117.5 km remains within 1% of its initial value in all cases. Cap profile RMS changes reach about 0.61–0.62 K, so survival of broad stratification does not mean the cap is unchanged.

## Videos, validation and cost

Two paired MP4s show every hourly frame from 0 through 120 hours, at six frames per video second. Each pair uses one fixed symmetric ζ/f scale over the entire movie, with actual timestamps and map heights. A weak-looking control on that shared scale still has nonzero motion. Initial overshoot sets a larger color range than the later flow. The maps are dry-model vorticity, not clouds or JIRAM radiance.

Public videos are compact viewing copies; original 1200-pixel renders and hourly numerical fields are retained on Unity. Use the fields, not compressed video pixels, for quantitative morphology measurements. The movie metadata records the encoding and color limits. All 121 frames were checked for timing, finite values and consistent paired grids, and both encoded videos decode without errors.

All cases completed at 120 hours, with finite diagnostics and prescribed source powers verified. Saved initial maps, energy profiles and mean temperature profiles match within each heated/control pair. Final relative potential-temperature inventory residuals are below 4 × 10⁻¹⁰; this does not constitute a full discrete total-energy closure. The four jobs used **3 hours 32 minutes 12 seconds** in total, sequentially on one L40S, including startup and validation. Peak reported host memory was about 3.0 GiB. See [resources](../gallery/convection120/resources.json).

## Implication for the next step

The current setup can sustain convection-driven rotating flow at practical cost. FFR morphology remains the open question. The next numerical check should target horizontal resolution at 30× forcing, comparing the same time window, before treating fine vorticity patterns as physical. Any subsequent beta or cyclonic-shear experiment should retain a matched no-beta/no-shear reference. No additional jobs were submitted with this report.
