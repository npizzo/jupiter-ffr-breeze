# Prescribed cooling and reduced stable-layer thermal mixing

This experiment corrects the forcing and transport limitations found in the first 48-hour controls. The dynamics remain dry anelastic in the same 1,200 × 1,200 × 120-km periodic box. There is no beta or imposed shear. Jupiter's rotation remains fixed at Ω = 1.7585 × 10⁻⁴ s⁻¹ and f = 2.88096 × 10⁻⁴ s⁻¹ at 55°N.

## Upper heat sink

Instead of relaxing temperature toward the initial reference profile, prescribe a volumetric energy loss through Breeze's `ρE` forcing interface. The sink is proportional to sin²(π(z−80 km)/40 km) between 80 and 120 km, zero elsewhere. Its weights are normalized on the actual grid so the column integral is exactly −5.4 W/m². It is independent of temperature and cannot become heating. Each heated run therefore has prescribed external enthalpy powers of +5.4 W/m² below and −5.4 W/m² aloft.

The heating-off pair keeps the upper sink active, isolating the bottom input. These are cooling-only controls, not unforced controls; the extra unforced case disables both inputs. Balanced enthalpy powers do not imply zero change in density-weighted potential temperature because the source conversion contains the height-dependent Exner function. The diagnostics explicitly test the correct potential-temperature inventory, rather than assuming it remains constant. This source accounting is not a proof of a closed discrete total-energy budget.

## Thermal diffusivity

Retain κ = 1,000 m²/s below 60 km, taper smoothly to κ = 1 m²/s between 60 and 80 km using a cubic smoothstep, and keep κ = 1 m²/s throughout the stable cap. The 20-km transition spans about 11 cells at the standard resolution. Momentum viscosity remains 1,000 m²/s so the revision isolates thermal mixing.

This is a prescribed sensitivity experiment, not a calibrated Jovian turbulence closure. The cap value prevents the former artificial thermal diffusion from dominating the intended heat input. Low diffusivity does not make subgrid-scale boundary layers resolved; the vertically refined run checks sensitivity but is not a full convergence study. The same explicit and numerical diffusion caveats apply as before.

## Cases

All cases run for 48 Earth hours. The four primary cases repeat the previous seed, unstable lower-layer reference profile, and perturbation parameters.

| Case | Bottom heating [W/m²] | Upper cooling [W/m²] | Perturbation parameter [K] | Grid |
|---|---:|---:|---:|---|
| controlled48_on_p020 | 5.4 | 5.4 | 0.2 | 128 × 128 × 64 |
| controlled48_off_p020 | 0 | 5.4 | 0.2 | 128 × 128 × 64 |
| controlled48_on_p002 | 5.4 | 5.4 | 0.02 | 128 × 128 × 64 |
| controlled48_off_p002 | 0 | 5.4 | 0.02 | 128 × 128 × 64 |
| controlled48_unforced | 0 | 0 | 0.02 | 128 × 128 × 64 |
| controlled48_refined | 5.4 | 5.4 | 0.02 | 128 × 128 × 128 |

The initial random field is not identical across resolutions: the same seed on a different grid produces a different spatial realization. Refined comparisons therefore concern bulk behavior, source balances, and cap preservation, not pointwise trajectory agreement. The same-resolution heated/unheated pairs do have identical initial perturbations.

## Checks and diagnostics

Before the six runs, short hydrostatic-rest, bottom-heating-only, fixed-cooling-only, and balanced-heating/cooling checks validate the source sign, column-integrated sink, and expected potential-temperature change. The job stops if a check fails. During every revised run the upper input must remain nonpositive; its integrated value is checked against the configured loss. Finite fields and potential-temperature inventory checks apply to every case.

New diagnostics record the actual upper power, net prescribed enthalpy input, stable-cap temperature-contrast retention, and mean-profile departure from the initial cap. Output profiles include the prescribed κ and volumetric sink. `relaxation_*` diagnostic columns are retained as compatibility aliases for the total upper source; canonical `upper_*` columns identify the fixed sink in these runs.

Submit `sbatch slurm/controlled48.sbatch`. One L40S GPU executes all cases sequentially in one Julia process. The allocation requests four CPUs, 8 GB host memory, and a two-hour ceiling, including the larger vertical-resolution case. Outputs use a new directory, `data/controlled48-v2`; the previous experiments remain intact.

The first 48 hours remain a diagnostic window. We will assess whether the corrected source stays controlled, the cap survives, and bulk convection depends on imposed forcing after its initial transient. A persistent flow alone would not establish FFR morphology or thermal equilibrium.

The first allocation received an A16 GPU and could not complete the campaign within its original 45-minute limit. Its partial output is preserved under `data/controlled48-v1`; the replacement uses identical physical parameters and a new output directory. The four preliminary validation cases passed before the interrupted first long case. Full campaign results are pending.
