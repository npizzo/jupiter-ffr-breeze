# Stronger-heating snapshots

All six runs completed. Bottom heating now drives vigorous convection and strongly changes the seeded rotating flow. Fine-scale vorticity remains resolution-sensitive; folded filamentary regions are not established.

[Results and interpretation](../../docs/strong48-results.md) · [Experiment](../../docs/strong48.md) · [Metrics](summary.json) · [Resources](resources.json)

## Evolution and control comparisons

Solid curves include bottom heating and upper cooling; dashed curves have identical upper cooling and initial turbulence, with no bottom heating.

![Energy and vertical speed](evolution.png)

## 100× heating, finer horizontal grid

Vertical sections show vertical velocity and potential-temperature anomaly at 6, 12, 24 and 48 hours. Each row uses one color scale. The line marks the original 80-km cap base.

![Fine-grid vertical snapshots](f100_h256_heated_snapshots.png)

The following maps are vertical vorticity divided by f near 20 km. Heated and control maps share the same scale; a faint control means its amplitude is much smaller, not zero.

![Fine-grid final comparison](f100_h256_final_vorticity.png)

![Fine-grid vorticity evolution](f100_h256_heated_vorticity.png)

![Fine-grid final vertical comparison](f100_h256_final_sections.png)

## 30× heating

![30× vertical snapshots](f030_heated_snapshots.png)

![30× vorticity evolution](f030_heated_vorticity.png)

![30× final vorticity comparison](f030_final_vorticity.png)

![30× final vertical comparison](f030_final_sections.png)

## 100× heating, original horizontal grid

![100× vertical snapshots](f100_heated_snapshots.png)

![100× vorticity evolution](f100_heated_vorticity.png)

![100× final vorticity comparison](f100_final_vorticity.png)

![100× final vertical comparison](f100_final_sections.png)

## Spectra and data

These are horizontal velocity spectra from one level near 20 km, not full-depth transfer budgets. The same horizontal shell index represents the same physical wavelength across these fixed-size domains.

![Final horizontal spectra](spectra.png)

Full five-minute diagnostic histories:

- [strong48_f030_heated](strong48_f030_heated_diagnostics.csv.gz)
- [strong48_f030_cooled](strong48_f030_cooled_diagnostics.csv.gz)
- [strong48_f100_heated](strong48_f100_heated_diagnostics.csv.gz)
- [strong48_f100_cooled](strong48_f100_cooled_diagnostics.csv.gz)
- [strong48_f100_heated_h256](strong48_f100_heated_h256_diagnostics.csv.gz)
- [strong48_f100_cooled_h256](strong48_f100_cooled_h256_diagnostics.csv.gz)
