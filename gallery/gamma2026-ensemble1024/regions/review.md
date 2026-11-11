# Interior organization survives excluding the trap ring

The polar background produces a consistent interior contrast across all three completed 1024² seeds. Removing the outer 4–5 Lγ annulus **from the measurement** preserves the difference from the flat control. This strengthens the evidence for interior sign organization; it does not establish that the dynamics is independent of the trap boundary.

These are unforced, two-dimensional barotropic QG runs through 75 Earth days, without convection. [Watch the synchronized polar/control movies in the gallery](https://npizzo.github.io/jupiter-ffr-breeze/#gamma75).

The statistic S is the negative-vorticity mean radius minus the positive-vorticity mean radius, weighted by each sign's vorticity magnitude inside the stated region. Positive S means positive vorticity is more centrally weighted **within that region**. The table reports the already-computed polar-minus-flat difference in the change of S from day 0 to day 75, in Lγ:

| Seed | Interior r < 3 Lγ | Interior r < 4 Lγ | Whole trap r < 5 Lγ |
|---|---:|---:|---:|
| 11 | +0.5953 | +0.5087 | +0.5180 |
| 22 | +0.5353 | +0.7187 | +0.6118 |
| 33 | +0.3954 | +0.5084 | +0.4551 |

All three prespecified regions retain a positive paired contrast. The **absolute** day-75 polar separation is also positive in both interiors: S = 0.504–0.597 Lγ inside r < 3 Lγ and 0.456–0.654 Lγ inside r < 4 Lγ. The flat controls are close to zero there: −0.004–0.109 Lγ and −0.064–0.020 Lγ, respectively. Including the outer annulus reverses the absolute sign in the polar cases, explaining why the earlier whole-trap statistic was negative despite the visible interior organization.

The ring remains quantitatively important. At day 75, the 4–5 Lγ annulus contains 59–62% of the domain's positive circulation and 42–53% of its relative enstrophy in the polar cases; the corresponding control ranges are 35–37% and 2.6–5.2%. These are fractions of the remaining domain totals, not growth rates or energy fractions. In the innermost r < Lγ region, the positive-circulation fraction is larger in every polar case than its matched control, but the radial-separation sign there varies between seeds. The results do not establish a universal ordering within that smallest region.

The allocated Unity observer processed all 76 daily frames in each case. Its analytical fixtures, regional partition/moment identities, and saved-field energy/relative-enstrophy checks passed. The existing total-PV reconstruction discrepancies remain unresolved at day 1 in flat seeds 11 and 33; the report therefore retains its explicit diagnostic-discrepancy status.

Fixed-region moments can change through entry and exit, changes in vorticity strength, and filtering, as well as rearrangement. They do not track individual vortices or measure an inverse cascade. Ring exclusion removes its direct statistical contribution, not its possible dynamical influence. The next comparison remains the already planned 2048² matched refinement; no new model target is selected from this diagnostic alone.

Source: completed Unity regional report, job 1423043; analytical-only check job 1423042. This interpretation uses the existing JSON summaries; no raw-field calculations were performed locally.
