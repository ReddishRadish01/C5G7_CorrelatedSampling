# CS prediction validation against direct Monte Carlo

- Reference depth: 83.6388953145 cm
- Histories per cycle: 500000
- Cycles: 100 inactive + 150 active
- CS source propagation: L10
- Independent odd seeds: direct reference 92239000101; every target and CS run uses another seed
- Overall: 3 / 5 probe steps passed all enabled checks

The primary residual is `CS delta-k - [direct k(target) - direct k(reference)]`. Its z-score uses the quadrature sum of the CS delta-k standard error and both independent direct-k standard errors.

| step (cm) | direct delta-k | CS delta-k | agreement z | direct sign | source tail | status |
|---:|---:|---:|---:|---|---|---|
| 0.03125 | -0.000228500 ± 0.000173768 | -0.000193058 ± 0.000014508 | 0.20 | UNRESOLVED | pass | PASS |
| 0.0625 | -0.000251500 ± 0.000162953 | -0.000396194 ± 0.000015379 | -0.88 | UNRESOLVED | fail | SOURCE_TAIL_NOT_CONVERGED |
| 0.125 | -0.000862600 ± 0.000150490 | -0.000793810 ± 0.000018576 | 0.45 | RESOLVED_EXPECTED | pass | PASS |
| 0.25 | -0.001388800 ± 0.000183023 | -0.001520002 ± 0.000065067 | -0.68 | RESOLVED_EXPECTED | fail | SOURCE_TAIL_NOT_CONVERGED |
| 0.5 | -0.003155800 ± 0.000161680 | -0.003107273 ± 0.000048232 | 0.29 | RESOLVED_EXPECTED | pass | PASS |

Interpretation: `UNRESOLVED` direct sign means that two independent direct runs cannot distinguish the small rod-worth signal from Monte Carlo noise at the selected confidence level. A positive observed direct delta-k in that case is not evidence that insertion physically increased reactivity. Increase histories or active cycles before interpreting the sign.

Rows share the same direct reference calculation. Per-row z-scores are valid, but the row residuals are correlated and must not be treated as independent samples in a naive global chi-square test.
