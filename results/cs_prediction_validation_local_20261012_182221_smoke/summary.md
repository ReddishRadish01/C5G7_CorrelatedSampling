# CS prediction validation against direct Monte Carlo

- Reference depth: 83.6388953145 cm
- Histories per cycle: 20000
- Cycles: 10 inactive + 12 active
- CS source propagation: L10
- Independent odd seeds: direct reference 92239000101; every target and CS run uses another seed
- Overall: 0 / 1 probe steps passed all enabled checks

The primary residual is `CS delta-k - [direct k(target) - direct k(reference)]`. Its z-score uses the quadrature sum of the CS delta-k standard error and both independent direct-k standard errors.

| step (cm) | direct delta-k | CS delta-k | agreement z | direct sign | source tail | status |
|---:|---:|---:|---:|---|---|---|
| 0.125 | -0.000384500 ± 0.004450374 | -0.000658330 ± 0.000394910 | -0.06 | UNRESOLVED | fail | SOURCE_TAIL_NOT_CONVERGED |

Interpretation: `UNRESOLVED` direct sign means that two independent direct runs cannot distinguish the small rod-worth signal from Monte Carlo noise at the selected confidence level. A positive observed direct delta-k in that case is not evidence that insertion physically increased reactivity. Increase histories or active cycles before interpreting the sign.

Rows share the same direct reference calculation. Per-row z-scores are valid, but the row residuals are correlated and must not be treated as independent samples in a naive global chi-square test.
