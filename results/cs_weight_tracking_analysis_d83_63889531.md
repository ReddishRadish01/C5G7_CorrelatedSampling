# CS weight tracking report

## Run conditions

- Reference rod depth: 83.6388953145 cm
- Perturbed rod depth: 83.7638953145 cm
- Source propagation level: L=10
- Particles per cycle: 10,000,000
- Cycles: 100 inactive + 200 active
- Active fission-score observations: 16,636,781,080
- Histogram sample observations: 259,966,102
- Invalid observations: 0

## Definitions

- `w`: reference physical particle weight.
- `w*_NPS = w exp(log likelihood)`: transport likelihood-ratio weight.
- `w*R_L`: source-aware perturbed particle weight at L=10.
- `w*_f = w*_NPS (nuSigmaF*/nuSigmaF)`: fission-score adjusted weight.
- Maxima and comparison counters use every observation. Percentiles use the deterministic 1/64 histogram sample.

## Active-cycle results

| Metric | p50 | p99 | p99.9 | p99.99 | Median cycle max | p95 cycle max | Exact max | Max / w | Fraction > w |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| w | 0.843991 | 1.00115 | 1.00115 | 1.00115 | 1.00001 | 1.00064 | 1.00115 | n.a. | n.a. |
| wstar_nps | 0.843515 | 1.09105 | 1.12033 | 2.07773 | 48.4585 | 124.921 | 221.981 | 311.997 | 0.28156% |
| wstar_source | 0.833612 | 1.11682 | 3.14109 | 6.02952 | 73.8373 | 194.367 | 314.54 | 314.374 | 81.11893% |
| wstar_fission | 0.843515 | 1.09105 | 1.12033 | 2.07773 | 48.4585 | 124.921 | 221.981 | 311.997 | n.a. |
| wstar_score | 0.833612 | 1.11682 | 3.14109 | 6.02952 | 73.8373 | 194.367 | 314.54 | 314.374 | 81.11893% |

## Source-aware upper-tail checks

| Threshold | Fraction of sampled fission-score observations above threshold |
|---:|---:|
| 1 | 9.234666% |
| 2 | 0.213912% |
| 3 | 0.111584% |
| 5 | 0.017671% |
| 6 | 0.010213% |
| 10 | 0.003839% |

Approximate diagnostic splitting load (histogram-bin approximation; not a transport-runtime prediction):

| Upper bound | Split cap | Trigger fraction | Extra children per observed score event |
|---:|---:|---:|---:|
| 3 | 4 | 0.098306% | 0.00110706 |
| 3 | 8 | 0.098306% | 0.00115098 |
| 6 | 4 | 0.008030% | 0.00011009 |
| 6 | 8 | 0.008030% | 0.00012100 |

## Findings

1. The physical reference weight remains near unity, but the CS weights have a long rare-event tail.
2. The source-aware tail is materially larger than the NPS-only tail. Splitting on `w*_NPS` alone would miss much of the perturbed-source contribution.
3. The absolute maximum is sample-size dependent and unstable. A percentile-based bound is more useful.
4. The measured NPS/fission distributions are identical: True. The source/score distributions are identical: True. This is consistent with `nuSigmaF*/nuSigmaF = 1` at contributing fission collisions for this CR-GT perturbation.
5. A first A/B implementation should test `W_upper = 3` with a split cap of 4 or 8, then compare variance of delta-k (or slope) per wall-clock time. `W_upper = 6` is a more conservative alternative.

## Unbiased implementation requirements

- Trigger on the selected source-aware effective particle weight `w exp(log likelihood) R_L`.
- For `n` children, divide only the common physical particle weight `w` by `n`; copy `csLogLikelihood` and every `R_L` unchanged. This preserves both reference and perturbed total weight.
- Assign independent RNG streams to the split children.
- Apply a hard split-count cap and check bank capacity before insertion.
- Do not clip `csLogLikelihood`, `R_L`, or the perturbed score; clipping would bias the estimator.
- Treat this run as bound selection. Claim variance reduction only after an otherwise identical split/no-split A/B run.
