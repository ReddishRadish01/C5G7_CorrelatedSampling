# CS prediction validation against independent direct Monte Carlo

- Cases: 6 reference depths; 4/6 passed all enabled checks
- Confidence threshold: 2 sigma
- Source propagation: L10
- Each row uses three independent odd-seed calculations: direct MC at d, direct MC at d+h, and CS from d to d+h.

Two complementary comparisons are reported. `delta z` compares the paired CS delta-k with the difference of the two direct MC k values. `target z` compares the standalone CS prediction k*(d+h) with the independent direct MC value at d+h.

| d (cm) | d+h (cm) | direct slope (/cm) | CS slope (/cm) | delta z | target z | source tail | status |
|---:|---:|---:|---:|---:|---:|---|---|
| 64.000000 | 64.125000 | -2.3712000e-03 ± 2.01e-04 | -3.0133563e-03 ± 2.47e-05 | -3.18 | -1.36 | pass | DELTA_K_DISAGREE |
| 72.000000 | 72.125000 | -3.7336000e-03 ± 2.51e-04 | -4.0968634e-03 ± 2.43e-05 | -1.44 | -1.33 | pass | PASS |
| 80.000000 | 80.125000 | -5.4808000e-03 ± 2.59e-04 | -5.5917108e-03 ± 3.30e-05 | -0.43 | -0.85 | pass | PASS |
| 83.638895 | 83.763895 | -6.1976000e-03 ± 2.42e-04 | -6.2914609e-03 ± 3.30e-05 | -0.38 | 0.24 | pass | PASS |
| 88.000000 | 88.125000 | -7.2336000e-03 ± 2.23e-04 | -6.8244248e-03 ± 3.74e-05 | 1.81 | 1.25 | pass | PASS |
| 92.500000 | 92.625000 | -6.9000000e-03 ± 1.97e-04 | -6.1670357e-03 ± 3.05e-05 | 3.67 | 3.72 | pass | DELTA_K_DISAGREE |

Interpretation: agreement within ±2 sigma is statistically compatible, not proof of exact equality. A source-tail failure is kept separate because an apparently good direct comparison can occur by chance while L10 has not stabilized.

The direct finite-difference slope contains the independent uncertainties of both endpoint calculations. The CS slope uses the paired delta-k uncertainty, which preserves the correlation benefit of CS.
