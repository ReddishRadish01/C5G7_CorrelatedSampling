# CS prediction validation against independent direct Monte Carlo

- Cases: 13 reference depths; 9/13 passed all enabled checks
- Confidence threshold: 2 sigma
- Source propagation: L10
- Each row uses three independent odd-seed calculations: direct MC at d, direct MC at d+h, and CS from d to d+h.

Two complementary comparisons are reported. `delta z` compares the paired CS delta-k with the difference of the two direct MC k values. `target z` compares the standalone CS prediction k*(d+h) with the independent direct MC value at d+h.

| d (cm) | d+h (cm) | direct slope (/cm) | CS slope (/cm) | delta z | target z | source tail | status |
|---:|---:|---:|---:|---:|---:|---|---|
| 0.000000 | 0.125000 | 5.1360000e-04 ± 2.49e-04 | 7.7733967e-08 ± 6.71e-08 | -2.07 | -1.37 | pass | MONOTONICITY_FAIL |
| 8.000000 | 8.125000 | 6.2400000e-05 ± 1.97e-04 | -9.0690701e-08 ± 2.38e-07 | -0.32 | -0.46 | pass | PASS |
| 16.000000 | 16.125000 | -1.7600000e-05 ± 2.18e-04 | -3.4049765e-05 ± 2.09e-06 | -0.08 | 0.54 | pass | PASS |
| 24.000000 | 24.125000 | -4.0320000e-04 ± 2.09e-04 | -3.6801598e-04 ± 1.14e-05 | 0.17 | 1.44 | pass | PASS |
| 32.000000 | 32.125000 | -4.9840000e-04 ± 1.70e-04 | -7.9589479e-04 ± 1.48e-05 | -1.74 | 0.85 | pass | PASS |
| 40.000000 | 40.125000 | -1.4400000e-03 ± 2.08e-04 | -1.1787910e-03 ± 1.49e-05 | 1.25 | 0.91 | pass | PASS |
| 48.000000 | 48.125000 | -1.6216000e-03 ± 2.51e-04 | -1.6460189e-03 ± 1.72e-05 | -0.10 | 0.91 | pass | PASS |
| 56.000000 | 56.125000 | -2.0784000e-03 ± 2.16e-04 | -2.2068132e-03 ± 1.89e-05 | -0.59 | -0.45 | pass | PASS |
| 96.000000 | 96.125000 | -4.7720000e-03 ± 1.67e-04 | -4.3138906e-03 ± 3.33e-05 | 2.69 | 1.12 | pass | DELTA_K_DISAGREE |
| 104.000000 | 104.125000 | -7.7520000e-04 ± 1.79e-04 | -7.7800383e-04 ± 1.23e-05 | -0.02 | 2.24 | pass | TARGET_K_DISAGREE |
| 112.000000 | 112.125000 | 2.8800000e-05 ± 2.04e-04 | -1.7408297e-15 ± 5.26e-16 | -0.14 | 0.57 | pass | PASS |
| 120.000000 | 120.125000 | 4.8800000e-05 ± 1.79e-04 | 1.2615170e-09 ± 1.29e-09 | -0.27 | -2.22 | pass | TARGET_K_DISAGREE |
| 128.395000 | 128.520000 | -1.0080000e-04 ± 1.76e-04 | -6.0890848e-13 ± 6.13e-13 | 0.57 | 0.30 | pass | PASS |

Interpretation: agreement within ±2 sigma is statistically compatible, not proof of exact equality. A source-tail failure is kept separate because an apparently good direct comparison can occur by chance while L10 has not stabilized.

The direct finite-difference slope contains the independent uncertainties of both endpoint calculations. The CS slope uses the paired delta-k uncertainty, which preserves the correlation benefit of CS.
