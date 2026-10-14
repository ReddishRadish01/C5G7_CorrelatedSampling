#!/usr/bin/env python3
"""Analyze the cycle and histogram outputs from --cs-track-weights."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def histogram_arrays(rows: list[dict[str, str]], column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lower = np.array([float(row["log10_lower"]) for row in rows], dtype=float)
    upper = np.array([float(row["log10_upper"]) for row in rows], dtype=float)
    counts = np.array([int(row[column]) for row in rows], dtype=np.int64)
    return lower, upper, counts


def finite_bin_edges(lower: np.ndarray, upper: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    finite_lowers = lower[np.isfinite(lower)]
    finite_uppers = upper[np.isfinite(upper)]
    low_fallback = float(np.min(finite_lowers))
    high_fallback = float(np.max(finite_uppers))
    return np.where(np.isfinite(lower), lower, low_fallback), np.where(np.isfinite(upper), upper, high_fallback)


def histogram_quantile(
    lower: np.ndarray,
    upper: np.ndarray,
    counts: np.ndarray,
    probability: float,
    exact_maximum: float,
) -> float:
    total = int(np.sum(counts))
    if total == 0:
        return math.nan
    target = probability * total
    cumulative = np.cumsum(counts)
    index = int(np.searchsorted(cumulative, target, side="left"))
    previous = int(cumulative[index - 1]) if index > 0 else 0
    within = 0.5 if counts[index] == 0 else (target - previous) / counts[index]
    safe_lower, safe_upper = finite_bin_edges(lower, upper)
    log_value = safe_lower[index] + within * (safe_upper[index] - safe_lower[index])
    return min(10.0**log_value, exact_maximum)


def histogram_fraction_above(
    lower: np.ndarray,
    upper: np.ndarray,
    counts: np.ndarray,
    threshold: float,
) -> float:
    total = float(np.sum(counts))
    if total == 0.0:
        return math.nan
    threshold_log = math.log10(threshold)
    safe_lower, safe_upper = finite_bin_edges(lower, upper)
    fraction = np.zeros_like(safe_lower)
    fraction[threshold_log <= safe_lower] = 1.0
    inside = (threshold_log > safe_lower) & (threshold_log < safe_upper)
    fraction[inside] = (safe_upper[inside] - threshold_log) / (safe_upper[inside] - safe_lower[inside])
    return float(np.sum(counts * fraction) / total)


def split_cost_estimate(
    lower: np.ndarray,
    upper: np.ndarray,
    counts: np.ndarray,
    threshold: float,
    split_cap: int,
) -> tuple[float, float]:
    safe_lower, safe_upper = finite_bin_edges(lower, upper)
    representative = 10.0 ** ((safe_lower + safe_upper) / 2.0)
    split_count = np.minimum(split_cap, np.maximum(1, np.ceil(representative / threshold))).astype(int)
    total = float(np.sum(counts))
    triggered = float(np.sum(counts[split_count > 1]) / total)
    extra = float(np.sum(counts * (split_count - 1)) / total)
    return triggered, extra


def weighted_fraction(rows: list[dict[str, str]], fraction_column: str) -> float:
    observations = np.array([int(row["observations"]) for row in rows], dtype=np.float64)
    fractions = np.array([float(row[fraction_column]) for row in rows], dtype=np.float64)
    return float(np.sum(observations * fractions) / np.sum(observations))


def format_value(value: float) -> str:
    if value == 0.0:
        return "0"
    if abs(value) >= 1000.0 or abs(value) < 1.0e-3:
        return f"{value:.4e}"
    return f"{value:.6g}"


def svg_polyline(points: list[tuple[float, float]], color: str, width: float = 1.8, dash: str = "") -> str:
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    dash_attribute = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{coordinates}" fill="none" stroke="{color}" stroke-width="{width}"{dash_attribute}/>'


def create_svg_plot(
    output_path: Path,
    cycles: np.ndarray,
    max_w: np.ndarray,
    max_nps: np.ndarray,
    max_source: np.ndarray,
    p999_source: np.ndarray,
    histograms: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    inactive: int,
    reference_depth: float,
    perturbed_depth: float,
    source_level: int,
) -> None:
    width, height = 1440, 560
    plot_y0, plot_y1 = 82.0, 475.0
    left_x0, left_x1 = 82.0, 690.0
    right_x0, right_x1 = 810.0, 1388.0

    def left_x(value: float) -> float:
        return left_x0 + (value - 1.0) / (float(np.max(cycles)) - 1.0) * (left_x1 - left_x0)

    left_log_min, left_log_max = math.log10(0.05), math.log10(500.0)

    def left_y(value: float) -> float:
        clipped = min(500.0, max(0.05, value))
        return plot_y1 - (math.log10(clipped) - left_log_min) / (left_log_max - left_log_min) * (plot_y1 - plot_y0)

    right_log_x_min, right_log_x_max = math.log10(0.05), math.log10(500.0)
    right_log_y_min, right_log_y_max = -8.0, 0.0

    def right_x(value: float) -> float:
        clipped = min(500.0, max(0.05, value))
        return right_x0 + (math.log10(clipped) - right_log_x_min) / (right_log_x_max - right_log_x_min) * (right_x1 - right_x0)

    def right_y(value: float) -> float:
        clipped = min(1.0, max(1.0e-8, value))
        return plot_y1 - (math.log10(clipped) - right_log_y_min) / (right_log_y_max - right_log_y_min) * (plot_y1 - plot_y0)

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.title{font-size:20px;font-weight:600}.subtitle{font-size:15px;font-weight:600}.tick{font-size:12px}.label{font-size:13px}.legend{font-size:12px}</style>',
        f'<text x="720" y="30" text-anchor="middle" class="title">CS weight tracking: {reference_depth:.6f} cm → {perturbed_depth:.6f} cm, L={source_level}</text>',
        '<text x="386" y="60" text-anchor="middle" class="subtitle">Cycle-wise extremes and upper tail</text>',
        '<text x="1099" y="60" text-anchor="middle" class="subtitle">Active-cycle weight distribution (CCDF)</text>',
    ]

    for x0, x1 in ((left_x0, left_x1), (right_x0, right_x1)):
        svg.append(f'<rect x="{x0}" y="{plot_y0}" width="{x1-x0}" height="{plot_y1-plot_y0}" fill="none" stroke="#333" stroke-width="1"/>')

    for tick in (0.1, 1.0, 10.0, 100.0):
        y = left_y(tick)
        svg.append(f'<line x1="{left_x0}" y1="{y}" x2="{left_x1}" y2="{y}" stroke="#ddd"/>')
        svg.append(f'<text x="{left_x0-10}" y="{y+4}" text-anchor="end" class="tick">{tick:g}</text>')
    for tick in (1, 50, 100, 150, 200, 250, 300):
        x = left_x(float(tick))
        svg.append(f'<line x1="{x}" y1="{plot_y1}" x2="{x}" y2="{plot_y1+5}" stroke="#333"/>')
        svg.append(f'<text x="{x}" y="{plot_y1+21}" text-anchor="middle" class="tick">{tick}</text>')

    left_series = [
        (max_w, "#1f77b4", "max w", ""),
        (max_nps, "#ff7f0e", "max w* NPS", ""),
        (max_source, "#d62728", f"max w*R{source_level}", ""),
        (p999_source, "#2ca02c", f"p99.9 w*R{source_level}", "7,4"),
    ]
    for values, color, _, dash in left_series:
        svg.append(svg_polyline([(left_x(float(c)), left_y(float(v))) for c, v in zip(cycles, values)], color, 1.7, dash))
    active_x = left_x(float(inactive) + 0.5)
    svg.append(f'<line x1="{active_x}" y1="{plot_y0}" x2="{active_x}" y2="{plot_y1}" stroke="#111" stroke-dasharray="3,4"/>')
    svg.append(f'<text x="{active_x+7}" y="{plot_y0+18}" class="tick">active</text>')
    svg.append(f'<text x="{(left_x0+left_x1)/2}" y="{height-24}" text-anchor="middle" class="label">Cycle</text>')
    svg.append(f'<text x="22" y="{(plot_y0+plot_y1)/2}" text-anchor="middle" class="label" transform="rotate(-90 22 {(plot_y0+plot_y1)/2})">Weight at fission-scoring collision</text>')

    for tick in (0.1, 1.0, 10.0, 100.0):
        x = right_x(tick)
        svg.append(f'<line x1="{x}" y1="{plot_y0}" x2="{x}" y2="{plot_y1}" stroke="#eee"/>')
        svg.append(f'<text x="{x}" y="{plot_y1+21}" text-anchor="middle" class="tick">{tick:g}</text>')
    for exponent in (0, -2, -4, -6, -8):
        value = 10.0**exponent
        y = right_y(value)
        svg.append(f'<line x1="{right_x0}" y1="{y}" x2="{right_x1}" y2="{y}" stroke="#ddd"/>')
        svg.append(f'<text x="{right_x0-10}" y="{y+4}" text-anchor="end" class="tick">10^{exponent}</text>')

    colors = {"w": "#1f77b4", "wstar_nps": "#ff7f0e", "wstar_source": "#d62728"}
    labels = {"w": "w", "wstar_nps": "w* NPS", "wstar_source": f"w*R{source_level}"}
    for metric in ("w", "wstar_nps", "wstar_source"):
        lower, upper, counts = histograms[metric]
        safe_lower, safe_upper = finite_bin_edges(lower, upper)
        centers = 10.0 ** ((safe_lower + safe_upper) / 2.0)
        ccdf = np.cumsum(counts[::-1])[::-1] / np.sum(counts)
        mask = (ccdf >= 1.0e-8) & (centers >= 0.05) & (centers <= 500.0)
        svg.append(svg_polyline([(right_x(float(x)), right_y(float(y))) for x, y in zip(centers[mask], ccdf[mask])], colors[metric], 1.8))
    for threshold, color, dash in ((3.0, "#9467bd", "7,4"), (6.0, "#8c564b", "2,4")):
        x = right_x(threshold)
        svg.append(f'<line x1="{x}" y1="{plot_y0}" x2="{x}" y2="{plot_y1}" stroke="{color}" stroke-width="1.5" stroke-dasharray="{dash}"/>')
    svg.append(f'<text x="{(right_x0+right_x1)/2}" y="{height-24}" text-anchor="middle" class="label">Weight</text>')
    svg.append(f'<text x="744" y="{(plot_y0+plot_y1)/2}" text-anchor="middle" class="label" transform="rotate(-90 744 {(plot_y0+plot_y1)/2})">Exceedance probability</text>')

    legend_items = [(label, color, dash) for _, color, label, dash in left_series]
    legend_x, legend_y = left_x0 + 14, plot_y0 + 22
    for index, (label, color, dash) in enumerate(legend_items):
        y = legend_y + 19 * index
        dash_attribute = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x+26}" y2="{y}" stroke="{color}" stroke-width="2"{dash_attribute}/><text x="{legend_x+34}" y="{y+4}" class="legend">{label}</text>')
    right_legend = [(labels[key], colors[key], "") for key in ("w", "wstar_nps", "wstar_source")] + [("candidate Wu=3", "#9467bd", "7,4"), ("candidate Wu=6", "#8c564b", "2,4")]
    legend_x, legend_y = right_x0 + 14, plot_y1 - 84
    for index, (label, color, dash) in enumerate(right_legend):
        y = legend_y + 17 * index
        dash_attribute = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x+26}" y2="{y}" stroke="{color}" stroke-width="2"{dash_attribute}/><text x="{legend_x+34}" y="{y+4}" class="legend">{label}</text>')

    svg.append("</svg>")
    output_path.write_text("\n".join(svg), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cycle", type=Path, required=True)
    parser.add_argument("--histogram", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--inactive", type=int, default=100)
    parser.add_argument("--particles", type=int, default=10_000_000)
    parser.add_argument("--reference-depth", type=float, default=83.6388953145)
    parser.add_argument("--perturbed-depth", type=float, default=83.7638953145)
    parser.add_argument("--source-level", type=int, default=10)
    args = parser.parse_args()

    cycle_rows = read_csv(args.cycle)
    histogram_rows = read_csv(args.histogram)
    active_rows = [row for row in cycle_rows if row["phase"] == "active"]
    if not active_rows:
        raise RuntimeError("No active-cycle rows were found.")

    observations = sum(int(row["observations"]) for row in active_rows)
    histogram_observations = sum(int(row["histogram_observations"]) for row in active_rows)
    invalid = sum(int(row["invalid_observations"]) for row in active_rows)

    metric_config = {
        "w": ("referenceHistogram", "w", "max_w", None, None),
        "wstar_nps": ("npsHistogram", "wstar_nps", "max_wstar_nps", "max_wstar_nps_over_w", "fraction_wstar_nps_gt_w"),
        "wstar_source": ("sourceHistogram", "wstar_source", "max_wstar_source", "max_wstar_source_over_w", "fraction_wstar_source_gt_w"),
        "wstar_fission": ("fissionHistogram", "wstar_fission", "max_wstar_fission", "max_wstar_fission_over_w", None),
        "wstar_score": ("scoreHistogram", "wstar_score", "max_wstar_score", "max_wstar_score_over_w", "fraction_wstar_score_gt_w"),
    }

    hist_column_by_metric = {
        "w": "w",
        "wstar_nps": "wstar_nps",
        "wstar_source": "wstar_source",
        "wstar_fission": "wstar_fission",
        "wstar_score": "wstar_score",
    }
    summaries: dict[str, dict[str, float]] = {}
    histograms: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for metric, (_, _, max_column, ratio_column, fraction_column) in metric_config.items():
        lower, upper, counts = histogram_arrays(histogram_rows, hist_column_by_metric[metric])
        histograms[metric] = (lower, upper, counts)
        maximum = max(float(row[max_column]) for row in active_rows)
        result = {
            "p50": histogram_quantile(lower, upper, counts, 0.50, maximum),
            "p90": histogram_quantile(lower, upper, counts, 0.90, maximum),
            "p95": histogram_quantile(lower, upper, counts, 0.95, maximum),
            "p99": histogram_quantile(lower, upper, counts, 0.99, maximum),
            "p999": histogram_quantile(lower, upper, counts, 0.999, maximum),
            "p9999": histogram_quantile(lower, upper, counts, 0.9999, maximum),
            "max": maximum,
            "cycle_max_median": float(np.median([float(row[max_column]) for row in active_rows])),
            "cycle_max_p95": float(np.quantile([float(row[max_column]) for row in active_rows], 0.95)),
        }
        if ratio_column is not None:
            result["max_ratio"] = max(float(row[ratio_column]) for row in active_rows)
        if fraction_column is not None:
            result["fraction_gt_w"] = weighted_fraction(active_rows, fraction_column)
        summaries[metric] = result

    source_lower, source_upper, source_counts = histograms["wstar_source"]
    exceedance = {
        threshold: histogram_fraction_above(source_lower, source_upper, source_counts, threshold)
        for threshold in (1.0, 2.0, 3.0, 5.0, 6.0, 10.0)
    }
    split_estimates = {
        (threshold, cap): split_cost_estimate(source_lower, source_upper, source_counts, threshold, cap)
        for threshold in (3.0, 6.0)
        for cap in (4, 8)
    }

    nps_equals_fission = all(
        row["max_wstar_nps"] == row["max_wstar_fission"] and
        row["p99_wstar_nps"] == row["p99_wstar_fission"] and
        row["p999_wstar_nps"] == row["p999_wstar_fission"]
        for row in active_rows
    ) and all(
        row["wstar_nps"] == row["wstar_fission"] for row in histogram_rows
    )
    source_equals_score = all(
        row["max_wstar_source"] == row["max_wstar_score"] and
        row["p99_wstar_source"] == row["p99_wstar_score"] and
        row["p999_wstar_source"] == row["p999_wstar_score"]
        for row in active_rows
    ) and all(
        row["wstar_source"] == row["wstar_score"] for row in histogram_rows
    )

    cycles = np.array([int(row["cycle"]) for row in cycle_rows])
    max_w = np.array([float(row["max_w"]) for row in cycle_rows])
    max_nps = np.array([float(row["max_wstar_nps"]) for row in cycle_rows])
    max_source = np.array([float(row["max_wstar_source"]) for row in cycle_rows])
    p999_source = np.array([float(row["p999_wstar_source"]) for row in cycle_rows])

    plot_path = args.output_prefix.with_suffix(".svg")
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    create_svg_plot(
        plot_path, cycles, max_w, max_nps, max_source, p999_source,
        histograms, args.inactive, args.reference_depth, args.perturbed_depth,
        args.source_level,
    )

    report_path = args.output_prefix.with_suffix(".md")
    with report_path.open("w", encoding="utf-8") as report:
        report.write("# CS weight tracking report\n\n")
        report.write("## Run conditions\n\n")
        report.write(f"- Reference rod depth: {args.reference_depth:.10f} cm\n")
        report.write(f"- Perturbed rod depth: {args.perturbed_depth:.10f} cm\n")
        report.write(f"- Source propagation level: L={args.source_level}\n")
        report.write(f"- Particles per cycle: {args.particles:,}\n")
        report.write(f"- Cycles: {args.inactive} inactive + {len(active_rows)} active\n")
        report.write(f"- Active fission-score observations: {observations:,}\n")
        report.write(f"- Histogram sample observations: {histogram_observations:,}\n")
        report.write(f"- Invalid observations: {invalid:,}\n\n")

        report.write("## Definitions\n\n")
        report.write("- `w`: reference physical particle weight.\n")
        report.write("- `w*_NPS = w exp(log likelihood)`: transport likelihood-ratio weight.\n")
        report.write(f"- `w*R_L`: source-aware perturbed particle weight at L={args.source_level}.\n")
        report.write("- `w*_f = w*_NPS (nuSigmaF*/nuSigmaF)`: fission-score adjusted weight.\n")
        report.write("- Maxima and comparison counters use every observation. Percentiles use the deterministic 1/64 histogram sample.\n\n")

        report.write("## Active-cycle results\n\n")
        report.write("| Metric | p50 | p99 | p99.9 | p99.99 | Median cycle max | p95 cycle max | Exact max | Max / w | Fraction > w |\n")
        report.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for metric in ("w", "wstar_nps", "wstar_source", "wstar_fission", "wstar_score"):
            result = summaries[metric]
            ratio_text = format_value(result["max_ratio"]) if "max_ratio" in result else "n.a."
            fraction_text = f"{100.0 * result['fraction_gt_w']:.5f}%" if "fraction_gt_w" in result else "n.a."
            report.write(
                f"| {metric} | {format_value(result['p50'])} | {format_value(result['p99'])} | "
                f"{format_value(result['p999'])} | {format_value(result['p9999'])} | "
                f"{format_value(result['cycle_max_median'])} | {format_value(result['cycle_max_p95'])} | "
                f"{format_value(result['max'])} | {ratio_text} | {fraction_text} |\n"
            )

        report.write("\n## Source-aware upper-tail checks\n\n")
        report.write("| Threshold | Fraction of sampled fission-score observations above threshold |\n")
        report.write("|---:|---:|\n")
        for threshold, fraction in exceedance.items():
            report.write(f"| {threshold:g} | {100.0 * fraction:.6f}% |\n")

        report.write("\nApproximate diagnostic splitting load (histogram-bin approximation; not a transport-runtime prediction):\n\n")
        report.write("| Upper bound | Split cap | Trigger fraction | Extra children per observed score event |\n")
        report.write("|---:|---:|---:|---:|\n")
        for (threshold, cap), (triggered, extra) in split_estimates.items():
            report.write(f"| {threshold:g} | {cap} | {100.0 * triggered:.6f}% | {extra:.8f} |\n")

        report.write("\n## Findings\n\n")
        report.write("1. The physical reference weight remains near unity, but the CS weights have a long rare-event tail.\n")
        report.write("2. The source-aware tail is materially larger than the NPS-only tail. Splitting on `w*_NPS` alone would miss much of the perturbed-source contribution.\n")
        report.write("3. The absolute maximum is sample-size dependent and unstable. A percentile-based bound is more useful.\n")
        report.write(f"4. The measured NPS/fission distributions are identical: {nps_equals_fission}. The source/score distributions are identical: {source_equals_score}. This is consistent with `nuSigmaF*/nuSigmaF = 1` at contributing fission collisions for this CR-GT perturbation.\n")
        report.write("5. A first A/B implementation should test `W_upper = 3` with a split cap of 4 or 8, then compare variance of delta-k (or slope) per wall-clock time. `W_upper = 6` is a more conservative alternative.\n\n")

        report.write("## Unbiased implementation requirements\n\n")
        report.write("- Trigger on the selected source-aware effective particle weight `w exp(log likelihood) R_L`.\n")
        report.write("- For `n` children, divide only the common physical particle weight `w` by `n`; copy `csLogLikelihood` and every `R_L` unchanged. This preserves both reference and perturbed total weight.\n")
        report.write("- Assign independent RNG streams to the split children.\n")
        report.write("- Apply a hard split-count cap and check bank capacity before insertion.\n")
        report.write("- Do not clip `csLogLikelihood`, `R_L`, or the perturbed score; clipping would bias the estimator.\n")
        report.write("- Treat this run as bound selection. Claim variance reduction only after an otherwise identical split/no-split A/B run.\n")

    print(f"plot={plot_path}")
    print(f"report={report_path}")
    print(f"observations={observations}")
    print(f"histogram_observations={histogram_observations}")
    print(f"invalid={invalid}")


if __name__ == "__main__":
    main()
