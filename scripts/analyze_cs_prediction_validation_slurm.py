#!/usr/bin/env python3
"""Postprocess the multi-GPU CS prediction validation without third-party packages."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from pathlib import Path
from typing import Iterable


NUMBER = r"[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?"


def last_match(text: str, pattern: str, group: str = "value") -> float:
    matches = list(re.finditer(pattern, text, re.MULTILINE))
    if not matches:
        raise ValueError(f"Could not parse pattern: {pattern}")
    return float(matches[-1].group(group))


def conservative_stderr(values: list[float]) -> tuple[float, float, float, int, int]:
    if len(values) < 2:
        return 0.0, 0.0, 0.0, 0, 0
    iid = statistics.stdev(values) / math.sqrt(len(values))
    batch_size = max(2, int(math.floor(math.sqrt(len(values)))))
    batch_count = len(values) // batch_size
    if batch_count < 2:
        return iid, iid, iid, 0, 0
    batch_means = [
        statistics.fmean(values[i * batch_size : (i + 1) * batch_size])
        for i in range(batch_count)
    ]
    batch_error = statistics.stdev(batch_means) / math.sqrt(batch_count)
    return max(iid, batch_error), iid, batch_error, batch_size, batch_count


def read_direct(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return {
        "requested_depth": last_match(text, rf"^Requested rod depth:\s*(?P<value>{NUMBER})\s+cm\s*$"),
        "actual_depth": last_match(text, rf"^Actual modeled rod depth:\s*(?P<value>{NUMBER})\s+cm\s*$"),
        "k": last_match(text, rf"^k_mean:\s*(?P<value>{NUMBER})\s*$"),
        "k_stderr": last_match(text, rf"^k_stderr\(mean\):\s*(?P<value>{NUMBER})\s*$"),
    }


def find_cs_history(work_dir: Path) -> Path:
    candidates = sorted((work_dir / "results").glob("cs_history_*.txt"))
    if len(candidates) != 1:
        raise ValueError(
            f"Expected one CS history in {work_dir / 'results'}, found {len(candidates)}"
        )
    return candidates[0]


def read_cs(path: Path, work_dir: Path, inactive: int, level: int) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    history_path = find_cs_history(work_dir)
    lines = history_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise ValueError(f"Empty CS history: {history_path}")
    header = lines[0].lstrip("# ").split()
    column = f"k_cs_L{level}"
    if column not in header:
        raise ValueError(f"Missing {column} in {history_path}")
    cycle_index = header.index("cycle")
    k_index = header.index(column)
    first_usable_cycle = max(inactive, level) + 1
    samples: list[float] = []
    for line in lines[1:]:
        fields = line.split()
        if len(fields) <= max(cycle_index, k_index):
            continue
        if int(fields[cycle_index]) >= first_usable_cycle:
            samples.append(float(fields[k_index]))
    predicted_k_stderr, pred_iid, pred_batch, batch_size, batch_count = conservative_stderr(samples)
    reported_count = int(last_match(text, r"^Usable active cycles:\s*(?P<value>[0-9]+)\s*$"))
    if len(samples) != reported_count:
        raise ValueError(
            f"CS history has {len(samples)} usable samples but summary reports {reported_count}: {path}"
        )
    return {
        "reference_depth": last_match(
            text, rf"^Reference rod depth \(cm\):\s*(?P<value>{NUMBER})\s*$"
        ),
        "target_depth": last_match(
            text, rf"^Perturbed rod depth \(cm\):\s*(?P<value>{NUMBER})\s*$"
        ),
        "paired_reference_k": last_match(text, rf"^Paired reference k:\s*(?P<value>{NUMBER})\s*$"),
        "predicted_target_k": last_match(
            text,
            rf"^Perturbed k, NPS \+ fission-source propagation.*:\s*(?P<value>{NUMBER})\s*$",
        ),
        "predicted_target_k_stderr": predicted_k_stderr,
        "predicted_target_k_stderr_iid": pred_iid,
        "predicted_target_k_stderr_batch": pred_batch,
        "predicted_target_k_batch_size": batch_size,
        "predicted_target_k_batch_count": batch_count,
        "delta_k": last_match(
            text,
            rf"^Total delta k:\s*(?P<value>{NUMBER})\s+\+/-\s+(?P<stderr>{NUMBER})\s*$",
        ),
        "delta_k_stderr": last_match(
            text,
            rf"^Total delta k:\s*(?P<value>{NUMBER})\s+\+/-\s+(?P<stderr>{NUMBER})\s*$",
            "stderr",
        ),
        "delta_rho_pcm": last_match(
            text,
            rf"^Total delta reactivity \(pcm\):\s*(?P<value>{NUMBER})\s+\+/-\s+(?P<stderr>{NUMBER})\s*$",
        ),
        "delta_rho_stderr_pcm": last_match(
            text,
            rf"^Total delta reactivity \(pcm\):\s*(?P<value>{NUMBER})\s+\+/-\s+(?P<stderr>{NUMBER})\s*$",
            "stderr",
        ),
        "last_increment": last_match(
            text,
            rf"^Last propagation increment,.*:\s*(?P<value>{NUMBER})\s+\+/-\s+(?P<stderr>{NUMBER})\s*$",
        ),
        "last_increment_stderr": last_match(
            text,
            rf"^Last propagation increment,.*:\s*(?P<value>{NUMBER})\s+\+/-\s+(?P<stderr>{NUMBER})\s*$",
            "stderr",
        ),
        "diagnostics": int(
            last_match(text, r"^Support/non-finite diagnostic count:\s*(?P<value>[0-9]+)\s*$")
        ),
    }


def quadrature(*values: float) -> float:
    return math.sqrt(sum(value * value for value in values))


def z_score(difference: float, stderr: float) -> float:
    if stderr == 0.0:
        return math.copysign(math.inf, difference) if difference else 0.0
    return difference / stderr


def sign_assessment(mean: float, stderr: float, expected_sign: int, sigma: float) -> str:
    low, high = mean - sigma * stderr, mean + sigma * stderr
    if expected_sign < 0:
        if high < 0.0:
            return "RESOLVED_EXPECTED"
        if low > 0.0:
            return "RESOLVED_WRONG"
    else:
        if low > 0.0:
            return "RESOLVED_EXPECTED"
        if high < 0.0:
            return "RESOLVED_WRONG"
    return "UNRESOLVED"


def agreement_label(z_value: float) -> str:
    if abs(z_value) <= 1.0:
        return "AGREE_1SIGMA"
    if abs(z_value) <= 2.0:
        return "AGREE_2SIGMA"
    return "DISAGREE_GT_2SIGMA"


def svg_plot(rows: list[dict[str, object]], path: Path, confidence_sigma: float) -> None:
    width, height = 1180, 760
    left, right = 105.0, 45.0
    top1, bottom1 = 80.0, 405.0
    top2, bottom2 = 520.0, 680.0
    x_values = [float(row["reference_depth_cm"]) for row in rows]
    x_min, x_max = min(x_values), max(x_values)
    if x_max == x_min:
        x_min -= 0.5
        x_max += 0.5
    slope_values: list[float] = []
    for row in rows:
        for prefix in ("direct", "cs"):
            mean = float(row[f"{prefix}_slope_per_cm"])
            error = float(row[f"{prefix}_slope_stderr_per_cm"])
            slope_values.extend((mean - error, mean + error))
    y_min, y_max = min(slope_values), max(slope_values)
    pad = max(1.0e-8, 0.1 * (y_max - y_min))
    y_min, y_max = y_min - pad, y_max + pad
    plot_width = width - left - right

    def xp(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def yp(value: float) -> float:
        return bottom1 - (value - y_min) / (y_max - y_min) * (bottom1 - top1)

    z_extent = max(3.0, confidence_sigma + 0.5)

    def zp(value: float) -> float:
        clipped = max(-z_extent, min(z_extent, value))
        return bottom2 - (clipped + z_extent) / (2.0 * z_extent) * (bottom2 - top2)

    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#1f2937}.title{font-size:24px;font-weight:700}.label{font-size:15px}.tick{font-size:12px;fill:#4b5563}.grid{stroke:#d1d5db;stroke-width:1}.axis{stroke:#111827;stroke-width:1.5}</style>',
        '<text x="105" y="38" class="title">CS slope and next-state prediction validation</text>',
    ]
    for index in range(6):
        x_value = x_min + (x_max - x_min) * index / 5.0
        x = xp(x_value)
        out.append(f'<line x1="{x:.2f}" y1="{top1}" x2="{x:.2f}" y2="{bottom2}" class="grid"/>')
        out.append(f'<text x="{x:.2f}" y="{bottom2 + 25}" text-anchor="middle" class="tick">{x_value:.3f}</text>')
    for index in range(6):
        y_value = y_min + (y_max - y_min) * index / 5.0
        y = yp(y_value)
        out.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" class="grid"/>')
        out.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" class="tick">{y_value:.3e}</text>')
    for z_value in (-confidence_sigma, 0.0, confidence_sigma):
        y = zp(z_value)
        dash = "6,5" if z_value else ""
        out.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" stroke="#9ca3af" stroke-dasharray="{dash}"/>')
        out.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" class="tick">{z_value:g}</text>')
    out.extend(
        [
            f'<line x1="{left}" y1="{top1}" x2="{left}" y2="{bottom1}" class="axis"/>',
            f'<line x1="{left}" y1="{bottom1}" x2="{width-right}" y2="{bottom1}" class="axis"/>',
            f'<line x1="{left}" y1="{top2}" x2="{left}" y2="{bottom2}" class="axis"/>',
            f'<line x1="{left}" y1="{bottom2}" x2="{width-right}" y2="{bottom2}" class="axis"/>',
            f'<text x="{width/2:.1f}" y="735" text-anchor="middle" class="label">Reference rod insertion depth (cm)</text>',
            '<text x="22" y="245" transform="rotate(-90 22 245)" text-anchor="middle" class="label">dk/dd (cm^-1)</text>',
            '<text x="22" y="600" transform="rotate(-90 22 600)" text-anchor="middle" class="label">CS target-k agreement z</text>',
            '<text x="755" y="62" class="label" fill="#1769aa">● Direct MC finite difference (±1σ)</text>',
            '<text x="755" y="84" class="label" fill="#d97706">■ CS paired estimate (±1σ)</text>',
        ]
    )
    for row in rows:
        x = xp(float(row["reference_depth_cm"]))
        for prefix, offset, color, shape in (
            ("direct", -5.0, "#1769aa", "circle"),
            ("cs", 5.0, "#d97706", "square"),
        ):
            mean = float(row[f"{prefix}_slope_per_cm"])
            error = float(row[f"{prefix}_slope_stderr_per_cm"])
            low, high, center = yp(mean - error), yp(mean + error), yp(mean)
            out.append(f'<line x1="{x+offset:.2f}" y1="{low:.2f}" x2="{x+offset:.2f}" y2="{high:.2f}" stroke="{color}"/>')
            if shape == "circle":
                out.append(f'<circle cx="{x+offset:.2f}" cy="{center:.2f}" r="5" fill="{color}"/>')
            else:
                out.append(f'<rect x="{x+offset-4.5:.2f}" y="{center-4.5:.2f}" width="9" height="9" fill="{color}"/>')
        z_value = float(row["predicted_target_k_agreement_z"])
        color = "#15803d" if abs(z_value) <= confidence_sigma else "#b91c1c"
        out.append(f'<circle cx="{x:.2f}" cy="{zp(z_value):.2f}" r="5.5" fill="{color}"/>')
    out.append(
        f'<text x="105" y="488" class="label">Lower panel: [CS predicted k(d+h) - direct MC k(d+h)] / combined standard error; dashed lines are ±{confidence_sigma:g}σ.</text>'
    )
    out.append("</svg>")
    path.write_text("\n".join(out), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--inactive", type=int, required=True)
    parser.add_argument("--cs-generations", type=int, required=True)
    parser.add_argument("--confidence-sigma", type=float, default=2.0)
    parser.add_argument("--expected-sign", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--source-tail-abs-pcm", type=float, default=1.0)
    parser.add_argument("--source-tail-relative-fraction", type=float, default=0.02)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    manifest_path = run_dir / "tasks.tsv"
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        tasks = list(csv.DictReader(handle, delimiter="\t"))
    cases: dict[int, dict[str, dict[str, str]]] = {}
    for task in tasks:
        cases.setdefault(int(task["case_index"]), {})[task["role"]] = task

    rows: list[dict[str, object]] = []
    for case_index in sorted(cases):
        case = cases[case_index]
        required = {"direct_ref", "direct_target", "cs"}
        if set(case) != required:
            raise ValueError(f"Case {case_index} has roles {sorted(case)}, expected {sorted(required)}")
        direct_ref = read_direct(Path(case["direct_ref"]["log_path"]))
        direct_target = read_direct(Path(case["direct_target"]["log_path"]))
        cs_task = case["cs"]
        cs = read_cs(
            Path(cs_task["log_path"]),
            Path(cs_task["work_dir"]),
            args.inactive,
            args.cs_generations,
        )
        step = cs["target_depth"] - cs["reference_depth"]
        direct_delta_k = direct_target["k"] - direct_ref["k"]
        direct_delta_k_stderr = quadrature(direct_ref["k_stderr"], direct_target["k_stderr"])
        delta_difference = cs["delta_k"] - direct_delta_k
        delta_agreement_stderr = quadrature(cs["delta_k_stderr"], direct_delta_k_stderr)
        delta_z = z_score(delta_difference, delta_agreement_stderr)
        target_difference = cs["predicted_target_k"] - direct_target["k"]
        target_agreement_stderr = quadrature(
            cs["predicted_target_k_stderr"], direct_target["k_stderr"]
        )
        target_z = z_score(target_difference, target_agreement_stderr)

        rho_ref = 1.0 - 1.0 / direct_ref["k"]
        rho_target = 1.0 - 1.0 / direct_target["k"]
        direct_delta_rho_pcm = (rho_target - rho_ref) * 1.0e5
        direct_delta_rho_stderr_pcm = quadrature(
            direct_ref["k_stderr"] / direct_ref["k"] ** 2,
            direct_target["k_stderr"] / direct_target["k"] ** 2,
        ) * 1.0e5
        rho_z = z_score(
            cs["delta_rho_pcm"] - direct_delta_rho_pcm,
            quadrature(cs["delta_rho_stderr_pcm"], direct_delta_rho_stderr_pcm),
        )

        source_tail_upper_pcm = (
            abs(cs["last_increment"])
            + args.confidence_sigma * cs["last_increment_stderr"]
        ) * 1.0e5
        source_tail_limit_pcm = max(
            args.source_tail_abs_pcm,
            args.source_tail_relative_fraction * abs(cs["delta_k"]) * 1.0e5,
        )
        source_converged = source_tail_upper_pcm <= source_tail_limit_pcm
        direct_sign = sign_assessment(
            direct_delta_k, direct_delta_k_stderr, args.expected_sign, args.confidence_sigma
        )
        cs_sign = sign_assessment(
            cs["delta_k"], cs["delta_k_stderr"], args.expected_sign, args.confidence_sigma
        )
        status = "PASS"
        if cs["diagnostics"] != 0:
            status = "INVALID_CS_DIAGNOSTIC"
        elif not source_converged:
            status = "SOURCE_TAIL_NOT_CONVERGED"
        elif direct_sign == "RESOLVED_WRONG" or cs_sign == "RESOLVED_WRONG":
            status = "MONOTONICITY_FAIL"
        elif abs(delta_z) > args.confidence_sigma:
            status = "DELTA_K_DISAGREE"
        elif abs(target_z) > args.confidence_sigma:
            status = "TARGET_K_DISAGREE"

        rows.append(
            {
                "case_index": case_index,
                "reference_depth_cm": cs["reference_depth"],
                "target_depth_cm": cs["target_depth"],
                "probe_step_cm": step,
                "direct_reference_seed": case["direct_ref"]["seed"],
                "direct_target_seed": case["direct_target"]["seed"],
                "cs_seed": case["cs"]["seed"],
                "direct_reference_k": direct_ref["k"],
                "direct_reference_k_stderr": direct_ref["k_stderr"],
                "direct_target_k": direct_target["k"],
                "direct_target_k_stderr": direct_target["k_stderr"],
                "cs_paired_reference_k": cs["paired_reference_k"],
                "cs_predicted_target_k": cs["predicted_target_k"],
                "cs_predicted_target_k_stderr": cs["predicted_target_k_stderr"],
                "predicted_minus_direct_target_k": target_difference,
                "predicted_target_k_agreement_stderr": target_agreement_stderr,
                "predicted_target_k_agreement_z": target_z,
                "predicted_target_k_agreement": agreement_label(target_z),
                "direct_delta_k": direct_delta_k,
                "direct_delta_k_stderr": direct_delta_k_stderr,
                "cs_delta_k": cs["delta_k"],
                "cs_delta_k_stderr": cs["delta_k_stderr"],
                "cs_minus_direct_delta_k": delta_difference,
                "delta_k_agreement_stderr": delta_agreement_stderr,
                "delta_k_agreement_z": delta_z,
                "delta_k_agreement": agreement_label(delta_z),
                "direct_slope_per_cm": direct_delta_k / step,
                "direct_slope_stderr_per_cm": direct_delta_k_stderr / step,
                "cs_slope_per_cm": cs["delta_k"] / step,
                "cs_slope_stderr_per_cm": cs["delta_k_stderr"] / step,
                "direct_delta_rho_pcm": direct_delta_rho_pcm,
                "direct_delta_rho_stderr_pcm": direct_delta_rho_stderr_pcm,
                "cs_delta_rho_pcm": cs["delta_rho_pcm"],
                "cs_delta_rho_stderr_pcm": cs["delta_rho_stderr_pcm"],
                "delta_rho_agreement_z": rho_z,
                "direct_sign_assessment": direct_sign,
                "cs_sign_assessment": cs_sign,
                "last_source_increment_pcm": abs(cs["last_increment"]) * 1.0e5,
                "source_tail_confidence_upper_pcm": source_tail_upper_pcm,
                "source_tail_limit_pcm": source_tail_limit_pcm,
                "source_tail_converged": source_converged,
                "cs_support_diagnostic_count": cs["diagnostics"],
                "status": status,
                "direct_reference_log": case["direct_ref"]["log_path"],
                "direct_target_log": case["direct_target"]["log_path"],
                "cs_log": case["cs"]["log_path"],
            }
        )

    output_csv = run_dir / "validation.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    passed = sum(row["status"] == "PASS" for row in rows)
    summary = [
        "# CS prediction validation against independent direct Monte Carlo",
        "",
        f"- Cases: {len(rows)} reference depths; {passed}/{len(rows)} passed all enabled checks",
        f"- Confidence threshold: {args.confidence_sigma:g} sigma",
        f"- Source propagation: L{args.cs_generations}",
        "- Each row uses three independent odd-seed calculations: direct MC at d, direct MC at d+h, and CS from d to d+h.",
        "",
        "Two complementary comparisons are reported. `delta z` compares the paired CS delta-k with the difference of the two direct MC k values. `target z` compares the standalone CS prediction k*(d+h) with the independent direct MC value at d+h.",
        "",
        "| d (cm) | d+h (cm) | direct slope (/cm) | CS slope (/cm) | delta z | target z | source tail | status |",
        "|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        tail = "pass" if row["source_tail_converged"] else "fail"
        summary.append(
            "| {reference_depth_cm:.6f} | {target_depth_cm:.6f} | "
            "{direct_slope_per_cm:.7e} ± {direct_slope_stderr_per_cm:.2e} | "
            "{cs_slope_per_cm:.7e} ± {cs_slope_stderr_per_cm:.2e} | "
            "{delta_k_agreement_z:.2f} | {predicted_target_k_agreement_z:.2f} | {tail} | {status} |".format(
                tail=tail, **row
            )
        )
    summary.extend(
        [
            "",
            "Interpretation: agreement within ±2 sigma is statistically compatible, not proof of exact equality. A source-tail failure is kept separate because an apparently good direct comparison can occur by chance while L10 has not stabilized.",
            "",
            "The direct finite-difference slope contains the independent uncertainties of both endpoint calculations. The CS slope uses the paired delta-k uncertainty, which preserves the correlation benefit of CS.",
        ]
    )
    (run_dir / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    svg_plot(rows, run_dir / "cs_validation.svg", args.confidence_sigma)

    print(f"Postprocessing complete: {passed}/{len(rows)} cases passed")
    print(f"  {output_csv}")
    print(f"  {run_dir / 'summary.md'}")
    print(f"  {run_dir / 'cs_validation.svg'}")


if __name__ == "__main__":
    main()

