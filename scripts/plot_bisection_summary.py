from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from PIL import Image, ImageDraw

from plot_newton_history import COLORS, draw_axes, font, scale, text_center


def load_summary(path: Path) -> list[dict[str, float | str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        source = list(csv.DictReader(handle))
    if not source:
        raise ValueError(f"No bisection evaluations found in {path}")

    rows: list[dict[str, float | str]] = []
    for row in source:
        rows.append({
            "evaluation": float(row["evaluation"]),
            "depth": float(row["requested_depth_cm"]),
            "actual_depth": float(row["actual_depth_cm"]),
            "k": float(row["k_mean"]),
            "stderr": float(row["k_stderr_conservative"]),
            "status": row["status"],
        })
    return rows


def load_result(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        for match in re.finditer(r"([A-Za-z0-9_]+)=([^\s]+)", line):
            values[match.group(1)] = match.group(2)
    return values


def reconstruct_brackets(rows: list[dict[str, float | str]], target: float):
    lower = None
    upper = None
    history: list[tuple[int, float, float]] = []
    for row in rows:
        depth = float(row["depth"])
        k_value = float(row["k"])
        evaluation = int(float(row["evaluation"]))
        if k_value >= target:
            lower = depth if lower is None else max(lower, depth)
        else:
            upper = depth if upper is None else min(upper, depth)
        if lower is not None and upper is not None:
            history.append((evaluation, lower, upper))
    return history


def build_plot(rows: list[dict[str, float | str]], result: dict[str, str], output: Path) -> None:
    image = Image.new("RGB", (1800, 1200), COLORS["white"])
    draw = ImageDraw.Draw(image)
    title = font(42, bold=True)
    subtitle = font(24)
    panel_title = font(27, bold=True)
    body = font(23)
    body_bold = font(23, bold=True)
    small = font(20)

    target = float(result.get("target_k", 1.0))
    evaluations = [int(float(row["evaluation"])) for row in rows]
    depths = [float(row["depth"]) for row in rows]
    k_values = [float(row["k"]) for row in rows]
    k_errors = [float(row["stderr"]) for row in rows]
    residuals = [abs(value - target) * 1.0e5 for value in k_values]
    errors_pcm = [value * 1.0e5 for value in k_errors]

    bracket_history = reconstruct_brackets(rows, target)
    if not bracket_history:
        raise ValueError("Could not reconstruct a valid supercritical/subcritical bracket")

    final_lower = float(result.get("final_lower_depth_cm", bracket_history[-1][1]))
    final_upper = float(result.get("final_upper_depth_cm", bracket_history[-1][2]))
    estimated_depth = float(result.get("estimated_critical_depth_cm", 0.5 * (final_lower + final_upper)))
    best_depth = float(result.get("best_sampled_depth_cm", depths[-1]))
    best_k = float(result.get("best_sampled_k", k_values[-1]))
    best_stderr = float(result.get("best_sampled_k_stderr", k_errors[-1]))
    best_residual = (best_k - target) * 1.0e5
    final_width = final_upper - final_lower

    draw.text((90, 55), "Bisection Critical Rod Position Search", font=title, fill=COLORS["ink"])
    draw.text((90, 110), "Independent k-effective evaluations   |   Target k-effective = 1.000000",
              font=subtitle, fill=COLORS["muted"])

    main_box = (135, 205, 1260, 675)
    draw.text((main_box[0], 160), "k-effective search sequence", font=panel_title, fill=COLORS["ink"])
    y_lo = min(k_values) - 0.0008
    y_hi = max(k_values) + 0.0008
    y_ticks = [0.998, 1.000, 1.002, 1.004, 1.006, 1.008, 1.010]
    x_map = lambda value: scale(value, evaluations[0], evaluations[-1], main_box[0], main_box[2])
    y_map = lambda value: scale(value, y_lo, y_hi, main_box[3], main_box[1])
    draw_axes(image, draw, main_box, evaluations, y_ticks, x_map, y_map,
              lambda value: f"{int(value)}", lambda value: f"{value:.3f}",
              "Bisection evaluation", "k-effective")

    target_y = y_map(target)
    draw.line([(main_box[0], target_y), (main_box[2], target_y)], fill=COLORS["red"], width=4)
    draw.rounded_rectangle((1010, target_y + 8, 1248, target_y + 38), radius=8,
                           fill=COLORS["white"], outline=COLORS["red"], width=2)
    text_center(draw, (1129, target_y + 23), "Critical target  k = 1", small, COLORS["red"])

    points = [(x_map(ev), y_map(kval)) for ev, kval in zip(evaluations, k_values)]
    draw.line(points, fill=COLORS["blue"], width=5)
    best_index = min(range(len(rows)), key=lambda i: abs(depths[i] - best_depth))
    for index, ((x, y), depth, kval, error) in enumerate(zip(points, depths, k_values, k_errors)):
        error_pixels = abs(y_map(kval + error) - y_map(kval))
        if error_pixels >= 1.0:
            draw.line([(x, y - error_pixels), (x, y + error_pixels)], fill=COLORS["ink"], width=2)
            draw.line([(x - 6, y - error_pixels), (x + 6, y - error_pixels)], fill=COLORS["ink"], width=2)
            draw.line([(x - 6, y + error_pixels), (x + 6, y + error_pixels)], fill=COLORS["ink"], width=2)
        marker_color = COLORS["green"] if index == best_index else COLORS["blue"]
        radius = 18 if index == best_index else 15
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=marker_color,
                     outline=COLORS["white"], width=4)
        text_center(draw, (x, y), str(evaluations[index]), font(17, bold=True), COLORS["white"])
        offset_y = -42 if index % 2 == 0 else 26
        if index == 0:
            offset_y = -42
        if index == best_index:
            offset_y = -48
        label = f"{depth:.4f} cm"
        bounds = draw.textbbox((0, 0), label, font=small)
        label_x = max(main_box[0], min(x - (bounds[2] - bounds[0]) / 2, main_box[2] - (bounds[2] - bounds[0])))
        draw.text((label_x, y + offset_y), label, font=small,
                  fill=COLORS["green"] if index == best_index else COLORS["ink"])

    card = (1310, 160, 1710, 675)
    draw.rounded_rectangle(card, radius=22, fill=COLORS["card"], outline=COLORS["grid"], width=2)
    draw.text((1350, 202), "Converged solution", font=panel_title, fill=COLORS["ink"])
    draw.text((1350, 266), f"{estimated_depth:.5f} cm", font=font(38, bold=True), fill=COLORS["green"])
    draw.text((1350, 326), f"Bracket: {final_lower:.6f}", font=body, fill=COLORS["ink"])
    draw.text((1350, 362), f"to {final_upper:.6f} cm", font=body, fill=COLORS["ink"])
    draw.text((1350, 402), f"Width = {final_width:.5f} cm", font=body, fill=COLORS["muted"])
    draw.line([(1350, 448), (1670, 448)], fill=COLORS["grid"], width=2)
    draw.text((1350, 478), f"Best sampled k = {best_k:.7f}", font=body_bold, fill=COLORS["ink"])
    draw.text((1350, 516), f"Residual = {best_residual:+.2f} pcm", font=body, fill=COLORS["muted"])
    draw.text((1350, 552), f"1 sigma = {best_stderr * 1e5:.2f} pcm", font=body, fill=COLORS["muted"])
    draw.text((1350, 592), "Stop: bracket width <= 0.01 cm", font=small, fill=COLORS["ink"])
    draw.rounded_rectangle((1350, 624, 1540, 658), radius=10, fill=COLORS["green_light"])
    text_center(draw, (1445, 641), "INSIDE 1 SIGMA", font(18, bold=True), COLORS["green"])

    bracket_box = (135, 795, 830, 1085)
    draw.text((bracket_box[0], 746), "Bisection bracket", font=panel_title, fill=COLORS["ink"])
    b_evals = [item[0] for item in bracket_history]
    lowers = [item[1] for item in bracket_history]
    uppers = [item[2] for item in bracket_history]
    bx_map = lambda value: scale(value, b_evals[0], b_evals[-1], bracket_box[0], bracket_box[2])
    bd_map = lambda value: scale(value, 81.9, 84.1, bracket_box[3], bracket_box[1])
    draw_axes(image, draw, bracket_box, b_evals, [82.0, 82.5, 83.0, 83.5, 84.0], bx_map, bd_map,
              lambda value: f"{int(value)}", lambda value: f"{value:.1f}",
              "Bisection evaluation", "")
    lower_points = [(bx_map(ev), bd_map(value)) for ev, value in zip(b_evals, lowers)]
    upper_points = [(bx_map(ev), bd_map(value)) for ev, value in zip(b_evals, uppers)]
    band = lower_points + list(reversed(upper_points))
    draw.polygon(band, fill=COLORS["blue_light"])
    draw.line(lower_points, fill=COLORS["green"], width=5)
    draw.line(upper_points, fill=COLORS["red"], width=5)
    estimate_y = bd_map(estimated_depth)
    draw.line([(bracket_box[0], estimate_y), (bracket_box[2], estimate_y)], fill=COLORS["ink"], width=3)
    for point in lower_points:
        draw.ellipse((point[0] - 6, point[1] - 6, point[0] + 6, point[1] + 6), fill=COLORS["green"])
    for point in upper_points:
        draw.ellipse((point[0] - 6, point[1] - 6, point[0] + 6, point[1] + 6), fill=COLORS["red"])
    draw.line([(455, 762), (495, 762)], fill=COLORS["green"], width=5)
    draw.text((505, 750), "lower", font=small, fill=COLORS["ink"])
    draw.line([(585, 762), (625, 762)], fill=COLORS["red"], width=5)
    draw.text((635, 750), "upper", font=small, fill=COLORS["ink"])
    draw.line([(715, 762), (755, 762)], fill=COLORS["ink"], width=3)
    draw.text((765, 750), "estimate", font=small, fill=COLORS["ink"])

    residual_box = (1005, 795, 1710, 1085)
    draw.text((residual_box[0], 746), "Residual versus 1 sigma", font=panel_title, fill=COLORS["ink"])
    r_map = lambda value: scale(math.log10(max(value, 0.1)), -1.0, 4.0, residual_box[3], residual_box[1])
    rx_map = lambda value: scale(value, evaluations[0], evaluations[-1], residual_box[0], residual_box[2])
    draw_axes(image, draw, residual_box, evaluations, [0.1, 1, 10, 100, 1000, 10000], rx_map, r_map,
              lambda value: f"{int(value)}", lambda value: "0.1" if value < 1 else f"{int(value):,}",
              "Bisection evaluation", "")
    residual_points = [(rx_map(ev), r_map(value)) for ev, value in zip(evaluations, residuals)]
    error_points = [(rx_map(ev), r_map(value)) for ev, value in zip(evaluations, errors_pcm)]
    draw.line(residual_points, fill=COLORS["orange"], width=5)
    draw.line(error_points, fill=COLORS["blue"], width=4)
    for point in residual_points:
        draw.ellipse((point[0] - 8, point[1] - 8, point[0] + 8, point[1] + 8), fill=COLORS["orange"])
    for point in error_points:
        draw.ellipse((point[0] - 6, point[1] - 6, point[0] + 6, point[1] + 6), fill=COLORS["blue"])
    best_point = residual_points[best_index]
    draw.ellipse((best_point[0] - 14, best_point[1] - 14, best_point[0] + 14, best_point[1] + 14),
                 outline=COLORS["green"], width=4)
    draw.text((best_point[0] - 145, best_point[1] - 42), "inside 1 sigma", font=small, fill=COLORS["green"])
    draw.line([(1375, 762), (1418, 762)], fill=COLORS["orange"], width=5)
    draw.text((1430, 750), "|k - 1|", font=small, fill=COLORS["ink"])
    draw.line([(1535, 762), (1578, 762)], fill=COLORS["blue"], width=4)
    draw.text((1590, 750), "1 sigma", font=small, fill=COLORS["ink"])

    draw.text((90, 1152), "Source: summary.csv and result.txt", font=font(18), fill=COLORS["muted"])
    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output, format="PNG", optimize=True, dpi=(200, 200))


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot a bisection critical-depth search summary.")
    parser.add_argument("summary", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build_plot(load_summary(args.summary), load_result(args.result), args.output)


if __name__ == "__main__":
    main()
