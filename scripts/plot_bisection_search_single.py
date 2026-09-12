from __future__ import annotations

import argparse
import csv
from pathlib import Path

from PIL import Image, ImageDraw

from plot_cs_newton_search_single import (
    COLORS,
    arrowhead,
    centered_text,
    dashed_line,
    load_font,
    scale,
    vertical_text,
)


def load_search_rows(path: Path, first_search_evaluation: int) -> list[dict[str, float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        source = list(csv.DictReader(handle))
    rows = []
    for row in source:
        evaluation = int(row["evaluation"])
        if evaluation < first_search_evaluation or row["status"].upper() != "OK":
            continue
        rows.append(
            {
                "iteration": float(len(rows) + 1),
                "source_evaluation": float(evaluation),
                "depth": float(row["requested_depth_cm"]),
                "k": float(row["k_mean"]),
            }
        )
    if not rows:
        raise ValueError(f"No bisection-search rows found in {path}")
    return rows


def build_plot(rows: list[dict[str, float]], output: Path) -> None:
    width, height = 2000, 1125
    image = Image.new("RGB", (width, height), COLORS["white"])
    draw = ImageDraw.Draw(image)

    title_font = load_font(52, bold=True)
    axis_font = load_font(30, bold=True)
    tick_font = load_font(25)
    label_font = load_font(25)
    point_font = load_font(21, bold=True)
    target_font = load_font(25, bold=True)
    note_font = load_font(22, bold=True)

    draw.text((105, 55), "Bisection Critical Rod Position Search", font=title_font, fill=COLORS["ink"])

    left, top, right, bottom = 175, 190, 1870, 930
    iterations = [row["iteration"] for row in rows]
    k_values = [row["k"] for row in rows]
    x_lo, x_hi = min(iterations) - 0.35, max(iterations) + 0.35
    y_lo = min(k_values) - 0.007
    y_hi = max(k_values) + 0.008
    x_map = lambda value: scale(value, x_lo, x_hi, left, right)
    y_map = lambda value: scale(value, y_lo, y_hi, bottom, top)

    y_ticks = [0.98, 1.00, 1.02, 1.04, 1.06, 1.08]
    for value in y_ticks:
        if not y_lo <= value <= y_hi:
            continue
        y = y_map(value)
        draw.line([(left, y), (right, y)], fill=COLORS["grid"], width=2)
        label = f"{value:.2f}"
        box = draw.textbbox((0, 0), label, font=tick_font)
        draw.text((left - 20 - (box[2] - box[0]), y - 14), label, font=tick_font, fill=COLORS["muted"])

    for value in iterations:
        x = x_map(value)
        draw.line([(x, top), (x, bottom)], fill=COLORS["grid"], width=2)
        label = str(int(value))
        box = draw.textbbox((0, 0), label, font=tick_font)
        draw.text((x - (box[2] - box[0]) / 2, bottom + 18), label, font=tick_font, fill=COLORS["muted"])

    draw.line([(left, top), (left, bottom), (right, bottom)], fill=COLORS["axis"], width=4)
    centered_text(draw, ((left + right) / 2, 1020), "Bisection iteration", axis_font, COLORS["ink"])
    vertical_text(image, (62, (top + bottom) / 2), "k-effective", axis_font, COLORS["ink"])

    target_y = y_map(1.0)
    dashed_line(draw, (left, target_y), (right, target_y), COLORS["red"], width=4)
    target_label = "Target  k = 1.000000"
    target_x = left + 18
    target_box = draw.textbbox((0, 0), target_label, font=target_font)
    draw.rounded_rectangle(
        (target_x - 10, target_y - 47, target_x + target_box[2] - target_box[0] + 10, target_y - 10),
        radius=9,
        fill=COLORS["white"],
    )
    draw.text((target_x, target_y - 45), target_label, font=target_font, fill=COLORS["red"])

    points = [(x_map(row["iteration"]), y_map(row["k"])) for row in rows]
    for start, end in zip(points[:-1], points[1:]):
        draw.line([start, end], fill=COLORS["blue"], width=7)
        arrowhead(draw, start, end, COLORS["blue"])

    best_index = min(range(len(rows)), key=lambda index: abs(rows[index]["k"] - 1.0))
    for index, (point, row) in enumerate(zip(points, rows)):
        x, y = point
        best = index == best_index
        marker_color = COLORS["green"] if best else COLORS["blue"]
        radius = 25 if best else 22
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=marker_color,
            outline=COLORS["white"],
            width=5,
        )
        centered_text(draw, (x, y), str(int(row["iteration"])), point_font, COLORS["white"])

        if index == 0:
            offset_y = 44
        elif index == 1:
            offset_y = 38
        elif index in (2, 3, 4):
            offset_y = -62
        else:
            offset_y = -62 if index % 2 == 0 else 36
        depth_label = f'{row["depth"]:.3f} cm'
        box = draw.textbbox((0, 0), depth_label, font=label_font)
        label_x = max(left, min(x - (box[2] - box[0]) / 2, right - (box[2] - box[0])))
        draw.text(
            (label_x, y + offset_y),
            depth_label,
            font=label_font,
            fill=marker_color if best else COLORS["ink"],
        )

    best_x, best_y = points[best_index]
    note = "Closest sampled point"
    note_box = draw.textbbox((0, 0), note, font=note_font)
    note_x = best_x - (note_box[2] - note_box[0]) / 2
    note_y = best_y - 108
    draw.rounded_rectangle(
        (note_x - 14, note_y - 7, note_x + note_box[2] - note_box[0] + 14, note_y + 34),
        radius=11,
        fill="#DDF3E9",
    )
    draw.text((note_x, note_y), note, font=note_font, fill=COLORS["green"])

    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output, format="PNG", optimize=True, dpi=(200, 200))


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a single bisection search-path graph.")
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--first-search-evaluation", type=int, default=4)
    args = parser.parse_args()
    build_plot(load_search_rows(args.summary, args.first_search_evaluation), args.output)


if __name__ == "__main__":
    main()
