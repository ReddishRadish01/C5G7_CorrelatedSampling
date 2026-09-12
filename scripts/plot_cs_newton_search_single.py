from __future__ import annotations

import argparse
import csv
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


COLORS = {
    "ink": "#183247",
    "muted": "#61748F",
    "axis": "#8DA1B5",
    "grid": "#D8E2EC",
    "blue": "#2F68B8",
    "green": "#16875E",
    "red": "#D14747",
    "white": "#FFFFFF",
}


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    filename = "arialbd.ttf" if bold else "arial.ttf"
    for candidate in (Path("C:/Windows/Fonts") / filename, Path(filename)):
        try:
            return ImageFont.truetype(str(candidate), size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def find_column(fieldnames: list[str], *candidates: str) -> str:
    normalized = {name.lower().replace("_", ""): name for name in fieldnames}
    for candidate in candidates:
        key = candidate.lower().replace("_", "")
        if key in normalized:
            return normalized[key]
    raise ValueError(f"Missing required column; expected one of {candidates}")


def load_history(path: Path) -> list[dict[str, float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        iteration_col = find_column(reader.fieldnames, "iteration")
        depth_col = find_column(reader.fieldnames, "requested_depth_cm", "RequestedDepthCm")
        k_col = find_column(reader.fieldnames, "k_current", "KCurrent")
        rows = [
            {
                "iteration": float(row[iteration_col]),
                "depth": float(row[depth_col]),
                "k": float(row[k_col]),
            }
            for row in reader
            if (row.get(iteration_col) or "").strip()
        ]
    if not rows:
        raise ValueError(f"No Newton iterations found in {path}")
    return rows


def scale(value: float, lo: float, hi: float, pixel_lo: float, pixel_hi: float) -> float:
    return pixel_lo + (value - lo) * (pixel_hi - pixel_lo) / (hi - lo)


def centered_text(
    draw: ImageDraw.ImageDraw,
    center: tuple[float, float],
    label: str,
    font: ImageFont.ImageFont,
    fill: str,
) -> None:
    box = draw.textbbox((0, 0), label, font=font)
    draw.text(
        (center[0] - (box[2] - box[0]) / 2, center[1] - (box[3] - box[1]) / 2 - box[1]),
        label,
        font=font,
        fill=fill,
    )


def vertical_text(
    image: Image.Image,
    center: tuple[float, float],
    label: str,
    font: ImageFont.ImageFont,
    fill: str,
) -> None:
    scratch = Image.new("RGBA", (500, 100), (0, 0, 0, 0))
    scratch_draw = ImageDraw.Draw(scratch)
    box = scratch_draw.textbbox((0, 0), label, font=font)
    crop = Image.new("RGBA", (box[2] - box[0] + 16, box[3] - box[1] + 16), (0, 0, 0, 0))
    ImageDraw.Draw(crop).text((8 - box[0], 8 - box[1]), label, font=font, fill=fill)
    rotated = crop.rotate(90, expand=True)
    image.paste(
        rotated,
        (int(center[0] - rotated.width / 2), int(center[1] - rotated.height / 2)),
        rotated,
    )


def dashed_line(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    fill: str,
    width: int,
    dash: int = 18,
    gap: int = 12,
) -> None:
    x1, y1 = start
    x2, y2 = end
    length = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
    if length == 0:
        return
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    position = 0.0
    while position < length:
        segment_end = min(position + dash, length)
        draw.line(
            [
                (x1 + ux * position, y1 + uy * position),
                (x1 + ux * segment_end, y1 + uy * segment_end),
            ],
            fill=fill,
            width=width,
        )
        position += dash + gap


def arrowhead(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    fill: str,
) -> None:
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = (dx * dx + dy * dy) ** 0.5
    if length < 1:
        return
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    cx = start[0] + dx * 0.72
    cy = start[1] + dy * 0.72
    size = 15
    draw.polygon(
        [
            (cx + size * ux, cy + size * uy),
            (cx - size * ux + 8 * px, cy - size * uy + 8 * py),
            (cx - size * ux - 8 * px, cy - size * uy - 8 * py),
        ],
        fill=fill,
    )


def build_plot(rows: list[dict[str, float]], output: Path) -> None:
    width, height = 2000, 1125
    image = Image.new("RGB", (width, height), COLORS["white"])
    draw = ImageDraw.Draw(image)

    title_font = load_font(52, bold=True)
    axis_font = load_font(30, bold=True)
    tick_font = load_font(25)
    label_font = load_font(27)
    point_font = load_font(23, bold=True)
    target_font = load_font(25, bold=True)

    draw.text((105, 55), "CS Newton-Raphson Critical Rod Position Search", font=title_font, fill=COLORS["ink"])

    left, top, right, bottom = 175, 190, 1870, 930
    iterations = [row["iteration"] for row in rows]
    k_values = [row["k"] for row in rows]
    depths = [row["depth"] for row in rows]

    x_lo, x_hi = min(iterations) - 0.35, max(iterations) + 0.35
    y_lo = min(0.925, min(k_values) - 0.008)
    y_hi = max(1.095, max(k_values) + 0.008)
    x_map = lambda value: scale(value, x_lo, x_hi, left, right)
    y_map = lambda value: scale(value, y_lo, y_hi, bottom, top)

    y_ticks = [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06, 1.08]
    for value in y_ticks:
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
    centered_text(draw, ((left + right) / 2, 1020), "Newton iteration", axis_font, COLORS["ink"])
    vertical_text(image, (62, (top + bottom) / 2), "k-effective", axis_font, COLORS["ink"])

    target_y = y_map(1.0)
    dashed_line(draw, (left, target_y), (right, target_y), COLORS["red"], width=4)
    target_label = "Target  k = 1.000000"
    target_box = draw.textbbox((0, 0), target_label, font=target_font)
    # Keep the target label away from the final depth annotation.
    target_x = x_map(iterations[-1] - 0.9)
    draw.rounded_rectangle(
        (target_x - 12, target_y - 47, right - 4, target_y - 10),
        radius=9,
        fill=COLORS["white"],
    )
    draw.text((target_x, target_y - 45), target_label, font=target_font, fill=COLORS["red"])

    points = [(x_map(row["iteration"]), y_map(row["k"])) for row in rows]
    for start, end in zip(points[:-1], points[1:]):
        draw.line([start, end], fill=COLORS["blue"], width=7)
        arrowhead(draw, start, end, COLORS["blue"])

    label_offsets = [(0, 48), (0, -60), (0, -60), (-55, 50), (0, -64)]
    if len(rows) != 5:
        label_offsets = [(0, -58 if index % 2 == 0 else 48) for index in range(len(rows))]

    for index, (point, row, offset) in enumerate(zip(points, rows, label_offsets)):
        x, y = point
        final = index == len(rows) - 1
        marker_color = COLORS["green"] if final else COLORS["blue"]
        radius = 25 if final else 22
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=marker_color,
            outline=COLORS["white"],
            width=5,
        )
        centered_text(draw, (x, y), str(int(row["iteration"])), point_font, COLORS["white"])
        depth_label = f'{row["depth"]:.3f} cm'
        box = draw.textbbox((0, 0), depth_label, font=label_font)
        label_x = x + offset[0] - (box[2] - box[0]) / 2
        label_y = y + offset[1]
        draw.text((label_x, label_y), depth_label, font=label_font, fill=marker_color if final else COLORS["ink"])

    final_x, final_y = points[-1]
    draw.rounded_rectangle(
        (final_x + 40, final_y + 34, final_x + 230, final_y + 78),
        radius=12,
        fill="#DDF3E9",
    )
    centered_text(draw, (final_x + 135, final_y + 56), "Converged", target_font, COLORS["green"])

    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output, format="PNG", optimize=True, dpi=(200, 200))


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a single CS Newton search-path graph.")
    parser.add_argument("history", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build_plot(load_history(args.history), args.output)


if __name__ == "__main__":
    main()
