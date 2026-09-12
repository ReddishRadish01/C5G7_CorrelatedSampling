from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


COLORS = {
    "ink": "#183247",
    "muted": "#64748B",
    "grid": "#D9E2EA",
    "axis": "#8CA0B3",
    "blue": "#2563B8",
    "blue_light": "#DCEBFA",
    "orange": "#E88A17",
    "green": "#17845B",
    "green_light": "#DDF3E9",
    "red": "#C84545",
    "card": "#F5F8FB",
    "white": "#FFFFFF",
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "arialbd.ttf" if bold else "arial.ttf"
    candidates = [Path("C:/Windows/Fonts") / name, Path(name)]
    for candidate in candidates:
        try:
            return ImageFont.truetype(str(candidate), size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def load_history(path: Path) -> list[dict[str, object]]:
    numeric = {
        "Iteration",
        "RequestedDepthCm",
        "KCurrent",
        "KStderr",
        "ResidualPcm",
        "RawNewtonDepthCm",
        "NextDepthCm",
    }
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No Newton iterations found in {path}")

    parsed: list[dict[str, object]] = []
    for row in rows:
        item: dict[str, object] = dict(row)
        for key in numeric:
            value = (row.get(key) or "").strip()
            item[key] = float(value) if value else math.nan
        parsed.append(item)
    return parsed


def scale(value: float, lo: float, hi: float, pixel_lo: float, pixel_hi: float) -> float:
    if hi == lo:
        return 0.5 * (pixel_lo + pixel_hi)
    return pixel_lo + (value - lo) * (pixel_hi - pixel_lo) / (hi - lo)


def text_center(draw: ImageDraw.ImageDraw, xy: tuple[float, float], label: str,
                face: ImageFont.ImageFont, fill: str) -> None:
    box = draw.textbbox((0, 0), label, font=face)
    width = box[2] - box[0]
    height = box[3] - box[1]
    draw.text((xy[0] - width / 2, xy[1] - height / 2 - box[1]), label, font=face, fill=fill)


def vertical_text(image: Image.Image, center: tuple[float, float], label: str,
                  face: ImageFont.ImageFont, fill: str) -> None:
    scratch = Image.new("RGBA", (1200, 120), (0, 0, 0, 0))
    scratch_draw = ImageDraw.Draw(scratch)
    box = scratch_draw.textbbox((0, 0), label, font=face)
    width = box[2] - box[0]
    height = box[3] - box[1]
    crop = Image.new("RGBA", (width + 12, height + 12), (0, 0, 0, 0))
    ImageDraw.Draw(crop).text((6 - box[0], 6 - box[1]), label, font=face, fill=fill)
    rotated = crop.rotate(90, expand=True)
    image.paste(rotated, (int(center[0] - rotated.width / 2), int(center[1] - rotated.height / 2)), rotated)


def arrow_segment(draw: ImageDraw.ImageDraw, start: tuple[float, float], end: tuple[float, float],
                  color: str, width: int = 5) -> None:
    draw.line([start, end], fill=color, width=width)
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if length < 1:
        return
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    tip = (end[0] - 10 * ux, end[1] - 10 * uy)
    left = (tip[0] - 14 * ux + 8 * px, tip[1] - 14 * uy + 8 * py)
    right = (tip[0] - 14 * ux - 8 * px, tip[1] - 14 * uy - 8 * py)
    draw.polygon([end, left, right], fill=color)


def draw_axes(image: Image.Image, draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int],
              x_ticks: list[float], y_ticks: list[float],
              x_map, y_map, x_formatter, y_formatter,
              x_label: str, y_label: str) -> None:
    left, top, right, bottom = box
    small = font(22)
    label_font = font(24, bold=True)
    for tick in y_ticks:
        y = y_map(tick)
        draw.line([(left, y), (right, y)], fill=COLORS["grid"], width=2)
        label = y_formatter(tick)
        bounds = draw.textbbox((0, 0), label, font=small)
        draw.text((left - 16 - (bounds[2] - bounds[0]), y - 12), label, font=small, fill=COLORS["muted"])
    for tick in x_ticks:
        x = x_map(tick)
        draw.line([(x, top), (x, bottom)], fill=COLORS["grid"], width=2)
        label = x_formatter(tick)
        bounds = draw.textbbox((0, 0), label, font=small)
        draw.text((x - (bounds[2] - bounds[0]) / 2, bottom + 12), label, font=small, fill=COLORS["muted"])
    draw.line([(left, top), (left, bottom), (right, bottom)], fill=COLORS["axis"], width=3)
    text_center(draw, ((left + right) / 2, bottom + 58), x_label, label_font, COLORS["ink"])
    if y_label:
        vertical_text(image, (left - 94, (top + bottom) / 2), y_label, label_font, COLORS["ink"])


def build_plot(rows: list[dict[str, object]], output: Path) -> None:
    image = Image.new("RGB", (1800, 1200), COLORS["white"])
    draw = ImageDraw.Draw(image)

    title = font(42, bold=True)
    subtitle = font(24)
    panel_title = font(27, bold=True)
    body = font(23)
    body_bold = font(23, bold=True)
    small = font(20)

    draw.text((90, 55), "Newton-Raphson Critical Rod Position Search", font=title, fill=COLORS["ink"])
    draw.text((90, 110), "Correlated-sampling derivative   |   Target k-effective = 1.000000",
              font=subtitle, fill=COLORS["muted"])

    depths = [float(row["RequestedDepthCm"]) for row in rows]
    k_values = [float(row["KCurrent"]) for row in rows]
    k_errors = [float(row["KStderr"]) for row in rows]
    iterations = [int(float(row["Iteration"])) for row in rows]
    residuals = [abs(float(row["ResidualPcm"])) for row in rows]
    k_errors_pcm = [value * 1.0e5 for value in k_errors]

    main_box = (135, 205, 1260, 675)
    draw.text((main_box[0], 160), "Search path in rod-depth space", font=panel_title, fill=COLORS["ink"])
    x_lo = math.floor(min(depths) - 2.0)
    x_hi = math.ceil(max(depths) + 2.0)
    y_lo = min(min(k_values) - 0.012, 0.925)
    y_hi = max(max(k_values) + 0.010, 1.095)
    x_ticks = [64, 70, 76, 82, 88, 94]
    y_ticks = [0.94, 0.97, 1.00, 1.03, 1.06, 1.09]
    x_map = lambda value: scale(value, x_lo, x_hi, main_box[0], main_box[2])
    y_map = lambda value: scale(value, y_lo, y_hi, main_box[3], main_box[1])
    draw_axes(image, draw, main_box, x_ticks, y_ticks, x_map, y_map,
              lambda value: f"{value:.0f}", lambda value: f"{value:.2f}",
              "Rod insertion depth (cm)", "k-effective")

    target_y = y_map(1.0)
    draw.line([(main_box[0], target_y), (main_box[2], target_y)], fill=COLORS["red"], width=4)
    draw.rounded_rectangle((1010, target_y + 7, 1248, target_y + 37), radius=8,
                           fill=COLORS["white"], outline=COLORS["red"], width=2)
    text_center(draw, (1129, target_y + 22), "Critical target  k = 1", small, COLORS["red"])

    points = [(x_map(depth), y_map(kval)) for depth, kval in zip(depths, k_values)]
    for start, end in zip(points[:-1], points[1:]):
        arrow_segment(draw, start, end, COLORS["blue"], width=5)

    for index, ((x, y), iteration, kval, error) in enumerate(zip(points, iterations, k_values, k_errors)):
        error_pixels = abs(y_map(kval + error) - y_map(kval))
        if error_pixels >= 1.0:
            draw.line([(x, y - error_pixels), (x, y + error_pixels)], fill=COLORS["ink"], width=2)
            draw.line([(x - 7, y - error_pixels), (x + 7, y - error_pixels)], fill=COLORS["ink"], width=2)
            draw.line([(x - 7, y + error_pixels), (x + 7, y + error_pixels)], fill=COLORS["ink"], width=2)
        marker_color = COLORS["green"] if index == len(points) - 1 else COLORS["blue"]
        radius = 20 if index == len(points) - 1 else 17
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=marker_color,
                     outline=COLORS["white"], width=4)
        text_center(draw, (x, y), str(iteration), font(19, bold=True), COLORS["white"])

    # Give the clustered root iterations distinct labels without obscuring the path.
    offsets = [(-8, -42), (-180, -42), (-100, -42), (45, 30), (70, -110)]
    for (x, y), depth, kval, offset in zip(points, depths, k_values, offsets):
        label = f"{depth:.3f} cm,  {kval:.6f}"
        draw.text((x + offset[0], y + offset[1]), label, font=small, fill=COLORS["ink"])

    card = (1310, 160, 1710, 675)
    draw.rounded_rectangle(card, radius=22, fill=COLORS["card"], outline=COLORS["grid"], width=2)
    draw.text((1350, 202), "Converged solution", font=panel_title, fill=COLORS["ink"])
    final = rows[-1]
    final_depth = float(final["RequestedDepthCm"])
    final_k = float(final["KCurrent"])
    final_se = float(final["KStderr"])
    final_residual = float(final["ResidualPcm"])
    draw.text((1350, 272), f"{final_depth:.5f} cm", font=font(38, bold=True), fill=COLORS["green"])
    draw.text((1350, 332), f"k = {final_k:.7f}", font=body_bold, fill=COLORS["ink"])
    draw.text((1350, 372), f"1 sigma = {final_se * 1e5:.2f} pcm", font=body, fill=COLORS["muted"])
    draw.text((1350, 412), f"Residual = {final_residual:.2f} pcm", font=body, fill=COLORS["muted"])
    draw.line([(1350, 468), (1670, 468)], fill=COLORS["grid"], width=2)
    draw.text((1350, 500), f"Iterations: {len(rows)}", font=body_bold, fill=COLORS["ink"])
    draw.text((1350, 544), "Stopping criterion:", font=body_bold, fill=COLORS["ink"])
    draw.text((1350, 580), "|k - 1| <= 1 sigma(k)", font=body, fill=COLORS["muted"])
    draw.rounded_rectangle((1350, 620, 1502, 655), radius=10, fill=COLORS["green_light"])
    text_center(draw, (1426, 637), str(final["Status"]), small, COLORS["green"])

    depth_box = (135, 795, 830, 1085)
    draw.text((depth_box[0], 746), "Rod-depth convergence (cm)", font=panel_title, fill=COLORS["ink"])
    d_lo = math.floor(min(depths) / 5.0) * 5.0
    d_hi = math.ceil(max(depths) / 5.0) * 5.0
    it_map = lambda value: scale(value, min(iterations), max(iterations), depth_box[0], depth_box[2])
    d_map = lambda value: scale(value, d_lo, d_hi, depth_box[3], depth_box[1])
    depth_ticks = list(range(int(d_lo), int(d_hi) + 1, 5))
    draw_axes(image, draw, depth_box, iterations, depth_ticks, it_map, d_map,
              lambda value: f"{int(value)}", lambda value: f"{int(value)}",
              "Newton iteration", "")
    depth_points = [(it_map(it), d_map(depth)) for it, depth in zip(iterations, depths)]
    draw.line(depth_points, fill=COLORS["blue"], width=5)
    for (x, y), depth in zip(depth_points, depths):
        draw.ellipse((x - 10, y - 10, x + 10, y + 10), fill=COLORS["blue"])
        draw.text((x + 12, y - 25), f"{depth:.2f}", font=small, fill=COLORS["ink"])
    final_y = d_map(final_depth)
    draw.line([(depth_box[0], final_y), (depth_box[2], final_y)], fill=COLORS["green"], width=3)

    residual_box = (1005, 795, 1710, 1085)
    draw.text((residual_box[0], 746), "Residual versus 1 sigma", font=panel_title, fill=COLORS["ink"])
    r_map = lambda value: scale(math.log10(max(value, 1.0)), 0.0, 4.0, residual_box[3], residual_box[1])
    r_it_map = lambda value: scale(value, min(iterations), max(iterations), residual_box[0], residual_box[2])
    draw_axes(image, draw, residual_box, iterations, [1, 10, 100, 1000, 10000], r_it_map, r_map,
              lambda value: f"{int(value)}", lambda value: f"{int(value):,}",
              "Newton iteration", "")
    residual_points = [(r_it_map(it), r_map(value)) for it, value in zip(iterations, residuals)]
    error_points = [(r_it_map(it), r_map(value)) for it, value in zip(iterations, k_errors_pcm)]
    draw.line(residual_points, fill=COLORS["orange"], width=5)
    draw.line(error_points, fill=COLORS["blue"], width=4)
    for point in residual_points:
        draw.ellipse((point[0] - 9, point[1] - 9, point[0] + 9, point[1] + 9), fill=COLORS["orange"])
    for point in error_points:
        draw.ellipse((point[0] - 7, point[1] - 7, point[0] + 7, point[1] + 7), fill=COLORS["blue"])
    final_residual_point = residual_points[-1]
    draw.ellipse((final_residual_point[0] - 14, final_residual_point[1] - 14,
                  final_residual_point[0] + 14, final_residual_point[1] + 14),
                 outline=COLORS["green"], width=4)
    draw.text((final_residual_point[0] - 150, final_residual_point[1] - 42),
              "inside 1 sigma", font=small, fill=COLORS["green"])
    draw.line([(1375, 762), (1418, 762)], fill=COLORS["orange"], width=5)
    draw.text((1430, 750), "|k - 1|", font=small, fill=COLORS["ink"])
    draw.line([(1535, 762), (1578, 762)], fill=COLORS["blue"], width=4)
    draw.text((1590, 750), "1 sigma", font=small, fill=COLORS["ink"])

    draw.text((90, 1152), "Source: newton_history_.csv", font=font(18), fill=COLORS["muted"])
    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output, format="PNG", optimize=True, dpi=(200, 200))


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot a correlated-sampling Newton search history.")
    parser.add_argument("history", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build_plot(load_history(args.history), args.output)


if __name__ == "__main__":
    main()
