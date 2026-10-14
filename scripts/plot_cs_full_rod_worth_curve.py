import argparse
import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]


def load_rows(paths):
    by_depth = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                depth = float(row["reference_depth_cm"])
                by_depth[depth] = {
                    "z": depth,
                    "k": float(row["cs_paired_reference_k"]),
                    "m": float(row["cs_slope_per_cm"]),
                    "m_se": float(row["cs_slope_stderr_per_cm"]),
                    "status": row["status"],
                    "target_z": float(row["target_depth_cm"]),
                    "target_k": float(row["cs_predicted_target_k"]),
                }
    rows = sorted(by_depth.values(), key=lambda item: item["z"])
    if not rows:
        raise RuntimeError("No validation rows were found.")

    # The final run contains the ARI target at 128.52 cm.  Add it as the final
    # curve endpoint when it extends beyond the last reference point.
    last = rows[-1]
    if last["target_z"] > last["z"]:
        rows.append({
            "z": last["target_z"],
            "k": last["target_k"],
            "m": last["m"],
            "m_se": last["m_se"],
            "status": last["status"],
            "target_z": last["target_z"],
            "target_k": last["target_k"],
            "endpoint_only": True,
        })
    return rows


def get_font(size, bold=False):
    name = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(str(Path("C:/Windows/Fonts") / name), size)


def hermite(a, b, t):
    az, ak, am = a["z"], a["k"], a["m"]
    bz, bk, bm = b["z"], b["k"], b["m"]
    h = bz - az
    t2, t3 = t * t, t * t * t
    value = ((2 * t3 - 3 * t2 + 1) * ak
             + (t3 - 2 * t2 + t) * h * am
             + (-2 * t3 + 3 * t2) * bk
             + (t3 - t2) * h * bm)
    derivative = (((6 * t2 - 6 * t) / h) * ak
                  + (3 * t2 - 4 * t + 1) * am
                  + ((-6 * t2 + 6 * t) / h) * bk
                  + (3 * t2 - 2 * t) * bm)
    return value, derivative


def reconstruct(rows):
    curve = []
    for index in range(len(rows) - 1):
        steps = max(40, int((rows[index + 1]["z"] - rows[index]["z"]) * 35))
        stop = steps + 1 if index == len(rows) - 2 else steps
        for step in range(stop):
            t = step / steps
            kval, slope = hermite(rows[index], rows[index + 1], t)
            zval = rows[index]["z"] + t * (rows[index + 1]["z"] - rows[index]["z"])
            worth = -slope / (kval * kval) * 1.0e5
            curve.append((zval, kval, slope, worth))
    return curve


def rounded_ceiling(value, quantum):
    return math.ceil(value / quantum) * quantum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = load_rows(args.csv)
    curve = reconstruct(rows)
    actual_rows = [row for row in rows if not row.get("endpoint_only")]

    critical_z = None
    for a, b in zip(curve, curve[1:]):
        if (a[1] - 1.0) * (b[1] - 1.0) <= 0 and b[1] != a[1]:
            fraction = (1.0 - a[1]) / (b[1] - a[1])
            critical_z = a[0] + fraction * (b[0] - a[0])
            break

    peak = max(curve, key=lambda point: point[3])
    min_worth = min(point[3] for point in curve)
    worth_ymin = min(-25.0, math.floor(min_worth / 25.0) * 25.0)
    worth_ymax = max(100.0, rounded_ceiling(peak[3] * 1.08, 100.0))

    width, height = 2000, 1440
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    blue = "#1769aa"
    light_blue = "#e5f0f8"
    orange = "#d55e00"
    grid = "#d7dce2"
    text = "#17202a"
    gray = "#59636e"
    border = "#89939e"

    draw.text((105, 38), "CS reconstruction of the full control-rod worth curve",
              fill=text, font=get_font(44, True))

    legend_y = 102
    draw.line((108, legend_y, 164, legend_y), fill=blue, width=7)
    draw.text((178, legend_y - 17), "Hermite reconstruction", fill=text, font=get_font(23))
    draw.ellipse((508, legend_y - 9, 526, legend_y + 9), fill=blue)
    draw.text((540, legend_y - 17), "CS reference point", fill=text, font=get_font(23))
    draw.line((805, legend_y, 861, legend_y), fill=orange, width=5)
    draw.text((875, legend_y - 17), "Local CS slope", fill=text, font=get_font(23))
    draw.ellipse((1138, legend_y - 10, 1158, legend_y + 10), fill="white", outline=orange, width=4)
    draw.text((1172, legend_y - 17), "2-sigma validation flag", fill=text, font=get_font(23))

    zmin, zmax = rows[0]["z"], rows[-1]["z"]
    xticks = [float(value) for value in range(0, 128, 16)]
    if zmax > 128.0:
        xticks.append(zmax)

    def panel(box, ymin, ymax, ylabel, title, yticks):
        x0, y0, x1, y1 = box

        def sx(value):
            return x0 + (value - zmin) / (zmax - zmin) * (x1 - x0)

        def sy(value):
            return y1 - (value - ymin) / (ymax - ymin) * (y1 - y0)

        draw.text((x0, y0 - 48), title, fill=text, font=get_font(28, True))
        for value in yticks:
            yy = sy(value)
            draw.line((x0, yy, x1, yy), fill=grid, width=2)
            label = f"{value:.3f}" if ymax < 2.0 else f"{value:.0f}"
            bounds = draw.textbbox((0, 0), label, font=get_font(20))
            draw.text((x0 - 16 - (bounds[2] - bounds[0]), yy - 11), label,
                      fill=gray, font=get_font(20))
        for value in xticks:
            xx = sx(value)
            draw.line((xx, y0, xx, y1), fill=grid, width=2)
            label = f"{value:.2f}" if abs(value - round(value)) > 1e-6 else f"{int(round(value))}"
            bounds = draw.textbbox((0, 0), label, font=get_font(19))
            draw.text((xx - (bounds[2] - bounds[0]) / 2, y1 + 12), label,
                      fill=gray, font=get_font(19))
        draw.rectangle(box, outline=border, width=2)
        axis_label = "Rod insertion depth, z (cm)"
        bounds = draw.textbbox((0, 0), axis_label, font=get_font(23))
        draw.text(((x0 + x1 - (bounds[2] - bounds[0])) / 2, y1 + 48), axis_label,
                  fill=text, font=get_font(23))
        draw.text((24, (y0 + y1) / 2 - 14), ylabel, fill=text, font=get_font(23))
        return sx, sy

    top = (170, 200, 1915, 700)
    sx, sy = panel(top, 0.88, 1.16, "k_eff", "Reconstructed k_eff(z)",
                   [0.90, 0.95, 1.00, 1.05, 1.10, 1.15])
    draw.line([(sx(point[0]), sy(point[1])) for point in curve],
              fill=blue, width=6, joint="curve")

    for row in actual_rows:
        half = 2.2
        lo = max(zmin, row["z"] - half)
        hi = min(zmax, row["z"] + half)
        draw.line((sx(lo), sy(row["k"] + row["m"] * (lo - row["z"])),
                   sx(hi), sy(row["k"] + row["m"] * (hi - row["z"]))),
                  fill=orange, width=4)
        x, y = sx(row["z"]), sy(row["k"])
        draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill=blue, outline="white", width=3)
        if row["status"] != "PASS":
            draw.ellipse((x - 13, y - 13, x + 13, y + 13), outline=orange, width=4)

    # Final ARI endpoint, predicted from the 128.395 cm reference calculation.
    endpoint = rows[-1]
    if endpoint.get("endpoint_only"):
        x, y = sx(endpoint["z"]), sy(endpoint["k"])
        draw.rectangle((x - 7, y - 7, x + 7, y + 7), fill="white", outline=blue, width=3)

    draw.line((top[0], sy(1.0), top[2], sy(1.0)), fill=gray, width=3)
    if critical_z is not None:
        cx = sx(critical_z)
        draw.line((cx, top[1], cx, top[3]), fill="#7d8792", width=3)
        draw.text((cx + 14, sy(1.0) - 43), f"k = 1 near {critical_z:.3f} cm",
                  fill=text, font=get_font(22))

    bottom = (170, 835, 1915, 1260)
    tick_step = 200.0 if worth_ymax >= 800 else 100.0
    worth_ticks = []
    tick = 0.0
    while tick <= worth_ymax + 1e-9:
        worth_ticks.append(tick)
        tick += tick_step
    sx2, sy2 = panel(bottom, worth_ymin, worth_ymax, "-d rho/dz (pcm/cm)",
                     "Differential rod worth implied by point values and local CS slopes",
                     worth_ticks)
    polygon = [(sx2(curve[0][0]), sy2(0.0))]
    polygon.extend((sx2(point[0]), sy2(point[3])) for point in curve)
    polygon.append((sx2(curve[-1][0]), sy2(0.0)))
    draw.polygon(polygon, fill=light_blue)
    draw.line([(sx2(point[0]), sy2(point[3])) for point in curve],
              fill=blue, width=6, joint="curve")
    draw.line((bottom[0], sy2(0.0), bottom[2], sy2(0.0)), fill=gray, width=2)

    for row in actual_rows:
        worth = -row["m"] / (row["k"] * row["k"]) * 1.0e5
        x, y = sx2(row["z"]), sy2(worth)
        draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill=blue, outline="white", width=3)
        if row["status"] != "PASS":
            draw.ellipse((x - 13, y - 13, x + 13, y + 13), outline=orange, width=4)

    px, py = sx2(peak[0]), sy2(peak[3])
    draw.ellipse((px - 9, py - 9, px + 9, py + 9), fill=orange)
    label_x = min(px + 25, bottom[2] - 430)
    label_y = max(bottom[1] + 22, py - 52)
    draw.text((label_x, label_y),
              f"peak about {peak[3]:.0f} pcm/cm at {peak[0]:.1f} cm",
              fill=text, font=get_font(22))

    note = ("Cubic Hermite interpolation between 19 calculated reference depths. "
            "Orange rings mark points that failed at least one enabled 2-sigma validation check.")
    draw.text((170, 1370), note, fill=gray, font=get_font(20))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
