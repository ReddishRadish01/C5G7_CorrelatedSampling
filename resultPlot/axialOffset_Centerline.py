import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


_C5G7_NON_FUEL_PIN_INDICES = frozenset({
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5), (8, 8), (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
})

_AXIAL_POWER_X_LABEL = (
    "Relative centerline axial power\n"
    "(height-weighted mean = 1)"
)


@dataclass
class AssemblyMeta:
    assembly_id: int
    start: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    length: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    dims: Tuple[int, int, int] = (0, 0, 0)


@dataclass
class SelectedScore:
    assembly_id: int
    z_low: float
    z_high: float
    pin_tally: float


@dataclass
class AxialBin:
    z_low: float
    z_high: float
    raw_sum: float = 0.0
    pin_count: int = 0
    mean_pin_score: float = 0.0
    linear_density: float = 0.0
    relative_power: float = 0.0

    @property
    def z_center(self) -> float:
        return 0.5 * (self.z_low + self.z_high)

    @property
    def height(self) -> float:
        return self.z_high - self.z_low


@dataclass
class AxialOffsetResult:
    bins: List[AxialBin]
    center_x: float
    center_y: float
    radius_cm: float
    midplane_z: float
    bottom_power: float
    top_power: float
    axial_offset: float
    fuel_assembly_count: int
    tally_type: str
    modeled_extent: str


def _mirrored_bin(axial_bin: AxialBin) -> AxialBin:
    return AxialBin(
        z_low=-axial_bin.z_high,
        z_high=-axial_bin.z_low,
        raw_sum=axial_bin.raw_sum,
        pin_count=axial_bin.pin_count,
        mean_pin_score=axial_bin.mean_pin_score,
        linear_density=axial_bin.linear_density,
        relative_power=axial_bin.relative_power,
    )


def _display_bins(result: AxialOffsetResult) -> List[AxialBin]:
    bins = list(result.bins)
    if result.modeled_extent == "top_half":
        bins.extend(_mirrored_bin(axial_bin) for axial_bin in result.bins)
    elif result.modeled_extent == "bottom_half":
        bins.extend(_mirrored_bin(axial_bin) for axial_bin in result.bins)
    return sorted(bins, key=lambda item: item.z_center)


def _half_labels(result: AxialOffsetResult) -> Tuple[str, str]:
    if result.modeled_extent == "top_half":
        return "Top half (modeled)", "Bottom half - reflective"
    if result.modeled_extent == "bottom_half":
        return "Top half - reflective", "Bottom half (modeled)"
    return "Top half", "Bottom half"


def _resolve_input_path(value: Optional[str]) -> Path:
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    if value:
        requested = Path(value).expanduser()
        candidates = [requested] if requested.is_absolute() else [
            script_dir / requested,
            requested,
            project_dir / "results" / requested,
            project_dir / requested,
        ]
        for candidate in candidates:
            if candidate.is_file():
                return candidate.resolve()
        checked = "\n  ".join(str(candidate) for candidate in candidates)
        raise FileNotFoundError(f"Input pinwise tally was not found. Checked:\n  {checked}")

    candidates = list((project_dir / "results").glob("pinwise_Tally_*.txt"))
    if not candidates:
        raise FileNotFoundError(f"No pinwise_Tally_*.txt files found in {project_dir / 'results'}")
    return max(candidates, key=lambda path: path.stat().st_mtime).resolve()


def _input_timestamp_token(path: Path) -> str:
    match = re.search(r"(\d{8}_\d{6})(?=\.[^.]+$|$)", path.name)
    if match:
        return match.group(1)
    return re.sub(r"[^A-Za-z0-9_-]+", "_", path.stem)


def _number_token(value: float) -> str:
    return f"{value:.3f}".replace("-", "m").replace(".", "p")


def _read_core_metadata(path: Path) -> Tuple[str, Tuple[float, float, float, float, float, float]]:
    tally_type = ""
    core_bounds = None
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith("TallyType "):
                tally_type = line.split(maxsplit=1)[1]
            elif line.startswith("CoreBounds(cm) "):
                parts = line.split()
                if len(parts) == 7:
                    core_bounds = tuple(float(value) for value in parts[1:])
            elif line.startswith("Assembly "):
                break
    if not tally_type:
        raise ValueError("TallyType is missing from the pinwise tally file")
    if core_bounds is None:
        raise ValueError("CoreBounds(cm) is missing from the pinwise tally file")
    return tally_type, core_bounds


def _select_centerline_scores(
    path: Path,
    center_x: float,
    center_y: float,
    radius_cm: float,
) -> Tuple[List[SelectedScore], int]:
    current: Optional[AssemblyMeta] = None
    fuel_score_by_assembly: Dict[int, float] = {}
    selected: List[SelectedScore] = []
    data_row = re.compile(
        r"^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
        r"([0-9eE.+-]+)\s+([0-9eE.+-]+)$"
    )

    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("Assembly "):
                current = AssemblyMeta(int(line.split()[1]))
                continue
            if current is None:
                continue
            if line.startswith("StartPos(cm) "):
                parts = line.split()
                current.start = tuple(float(value) for value in parts[1:4])
                continue
            if line.startswith("Length(cm) "):
                parts = line.split()
                current.length = tuple(float(value) for value in parts[1:4])
                continue
            if line.startswith("Dims "):
                parts = line.split()
                current.dims = tuple(int(value) for value in parts[1:4])
                continue

            match = data_row.match(line)
            if not match:
                continue
            assembly_id, z_index, j, i = (int(match.group(index)) for index in range(1, 5))
            pin_tally = float(match.group(5))
            if assembly_id != current.assembly_id or current.dims[:2] != (17, 17):
                continue
            if (j, i) in _C5G7_NON_FUEL_PIN_INDICES:
                continue
            if not math.isfinite(pin_tally) or pin_tally < 0.0:
                continue

            if pin_tally > 0.0:
                fuel_score_by_assembly[assembly_id] = (
                    fuel_score_by_assembly.get(assembly_id, 0.0) + pin_tally
                )

            nx, ny, nz = current.dims
            if nx <= 0 or ny <= 0 or nz <= 0:
                continue
            dx = current.length[0] / nx
            dy = current.length[1] / ny
            dz = current.length[2] / nz
            x_center = current.start[0] + (i + 0.5) * dx
            y_center = current.start[1] + (j + 0.5) * dy
            if math.hypot(x_center - center_x, y_center - center_y) > radius_cm + 1.0e-12:
                continue

            z_low = current.start[2] + z_index * dz
            selected.append(SelectedScore(assembly_id, z_low, z_low + dz, pin_tally))

    fuel_assembly_ids = {
        assembly_id
        for assembly_id, score in fuel_score_by_assembly.items()
        if score > 0.0
    }
    selected = [score for score in selected if score.assembly_id in fuel_assembly_ids]
    return selected, len(fuel_assembly_ids)


def calculate_axial_offset(
    path: Path,
    center_x: Optional[float] = None,
    center_y: Optional[float] = None,
    radius_cm: float = 2.0,
) -> AxialOffsetResult:
    tally_type, bounds = _read_core_metadata(path)
    if tally_type != "PinwiseFissionRateCollisionEstimator":
        raise ValueError(
            "Axial power requires TallyType PinwiseFissionRateCollisionEstimator; "
            f"found {tally_type}"
        )
    _xmin, _xmax, _ymin, _ymax, _zmin, _zmax = bounds
    center_x = 0.0 if center_x is None else center_x
    center_y = 0.0 if center_y is None else center_y
    if radius_cm <= 0.0:
        raise ValueError("radius_cm must be positive")

    selected, fuel_assembly_count = _select_centerline_scores(
        path, center_x, center_y, radius_cm
    )
    if not selected:
        raise RuntimeError(
            f"No fuel pins were selected within {radius_cm:g} cm of "
            f"({center_x:g}, {center_y:g})"
        )

    grouped: Dict[Tuple[int, int], AxialBin] = {}
    for score in selected:
        key = (
            int(round(score.z_low * 1.0e9)),
            int(round(score.z_high * 1.0e9)),
        )
        axial_bin = grouped.get(key)
        if axial_bin is None:
            axial_bin = AxialBin(score.z_low, score.z_high)
            grouped[key] = axial_bin
        axial_bin.raw_sum += score.pin_tally
        axial_bin.pin_count += 1

    bins = sorted(grouped.values(), key=lambda item: item.z_center)
    for axial_bin in bins:
        if axial_bin.pin_count <= 0 or axial_bin.height <= 0.0:
            raise RuntimeError("Invalid axial bin encountered")
        axial_bin.mean_pin_score = axial_bin.raw_sum / axial_bin.pin_count
        axial_bin.linear_density = axial_bin.mean_pin_score / axial_bin.height

    total_height = sum(axial_bin.height for axial_bin in bins)
    mean_density = (
        sum(axial_bin.linear_density * axial_bin.height for axial_bin in bins)
        / total_height
    )
    if mean_density <= 0.0:
        raise RuntimeError("Centerline axial fission-rate sum is zero")
    for axial_bin in bins:
        axial_bin.relative_power = axial_bin.linear_density / mean_density

    midplane_z = 0.0
    axial_tolerance = max(1.0e-9, min(axial_bin.height for axial_bin in bins) * 1.0e-6)
    fuel_z_min = min(axial_bin.z_low for axial_bin in bins)
    fuel_z_max = max(axial_bin.z_high for axial_bin in bins)
    if fuel_z_min >= midplane_z - axial_tolerance:
        modeled_extent = "top_half"
    elif fuel_z_max <= midplane_z + axial_tolerance:
        modeled_extent = "bottom_half"
    else:
        modeled_extent = "full_core"

    if modeled_extent == "top_half":
        top_power = sum(axial_bin.mean_pin_score for axial_bin in bins)
        bottom_power = top_power
    elif modeled_extent == "bottom_half":
        bottom_power = sum(axial_bin.mean_pin_score for axial_bin in bins)
        top_power = bottom_power
    else:
        bottom_power = 0.0
        top_power = 0.0
        for axial_bin in bins:
            height = axial_bin.height
            bottom_fraction = max(
                0.0,
                min(axial_bin.z_high, midplane_z) - axial_bin.z_low,
            ) / height
            top_fraction = max(
                0.0,
                axial_bin.z_high - max(axial_bin.z_low, midplane_z),
            ) / height
            bottom_power += axial_bin.mean_pin_score * bottom_fraction
            top_power += axial_bin.mean_pin_score * top_fraction

    total_power = bottom_power + top_power
    if total_power <= 0.0:
        raise RuntimeError("Top and bottom centerline powers sum to zero")
    axial_offset = (top_power - bottom_power) / total_power

    return AxialOffsetResult(
        bins=bins,
        center_x=center_x,
        center_y=center_y,
        radius_cm=radius_cm,
        midplane_z=midplane_z,
        bottom_power=bottom_power,
        top_power=top_power,
        axial_offset=axial_offset,
        fuel_assembly_count=fuel_assembly_count,
        tally_type=tally_type,
        modeled_extent=modeled_extent,
    )


def _write_csv(path: Path, result: AxialOffsetResult) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "z_center_cm",
            "z_low_cm",
            "z_high_cm",
            "selected_fuel_pin_count",
            "raw_fission_rate_sum",
            "mean_pin_fission_rate",
            "linear_fission_rate_density",
            "relative_axial_power",
        ])
        for axial_bin in result.bins:
            writer.writerow([
                f"{axial_bin.z_center:.17g}",
                f"{axial_bin.z_low:.17g}",
                f"{axial_bin.z_high:.17g}",
                axial_bin.pin_count,
                f"{axial_bin.raw_sum:.17g}",
                f"{axial_bin.mean_pin_score:.17g}",
                f"{axial_bin.linear_density:.17g}",
                f"{axial_bin.relative_power:.17g}",
            ])


def _write_summary(path: Path, infile: Path, result: AxialOffsetResult) -> None:
    maximum = max(result.bins, key=lambda item: item.relative_power)
    minimum = min(result.bins, key=lambda item: item.relative_power)
    text = (
        "CenterlineAxialOffsetSummary\n"
        f"InputFile {infile}\n"
        f"TallyType {result.tally_type}\n"
        "FuelSelection c5g7_17x17_geometry_mask\n"
        f"CenterXY(cm) {result.center_x:.17g} {result.center_y:.17g}\n"
        f"SelectionRadius(cm) {result.radius_cm:.17g}\n"
        f"FuelAssemblyCount {result.fuel_assembly_count}\n"
        f"AxialBinCount {len(result.bins)}\n"
        f"ModeledAxialExtent {result.modeled_extent}\n"
        f"ReflectiveReconstruction {'yes' if result.modeled_extent != 'full_core' else 'no'}\n"
        "Normalization axial_height_weighted_mean_equals_1\n"
        "AxialOffsetDefinition (top-bottom)/(top+bottom)\n"
        f"MidplaneZ(cm) {result.midplane_z:.17g}\n"
        f"BottomPower {result.bottom_power:.17g}\n"
        f"TopPower {result.top_power:.17g}\n"
        f"AxialOffset {result.axial_offset:.17g}\n"
        f"AxialOffsetPercent {100.0 * result.axial_offset:.17g}\n"
        f"MinimumRelativePower {minimum.relative_power:.17g} ZCenter(cm) {minimum.z_center:.17g}\n"
        f"MaximumRelativePower {maximum.relative_power:.17g} ZCenter(cm) {maximum.z_center:.17g}\n"
    )
    path.write_text(text, encoding="utf-8")


def _plot_curve_points(
    display_bins: List[AxialBin],
    midplane_z: float,
) -> List[Tuple[float, float]]:
    """Return (z, power) points with a shared, interpolated midplane point."""
    points = [
        (axial_bin.z_center, axial_bin.relative_power)
        for axial_bin in display_bins
    ]
    if any(math.isclose(z, midplane_z, rel_tol=0.0, abs_tol=1.0e-12) for z, _ in points):
        return points

    for index in range(len(points) - 1):
        z0, power0 = points[index]
        z1, power1 = points[index + 1]
        if z0 < midplane_z < z1:
            fraction = (midplane_z - z0) / (z1 - z0)
            midplane_power = power0 + fraction * (power1 - power0)
            points.insert(index + 1, (midplane_z, midplane_power))
            break
    return points


def _render_matplotlib(
    path: Path,
    result: AxialOffsetResult,
    font_scale: float = 1.0,
    show: bool = False,
) -> None:
    import matplotlib.pyplot as plt

    display_bins = _display_bins(result)
    top_label, bottom_label = _half_labels(result)
    curve_points = _plot_curve_points(display_bins, result.midplane_z)
    z = [point[0] for point in curve_points]
    relative = [point[1] for point in curve_points]
    maximum = max(result.bins, key=lambda item: item.relative_power)

    fig, ax = plt.subplots(figsize=(7.5, 9.5))
    ax.plot(relative, z, color="black", linewidth=1.6)
    ax.fill_betweenx(
        z, 0.0, relative,
        where=[value >= result.midplane_z for value in z],
        color="tab:red", alpha=0.24, label=top_label,
    )
    ax.fill_betweenx(
        z, 0.0, relative,
        where=[value <= result.midplane_z for value in z],
        color="tab:blue", alpha=0.24, label=bottom_label,
    )
    ax.axhline(result.midplane_z, color="0.25", linestyle="--", linewidth=1.0)
    ax.axvline(1.0, color="0.55", linestyle=":", linewidth=1.0)
    ax.scatter([maximum.relative_power], [maximum.z_center], color="black", s=28, zorder=3)
    ax.annotate(
        f"Max {maximum.relative_power:.3f} at z={maximum.z_center:.2f} cm",
        (maximum.relative_power, maximum.z_center),
        xytext=(8, 8), textcoords="offset points",
        fontsize=10.0 * font_scale,
    )
    ax.set_xlabel(
        _AXIAL_POWER_X_LABEL,
        fontsize=12.0 * font_scale,
        labelpad=8.0 * font_scale,
    )
    ax.set_ylabel("z position (cm)", fontsize=12.0 * font_scale)
    ax.set_title(
        f"Centerline Axial Power | AO={100.0 * result.axial_offset:+.3f}%\n"
        f"axis=({result.center_x:.3f}, {result.center_y:.3f}, z) cm, "
        f"radius={result.radius_cm:.3f} cm",
        fontsize=14.0 * font_scale,
    )
    ax.tick_params(axis="both", labelsize=10.0 * font_scale)
    ax.grid(True, linewidth=0.5, alpha=0.35)
    ax.legend(loc="best", fontsize=10.0 * font_scale)
    fig.tight_layout(pad=1.2)
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.18)
    if show:
        plt.show()
    plt.close(fig)


def _render_pillow(
    path: Path,
    result: AxialOffsetResult,
    font_scale: float = 1.0,
) -> None:
    from PIL import Image, ImageDraw, ImageFont

    def scaled(value: float) -> int:
        return max(1, int(round(value * font_scale)))

    width, height = 1400, 1700
    left, right, top, bottom = (
        scaled(190), scaled(90), scaled(175), scaled(170)
    )
    plot_width = width - left - right
    plot_height = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image, "RGBA")

    font_path = Path("C:/Windows/Fonts/arial.ttf")
    bold_path = Path("C:/Windows/Fonts/arialbd.ttf")
    font = ImageFont.truetype(str(font_path), scaled(28)) if font_path.is_file() else ImageFont.load_default()
    small = ImageFont.truetype(str(font_path), scaled(23)) if font_path.is_file() else ImageFont.load_default()
    title_font = ImageFont.truetype(str(bold_path), scaled(34)) if bold_path.is_file() else font

    display_bins = _display_bins(result)
    top_label, bottom_label = _half_labels(result)
    zmin = min(axial_bin.z_low for axial_bin in display_bins)
    zmax = max(axial_bin.z_high for axial_bin in display_bins)
    xmax = max(axial_bin.relative_power for axial_bin in display_bins) * 1.12
    xmax = max(1.2, xmax)

    def px(value: float) -> float:
        return left + value / xmax * plot_width

    def py(value: float) -> float:
        return top + (zmax - value) / (zmax - zmin) * plot_height

    for tick in range(6):
        value = xmax * tick / 5.0
        x = px(value)
        draw.line((x, top, x, top + plot_height), fill=(190, 190, 190, 120), width=1)
        label = f"{value:.2f}"
        draw.text(
            (x - scaled(24), top + plot_height + scaled(20)),
            label, fill="black", font=small,
        )
    for tick in range(7):
        value = zmin + (zmax - zmin) * tick / 6.0
        y = py(value)
        draw.line((left, y, left + plot_width, y), fill=(190, 190, 190, 120), width=1)
        draw.text(
            (left - scaled(115), y - scaled(14)),
            f"{value:.1f}", fill="black", font=small,
        )

    curve_points = _plot_curve_points(display_bins, result.midplane_z)
    points = [(px(relative_power), py(z_center)) for z_center, relative_power in curve_points]
    baseline = px(0.0)
    bottom_points = [
        point
        for point, (z_center, _relative_power) in zip(points, curve_points)
        if z_center <= result.midplane_z
    ]
    top_points = [
        point
        for point, (z_center, _relative_power) in zip(points, curve_points)
        if z_center >= result.midplane_z
    ]
    if bottom_points:
        draw.polygon(
            [(baseline, bottom_points[0][1]), *bottom_points, (baseline, bottom_points[-1][1])],
            fill=(55, 105, 190, 70),
        )
    if top_points:
        draw.polygon(
            [(baseline, top_points[0][1]), *top_points, (baseline, top_points[-1][1])],
            fill=(210, 70, 60, 70),
        )
    draw.line(points, fill=(25, 25, 25, 255), width=4, joint="curve")

    mid_y = py(result.midplane_z)
    for x in range(left, left + plot_width, 22):
        draw.line((x, mid_y, min(x + 12, left + plot_width), mid_y), fill=(70, 70, 70, 255), width=2)
    mean_x = px(1.0)
    for y in range(top, top + plot_height, 22):
        draw.line((mean_x, y, mean_x, min(y + 12, top + plot_height)), fill=(105, 105, 105, 220), width=2)

    draw.rectangle((left, top, left + plot_width, top + plot_height), outline="black", width=3)
    maximum = max(result.bins, key=lambda item: item.relative_power)
    max_point = (px(maximum.relative_power), py(maximum.z_center))
    draw.ellipse(
        (max_point[0] - 7, max_point[1] - 7, max_point[0] + 7, max_point[1] + 7),
        fill="black",
    )
    maximum_label = f"Max {maximum.relative_power:.3f}, z={maximum.z_center:.2f} cm"
    label_box = draw.textbbox((0, 0), maximum_label, font=small)
    label_width = label_box[2] - label_box[0]
    label_x = min(
        max_point[0] + scaled(14),
        left + plot_width - label_width - scaled(12),
    )
    draw.text(
        (label_x, max_point[1] - scaled(36)),
        maximum_label, fill="black", font=small,
    )

    title = f"Centerline Axial Power   AO = {100.0 * result.axial_offset:+.3f}%"
    subtitle = (
        f"axis=({result.center_x:.3f}, {result.center_y:.3f}, z) cm, "
        f"radius={result.radius_cm:.3f} cm"
    )
    draw.text((left, scaled(45)), title, fill="black", font=title_font)
    draw.text((left, scaled(100)), subtitle, fill="black", font=font)
    legend_x = left + plot_width - scaled(340)
    legend_y = top + scaled(28)
    draw.rectangle(
        (legend_x, legend_y, legend_x + scaled(32), legend_y + scaled(22)),
        fill=(210, 70, 60, 90),
    )
    draw.text(
        (legend_x + scaled(44), legend_y - scaled(4)),
        top_label, fill="black", font=small,
    )
    draw.rectangle(
        (
            legend_x, legend_y + scaled(38),
            legend_x + scaled(32), legend_y + scaled(60),
        ),
        fill=(55, 105, 190, 90),
    )
    draw.text(
        (legend_x + scaled(44), legend_y + scaled(34)),
        bottom_label, fill="black", font=small,
    )
    x_label = _AXIAL_POWER_X_LABEL
    x_label_spacing = scaled(6)
    x_label_box = draw.multiline_textbbox(
        (0, 0),
        x_label,
        font=font,
        spacing=x_label_spacing,
        align="center",
    )
    x_label_width = x_label_box[2] - x_label_box[0]
    x_label_height = x_label_box[3] - x_label_box[1]
    draw.multiline_text(
        (
            left + (plot_width - x_label_width) / 2,
            height - scaled(40) - x_label_height - x_label_box[1],
        ),
        x_label,
        fill="black",
        font=font,
        spacing=x_label_spacing,
        align="center",
    )
    y_label_box = draw.textbbox((0, 0), "z position (cm)", font=font)
    y_label_width = y_label_box[2] - y_label_box[0]
    y_label_height = y_label_box[3] - y_label_box[1]
    y_label_image = Image.new(
        "RGBA",
        (y_label_width + scaled(12), y_label_height + scaled(12)),
        (255, 255, 255, 0),
    )
    y_label_draw = ImageDraw.Draw(y_label_image)
    y_label_draw.text(
        (scaled(6), scaled(6) - y_label_box[1]),
        "z position (cm)", fill="black", font=font,
    )
    y_label_image = y_label_image.rotate(90, expand=True)
    image.paste(
        y_label_image,
        (scaled(20), top + (plot_height - y_label_image.height) // 2),
        y_label_image,
    )
    image.save(path, dpi=(300, 300))


def _render_plot(
    path: Path,
    result: AxialOffsetResult,
    backend: str,
    font_scale: float = 1.0,
    show: bool = False,
) -> str:
    if backend in ("auto", "matplotlib"):
        try:
            _render_matplotlib(path, result, font_scale=font_scale, show=show)
            return "matplotlib"
        except ModuleNotFoundError:
            if backend == "matplotlib":
                raise
    _render_pillow(path, result, font_scale=font_scale)
    return "pillow"


def main(
    default_infile: Optional[str] = None,
    default_font_scale: float = 1.0,
) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate and plot the C5G7 centerline axial fission-rate distribution "
            "and axial offset from a full pinwise tally."
        )
    )
    parser.add_argument(
        "--infile",
        default=default_infile,
        help=(
            "Pinwise tally path. Defaults to PINWISE_INPUT_FILENAME at the bottom "
            "of this script, or the newest results/pinwise_Tally_*.txt when blank."
        ),
    )
    parser.add_argument("--out-dir", default="resultPlot/plots_AxialOffset")
    parser.add_argument(
        "--center-x", type=float, default=None,
        help="Centerline x position in cm. Defaults to the symmetry axis x=0.",
    )
    parser.add_argument(
        "--center-y", type=float, default=None,
        help="Centerline y position in cm. Defaults to the symmetry axis y=0.",
    )
    parser.add_argument(
        "--radius-cm",
        type=float,
        default=2.0,
        help="Radius around the centerline used to average nearby fuel pins",
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "matplotlib", "pillow"],
        default="auto",
        help="PNG renderer. Auto falls back to Pillow when matplotlib is unavailable.",
    )
    parser.add_argument(
        "--font-scale",
        type=float,
        default=default_font_scale,
        help="Scale all plot fonts together; for example 1.3 makes them 30 percent larger.",
    )
    parser.add_argument("--show", action="store_true", help="Show the plot when matplotlib is used")
    args = parser.parse_args()
    if args.font_scale <= 0.0:
        parser.error("--font-scale must be greater than zero")

    infile = _resolve_input_path(args.infile)
    result = calculate_axial_offset(
        infile,
        center_x=args.center_x,
        center_y=args.center_y,
        radius_cm=args.radius_cm,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    timestamp = _input_timestamp_token(infile)
    location = (
        f"x{_number_token(result.center_x)}_y{_number_token(result.center_y)}_"
        f"r{_number_token(result.radius_cm)}"
    )
    stem = f"axial_offset_centerline_{location}_input{timestamp}"
    png_path = (out_dir / f"{stem}.png").resolve()
    csv_path = (out_dir / f"{stem}.csv").resolve()
    summary_path = (out_dir / f"{stem}.summary.txt").resolve()

    _write_csv(csv_path, result)
    _write_summary(summary_path, infile, result)
    backend = _render_plot(
        png_path,
        result,
        args.backend,
        font_scale=args.font_scale,
        show=args.show,
    )

    print(f"Input: {infile}")
    print(f"Centerline selection: ({result.center_x:.6f}, {result.center_y:.6f}) cm")
    print(f"Selection radius: {result.radius_cm:.6f} cm")
    print(f"Axial bins: {len(result.bins)}")
    print(f"Modeled axial extent: {result.modeled_extent}")
    print(f"Axial Offset: {result.axial_offset:+.8f} ({100.0 * result.axial_offset:+.5f}%)")
    print(f"PNG ({backend}): {png_path}")
    print(f"CSV: {csv_path}")
    print(f"Summary: {summary_path}")

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # USER SETTING: put the pinwise tally beside this Python file and
    # change only the filename below.
    # ------------------------------------------------------------------
    PINWISE_INPUT_FILENAME = "pinwise_Tally_CriticalSearch_Binary_20260730_044025.txt"
    FONT_SCALE = 1.30

    main(
        default_infile=PINWISE_INPUT_FILENAME,
        default_font_scale=FONT_SCALE,
    )
