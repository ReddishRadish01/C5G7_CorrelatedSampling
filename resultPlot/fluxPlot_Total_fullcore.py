import os
import re
import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Tuple, Optional, List

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.colors import PowerNorm, Normalize


@dataclass
class AssemblyBlock:
    a: int
    start: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    length: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    dims: Tuple[int, int, int] = (0, 0, 0)
    data: Dict[Tuple[int, int, int], Tuple[float, float]] = field(default_factory=dict)


@dataclass
class CycleBlock:
    cycle_no: Optional[int] = None
    keff: Optional[float] = None
    core_size: Optional[Tuple[float, float, float]] = None
    z_slice_position: Optional[float] = None
    assembly_no: Optional[int] = None
    assemblies: Dict[int, AssemblyBlock] = field(default_factory=dict)


_core_re = re.compile(r"^\s*CoreSize\(cm\)\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+)\s*$")
_z_slice_re = re.compile(r"^\s*ZSlicePosition\(cm\)\s+([0-9eE\.\+\-]+)\s*$")
_asmno_re = re.compile(r"^\s*AssemblyNo\s+(\d+)\s*$")
_cycle_re = re.compile(r"^\s*(\d+)\s+th\s+cycle:\s*$")
_keff_re = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s+multiplication\s+Factor\s*$")
_asm_re = re.compile(r"^\s*Assembly\s+(\d+)\s*$")
_start_re = re.compile(r"^\s*StartPos\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+)\s*$")
_len_re = re.compile(r"^\s*Length\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+)\s*$")
_dims_re = re.compile(r"^\s*Dims\s+(\d+)\s+(\d+)\s+(\d+)\s*$")
_format_re = re.compile(r"^\s*Format\s+a\s+k\s+j\s+i\s+pinTally\s+modTally\s*$")


def parse_single_cycle_file(path: str) -> CycleBlock:
    cyc = CycleBlock()
    cur_asm: Optional[AssemblyBlock] = None
    in_data = False

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            m = _core_re.match(line)
            if m:
                cyc.core_size = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
                continue

            m = _z_slice_re.match(line)
            if m:
                cyc.z_slice_position = float(m.group(1))
                continue

            m = _asmno_re.match(line)
            if m:
                cyc.assembly_no = int(m.group(1))
                continue

            m = _cycle_re.match(line)
            if m:
                cyc.cycle_no = int(m.group(1))
                continue

            m = _keff_re.match(line)
            if m:
                cyc.keff = float(m.group(1))
                continue

            m = _asm_re.match(line)
            if m:
                a = int(m.group(1))
                cur_asm = cyc.assemblies.get(a)
                if cur_asm is None:
                    cur_asm = AssemblyBlock(a=a)
                    cyc.assemblies[a] = cur_asm
                in_data = False
                continue

            if cur_asm is None:
                continue

            m = _start_re.match(line)
            if m:
                cur_asm.start = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
                continue

            m = _len_re.match(line)
            if m:
                cur_asm.length = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
                continue

            m = _dims_re.match(line)
            if m:
                cur_asm.dims = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                continue

            if _format_re.match(line):
                in_data = True
                continue

            if in_data:
                parts = line.split()
                if len(parts) != 6:
                    continue
                try:
                    a = int(parts[0])
                    k = int(parts[1])
                    j = int(parts[2])
                    i = int(parts[3])
                    pin = float(parts[4])
                    mod = float(parts[5])
                except ValueError:
                    continue

                asm = cyc.assemblies.get(a)
                if asm is None:
                    asm = AssemblyBlock(a=a)
                    cyc.assemblies[a] = asm
                asm.data[(k, j, i)] = (pin, mod)

    return cyc


def _z_cell_bounds(asm: AssemblyBlock, z_cell: int) -> Tuple[float, float]:
    nz = asm.dims[2]
    if nz <= 0 or asm.length[2] <= 0.0:
        raise ValueError(f"Assembly {asm.a} has invalid axial dimensions")
    dz = asm.length[2] / nz
    z_low = asm.start[2] + z_cell * dz
    return z_low, z_low + dz


def _contains_z_position(asm: AssemblyBlock, z_cell: int, z_position_cm: float) -> bool:
    z_low, z_high = _z_cell_bounds(asm, z_cell)
    tolerance = 1.0e-10
    return z_position_cm >= z_low - tolerance and z_position_cm < z_high - tolerance


def _records_for_z_position(cycle: CycleBlock, z_position_cm: float):
    records = []
    for asm in cycle.assemblies.values():
        for (z_cell, j, i), (pin, mod) in asm.data.items():
            if _contains_z_position(asm, z_cell, z_position_cm):
                records.append((asm, z_cell, j, i, pin, mod))
    return records


def _collect_values_for_z_position(records, tally_mode: str) -> np.ndarray:
    vals = []
    for _asm, _z_cell, _j, _i, pin, mod in records:
        if tally_mode in ("total", "pin"):
            vals.append(pin)
        if tally_mode in ("total", "mod"):
            vals.append(mod)
    return np.asarray(vals, dtype=float)


_C5G7_NON_FUEL_PIN_INDICES = frozenset({
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5), (8, 8), (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
})


def _c5g7_fuel_assembly_ids(records):
    fuel_score_by_assembly = {}
    for asm, _z_cell, j, i, pin, _mod in records:
        if asm.dims[:2] != (17, 17) or (j, i) in _C5G7_NON_FUEL_PIN_INDICES:
            continue
        if np.isfinite(pin) and pin > 0.0:
            fuel_score_by_assembly[asm.a] = fuel_score_by_assembly.get(asm.a, 0.0) + pin
    return {
        assembly_id
        for assembly_id, score in fuel_score_by_assembly.items()
        if score > 0.0
    }


def _is_c5g7_fuel_pin(asm, j: int, i: int, fuel_assembly_ids) -> bool:
    return (
        asm.a in fuel_assembly_ids
        and asm.dims[:2] == (17, 17)
        and (j, i) not in _C5G7_NON_FUEL_PIN_INDICES
    )


def _relative_records_to_pin_mean(records, raw_pin_mean: float):
    fuel_assembly_ids = _c5g7_fuel_assembly_ids(records)
    return [
        (
            asm,
            z_cell,
            j,
            i,
            pin / raw_pin_mean
            if _is_c5g7_fuel_pin(asm, j, i, fuel_assembly_ids)
            else 0.0,
            mod / raw_pin_mean,
        )
        for asm, z_cell, j, i, pin, mod in records
    ]


def _pin_power_statistics(records):
    fuel_assembly_ids = _c5g7_fuel_assembly_ids(records)
    samples = []
    for asm, z_cell, j, i, pin, _mod in records:
        if not _is_c5g7_fuel_pin(asm, j, i, fuel_assembly_ids):
            continue
        if not np.isfinite(pin) or pin < 0.0:
            continue
        nx, ny, _nz = asm.dims
        dx = asm.length[0] / nx
        dy = asm.length[1] / ny
        x_center = asm.start[0] + (i + 0.5) * dx
        y_center = asm.start[1] + (j + 0.5) * dy
        samples.append((pin, asm.a, z_cell, j, i, x_center, y_center))

    if not samples:
        return None

    raw_values = np.asarray([sample[0] for sample in samples], dtype=float)
    raw_mean = float(np.mean(raw_values))
    min_sample = min(samples, key=lambda sample: sample[0])
    max_sample = max(samples, key=lambda sample: sample[0])
    return {
        "fuel_pin_count": len(samples),
        "raw_mean": raw_mean,
        "raw_min": float(min_sample[0]),
        "raw_max": float(max_sample[0]),
        "normalized_min": float(min_sample[0] / raw_mean),
        "normalized_max": float(max_sample[0] / raw_mean),
        "min_location": min_sample[1:],
        "max_location": max_sample[1:],
    }


def _infer_quarter_size(cycle: CycleBlock) -> Tuple[float, float]:
    if cycle.core_size is not None:
        return cycle.core_size[0], cycle.core_size[1]

    xmax = 0.0
    ymax = 0.0
    for asm in cycle.assemblies.values():
        sx, sy, _ = asm.start
        lx, ly, _ = asm.length
        xmax = max(xmax, sx + lx)
        ymax = max(ymax, sy + ly)
    return xmax, ymax


def plot_cycle_z_position(
    cycle: CycleBlock,
    z_position_cm: float,
    pin_radius_cm: float = 0.54,
    outpath: Optional[str] = None,
    show: bool = False,
    full_core_xy: bool = True,
    cmap_name: str = "turbo",
    gamma: float = 1.5,
    tally_mode: str = "total",
):
    records = _records_for_z_position(cycle, z_position_cm)
    statistics = _pin_power_statistics(records)
    if statistics is None or not np.isfinite(statistics["raw_mean"]) or statistics["raw_mean"] <= 0.0:
        raise RuntimeError(
            f"No positive pin-wise fission-rate tallies were found at "
            f"z position {z_position_cm:.6f} cm"
        )

    raw_pin_mean = statistics["raw_mean"]
    relative_records = _relative_records_to_pin_mean(records, raw_pin_mean)
    vals = _collect_values_for_z_position(relative_records, tally_mode)
    if vals.size == 0:
        raise RuntimeError(
            f"No tally cells containing z position {z_position_cm:.6f} cm were found.\n"
            "For a newly generated file, run c5g7 with --tally-z set to this position. "
            "Older files contain only the z positions that were dumped during that run."
        )

    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        raise RuntimeError(f"All values are NaN/Inf at z position {z_position_cm:.6f} cm")

    # Keep zero visible, but avoid one huge outlier destroying contrast.
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if vmax <= vmin:
        vmax = vmin + 1.0

    cmap = plt.get_cmap(cmap_name)
    if abs(gamma - 1.0) < 1.0e-12:
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)

    Lx, Ly = _infer_quarter_size(cycle)
    if Lx <= 0.0 or Ly <= 0.0:
        raise RuntimeError("Could not infer core x/y size from tally file.")

    fig_size = (10.5, 9.5) if full_core_xy else (9.0, 9.0)
    fig, ax = plt.subplots(figsize=fig_size)

    def add_rect_and_pin(x0, y0, dx, dy, cx, cy, pin, mod):
        if full_core_xy:
            # Plot full core shifted to [0, 2Lx] x [0, 2Ly].
            # Original quarter-core [0,Lx]x[0,Ly] is drawn in the upper-right quadrant.
            rects = [
                (Lx + x0,      Ly + y0),       # +x, +y
                (Lx - x0 - dx, Ly + y0),       # -x, +y mirrored
                (Lx + x0,      Ly - y0 - dy),  # +x, -y mirrored
                (Lx - x0 - dx, Ly - y0 - dy),  # -x, -y mirrored
            ]
            centers = [
                (Lx + cx, Ly + cy),
                (Lx - cx, Ly + cy),
                (Lx + cx, Ly - cy),
                (Lx - cx, Ly - cy),
            ]
        else:
            rects = [(x0, y0)]
            centers = [(cx, cy)]

        if tally_mode in ("total", "mod") and np.isfinite(mod):
            for rx0, ry0 in rects:
                ax.add_patch(Rectangle(
                    (rx0, ry0),
                    dx,
                    dy,
                    linewidth=0.02,
                    edgecolor="k",
                    facecolor=cmap(norm(mod)),
                ))

        if tally_mode in ("total", "pin") and pin != 0.0 and np.isfinite(pin):
            for pcx, pcy in centers:
                ax.add_patch(Circle(
                    (pcx, pcy),
                    radius=pin_radius_cm,
                    linewidth=0.1,
                    alpha=0.75,
                    edgecolor="k",
                    facecolor=cmap(norm(pin)),
                ))

    for asm, _z_cell, j, i, pin, mod in relative_records:
        nx, ny, _nz = asm.dims
        sx, sy, _sz = asm.start
        lx, ly, _lz = asm.length
        dx = lx / nx
        dy = ly / ny

        x0 = sx + i * dx
        y0 = sy + j * dy
        cx = x0 + 0.5 * dx
        cy = y0 + 0.5 * dy
        add_rect_and_pin(x0, y0, dx, dy, cx, cy, pin, mod)

    ax.set_aspect("equal", adjustable="box")

    selected_bounds = sorted({
        tuple(round(value, 12) for value in _z_cell_bounds(asm, z_cell))
        for asm, z_cell, _j, _i, _pin, _mod in records
    })
    bounds_text = ", ".join(f"[{low:g}, {high:g})" for low, high in selected_bounds)
    title = (
        f"Cycle {cycle.cycle_no if cycle.cycle_no is not None else ''} | "
        f"z position slice {z_position_cm:.3f} cm ({bounds_text})"
    )
    if cycle.keff is not None:
        title += f" | keff={cycle.keff:.6f}"
    if full_core_xy:
        title += " | mirrored full core"
    title += f" | raw pin mean={raw_pin_mean:.6g}"
    ax.set_title(title)

    if full_core_xy:
        ax.set_xlim(0.0, 2.0 * Lx)
        ax.set_ylim(0.0, 2.0 * Ly)
        ax.axvline(Lx, color="k", linewidth=0.8)
        ax.axhline(Ly, color="k", linewidth=0.8)
        ax.set_xlabel("x (cm), mirrored full core")
        ax.set_ylabel("y (cm), mirrored full core")
    else:
        ax.set_xlim(0.0, Lx)
        ax.set_ylim(0.0, Ly)
        ax.set_xlabel("x (cm), quarter core")
        ax.set_ylabel("y (cm), quarter core")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(
        f"Relative Fission Rate ({tally_mode}, z position slice pin mean = 1, "
        f"gamma={gamma:g})"
    )

    plt.tight_layout()

    if outpath is not None:
        outdir = os.path.dirname(outpath)
        if outdir:
            os.makedirs(outdir, exist_ok=True)
        plt.savefig(outpath, dpi=1000)
        print(f"saved: {outpath}")

    if show:
        plt.show()

    plt.close(fig)
    return statistics


def resolve_input_path(value: str) -> Path:
    requested = Path(value).expanduser()
    script_dir = Path(__file__).resolve().parent
    candidates = [
        requested,
        script_dir / requested,
        script_dir.parent / "results" / requested,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    checked = "\n  ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(f"Input tally file not found. Checked:\n  {checked}")


def _z_filename_token(z_position_cm: float) -> str:
    token = f"{z_position_cm:.3f}".replace("-", "m").replace(".", "p")
    return f"z{token}cm"


def _input_timestamp_token(path: Path) -> str:
    match = re.search(r"(\d{8}_\d{6})(?=\.[^.]+$|$)", path.name)
    if match:
        return match.group(1)
    return re.sub(r"[^A-Za-z0-9_-]+", "_", path.stem)


def write_pin_power_statistics(path: Path, z_position_cm: float, statistics) -> None:
    if statistics is None:
        text = (
            f"ZPosition(cm) {z_position_cm:.17g}\n"
            "Scope z_position_plane_only_not_axially_integrated_core\n"
            "BenchmarkComparable no\n"
            "Status no_positive_fuel_pin_tallies\n"
        )
    else:
        min_asm, min_z_cell, min_j, min_i, min_x, min_y = statistics["min_location"]
        max_asm, max_z_cell, max_j, max_i, max_x, max_y = statistics["max_location"]
        text = (
            f"ZPosition(cm) {z_position_cm:.17g}\n"
            "Scope z_position_plane_only_not_axially_integrated_core\n"
            "BenchmarkComparable no\n"
            "FuelSelection c5g7_17x17_geometry_mask\n"
            "Normalization z_position_plane_mean_fuel_pin_power_equals_1\n"
            f"FuelPinCount {statistics['fuel_pin_count']}\n"
            f"RawMean {statistics['raw_mean']:.17g}\n"
            f"MinimumRaw {statistics['raw_min']:.17g}\n"
            f"ZPositionMinimumNormalized {statistics['normalized_min']:.17g}\n"
            f"MinimumLocation assembly={min_asm} axial_cell={min_z_cell} "
            f"j={min_j} i={min_i} x_cm={min_x:.17g} y_cm={min_y:.17g}\n"
            f"MaximumRaw {statistics['raw_max']:.17g}\n"
            f"ZPositionMaximumNormalized {statistics['normalized_max']:.17g}\n"
            f"MaximumLocation assembly={max_asm} axial_cell={max_z_cell} "
            f"j={max_j} i={max_i} x_cm={max_x:.17g} y_cm={max_y:.17g}\n"
        )
    path.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Plot C5G7 fission-rate tally, optionally mirrored from quarter-core to full-core.")
    parser.add_argument("--infile", default="flux_Tally_reference2003_20260727_201248.txt", help="Input tally text file")
    parser.add_argument("--out-dir", default="plots_Total", help="Output directory")
    parser.add_argument(
        "--z-position",
        type=float,
        default=None,
        help="Global z position in cm. Defaults to ZSlicePosition(cm) in the file, or 0.63 cm for an older file.",
    )
    parser.add_argument("--quarter", action="store_true", help="Plot quarter core only. Default is mirrored full core.")
    parser.add_argument("--cmap", default="turbo", help="Matplotlib colormap name, e.g. turbo, plasma, inferno")
    parser.add_argument("--gamma", type=float, default=1.5, help="PowerNorm gamma. Use 1.0 for linear.")
    parser.add_argument("--pin-radius", type=float, default=0.54, help="Pin radius in cm")
    parser.add_argument("--mode", choices=["total", "pin", "mod"], default="total", help="Tally values used for plotting/color scaling")
    parser.add_argument("--show", action="store_true", help="Show matplotlib window")
    args = parser.parse_args()

    infile = resolve_input_path(args.infile)
    cycle = parse_single_cycle_file(str(infile))
    z_position_cm = (
        args.z_position
        if args.z_position is not None
        else cycle.z_slice_position if cycle.z_slice_position is not None
        else 0.63
    )
    full_core = not args.quarter

    suffix = "fullcore" if full_core else "quarter"
    z_token = _z_filename_token(z_position_cm)
    input_timestamp = _input_timestamp_token(infile)
    outname = (
        f"fission_rate_{suffix}_cycle{cycle.cycle_no}_{z_token}_"
        f"input{input_timestamp}_{args.cmap}_g{args.gamma:g}.png"
    )
    outpath = Path(args.out_dir) / outname

    statistics = plot_cycle_z_position(
        cycle,
        z_position_cm=z_position_cm,
        pin_radius_cm=args.pin_radius,
        outpath=str(outpath),
        show=args.show,
        full_core_xy=full_core,
        cmap_name=args.cmap,
        gamma=args.gamma,
        tally_mode=args.mode,
    )
    statistics_path = outpath.with_suffix(".pin_power.txt")
    write_pin_power_statistics(statistics_path, z_position_cm, statistics)
    print(f"loaded: {infile}")
    print(f"pin power statistics: {statistics_path.resolve()}")


if __name__ == "__main__":
    main()
