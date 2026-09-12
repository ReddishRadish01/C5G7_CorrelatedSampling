import argparse
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.colors import Normalize, PowerNorm
from fluxPlot_Total_fullcore import (
    _input_timestamp_token,
    _pin_power_statistics,
    _relative_records_to_pin_mean,
    _z_filename_token,
    write_pin_power_statistics,
)


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
    core_bounds: Optional[Tuple[float, float, float, float, float, float]] = None
    z_slice_position: Optional[float] = None
    assembly_no: Optional[int] = None
    assemblies: Dict[int, AssemblyBlock] = field(default_factory=dict)


_number = r"([0-9eE\.\+\-]+)"
_core_size_re = re.compile(
    rf"^\s*CoreSize\(cm\)\s+{_number}\s+{_number}\s+{_number}\s*$"
)
_core_bounds_re = re.compile(
    rf"^\s*CoreBounds\(cm\)\s+{_number}\s+{_number}\s+{_number}\s+"
    rf"{_number}\s+{_number}\s+{_number}\s*$"
)
_z_slice_re = re.compile(rf"^\s*ZSlicePosition\(cm\)\s+{_number}\s*$")
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

            m = _core_size_re.match(line)
            if m:
                cyc.core_size = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
                continue

            m = _core_bounds_re.match(line)
            if m:
                cyc.core_bounds = tuple(float(m.group(i)) for i in range(1, 7))
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


def minmax01(values: np.ndarray) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v)
    out = np.zeros_like(v, dtype=float)
    if not np.any(finite):
        return out
    vmin = np.min(v[finite])
    vmax = np.max(v[finite])
    if vmax <= vmin:
        out[finite] = 1.0
        return out
    out[finite] = (v[finite] - vmin) / (vmax - vmin)
    return out


def _z_cell_bounds(asm: AssemblyBlock, z_cell: int) -> Tuple[float, float]:
    nz = asm.dims[2]
    if nz <= 0 or asm.length[2] <= 0.0:
        raise ValueError(f"Assembly {asm.a} has invalid axial dimensions")
    dz = asm.length[2] / nz
    z_low = asm.start[2] + z_cell * dz
    return z_low, z_low + dz


def _records_for_z_position(cycle: CycleBlock, z_position_cm: float):
    records = []
    tolerance = 1.0e-10
    for asm in cycle.assemblies.values():
        for (z_cell, j, i), (pin, mod) in asm.data.items():
            z_low, z_high = _z_cell_bounds(asm, z_cell)
            if z_position_cm >= z_low - tolerance and z_position_cm < z_high - tolerance:
                records.append((asm, z_cell, j, i, pin, mod))
    return records


def plot_cycle_z_position(
    cycle: CycleBlock,
    z_position_cm: float,
    pin_radius_cm: float = 0.54,
    outpath: Optional[str] = None,
    show: bool = False,
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
    vals = []
    for _asm, _z_cell, _j, _i, pin, mod in relative_records:
        vals.append(pin)
        vals.append(mod)

    if len(vals) == 0:
        raise RuntimeError(
            f"No tally cells containing z position {z_position_cm:.6f} cm were found. "
            "Run c5g7 with --tally-z set to this position when creating a new tally file."
        )

    vals = np.asarray(vals, dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        raise RuntimeError(f"All values are NaN/Inf at z position {z_position_cm:.6f} cm")

    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if vmax <= vmin:
        vmax = vmin + 1.0
    #norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("turbo")
    norm = PowerNorm(gamma=1.6, vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=(9, 9))

    for asm, _z_cell, j, i, pin, mod in relative_records:
        nx, ny, _nz = asm.dims
        sx, sy, _ = asm.start
        lx, ly, _ = asm.length
        dx = lx / nx
        dy = ly / ny

        x0 = sx + i * dx
        y0 = sy + j * dy
        cx = x0 + 0.5 * dx
        cy = y0 + 0.5 * dy

        sq = Rectangle(
            (x0, y0),
            dx,
            dy,
            linewidth=0.02,
            edgecolor="k",
            facecolor=cmap(norm(mod)),
        )
        ax.add_patch(sq)

        if pin != 0.0 and np.isfinite(pin):
            circ = Circle(
                (cx, cy),
                radius=pin_radius_cm,
                linewidth=0.1,
                alpha=0.5,
                edgecolor="k",
                facecolor=cmap(norm(pin)),
            )
            ax.add_patch(circ)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")

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
    title += f" | raw pin mean={raw_pin_mean:.6g}"
    ax.set_title(title)

    #sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=Normalize(vmin=0.0, vmax=1.0))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Relative Fission Rate (z position slice pin mean = 1)")

    if cycle.core_bounds is not None:
        xmin, xmax, ymin, ymax, _, _ = cycle.core_bounds
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    elif cycle.core_size is not None:
        cx, cy, _ = cycle.core_size
        ax.set_xlim(0.0, cx)
        ax.set_ylim(0.0, cy)
    else:
        ax.autoscale_view()

    plt.tight_layout()

    if outpath is not None:
        os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
        plt.savefig(outpath, dpi=600)

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


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot a C5G7 pin/mod tally slice.")
    parser.add_argument(
        "--infile",
        default="flux_Tally_reference2003_20260729_145023.txt",
        help="Tally file path or a filename located in resultPlot/ or results/.",
    )
    parser.add_argument("--out-dir", default="plots_Total", help="Output directory.")
    parser.add_argument(
        "--z-position",
        type=float,
        default=None,
        help="Global z position in cm. Defaults to ZSlicePosition(cm) in the file, or 0.63 cm for an older file.",
    )
    parser.add_argument("--show", action="store_true", help="Show the plot window.")
    args = parser.parse_args()

    infile = resolve_input_path(args.infile)
    cyc = parse_single_cycle_file(str(infile))
    if not cyc.assemblies:
        raise RuntimeError(f"No assembly tally blocks were parsed from {infile}")

    z_position_cm = (
        args.z_position
        if args.z_position is not None
        else cyc.z_slice_position if cyc.z_slice_position is not None
        else 0.63
    )
    input_timestamp = _input_timestamp_token(infile)
    outpath = Path(args.out_dir) / (
        f"power01_cycle{cyc.cycle_no}_{_z_filename_token(z_position_cm)}_"
        f"input{input_timestamp}.png"
    )
    statistics = plot_cycle_z_position(
        cyc,
        z_position_cm=z_position_cm,
        outpath=str(outpath),
        show=args.show,
    )
    statistics_path = outpath.with_suffix(".pin_power.txt")
    write_pin_power_statistics(statistics_path, z_position_cm, statistics)
    print(f"Loaded: {infile}")
    print(f"Saved:  {outpath.resolve()}")
    print(f"Pin power statistics: {statistics_path.resolve()}")


if __name__ == "__main__":
    main()
