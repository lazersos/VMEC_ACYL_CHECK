#!/usr/bin/env python3
"""Plot the CSV produced by compare_wout at one toroidal plane.

Examples:
  python plot_comparison.py comparison.csv --phi-index 1 -o phi1.png
  python plot_comparison.py comparison.csv --phi-index 8 --component magnitude
  python plot_comparison.py comparison.csv --list-phi

Indices are ONE BASED over the sorted phi planes present in the CSV.
The duplicated periodic endpoint is absent from compare_wout's CSV.
"""

import argparse
import csv
from pathlib import Path

import numpy as np


def uniform_axis(values):
    """Reconstruct the uniform grid, retaining gaps as empty grid cells."""
    axis = np.unique(values)
    if axis.size < 2:
        raise ValueError("Need at least two distinct R and Z coordinates to plot a plane.")
    step = np.min(np.diff(axis))
    count = int(round((axis[-1] - axis[0]) / step)) + 1
    if count > 10 * axis.size:
        raise ValueError("Cannot infer a reasonable uniform R/Z grid from the CSV.")
    grid = np.linspace(axis[0], axis[-1], count)
    indices = np.rint((axis - grid[0]) / (grid[1] - grid[0])).astype(int)
    if not np.allclose(axis, grid[indices], rtol=1e-10, atol=1e-12):
        raise ValueError("R and Z coordinates must lie on a uniform grid.")
    return grid


def load_plane(filename, phi_index, list_phi=False):
    names = ["R", "phi", "Z", "Br_getBcyl", "Bphi_getBcyl", "Bz_getBcyl",
             "Br_curlA", "Bphi_curlA", "Bz_curlA"]
    with Path(filename).open(newline="") as stream:
        header = [name.strip() for name in next(csv.reader(stream))]
    missing = set(names) - set(header)
    if missing:
        raise ValueError("Missing CSV columns: " + ", ".join(sorted(missing)))
    data = np.loadtxt(filename, delimiter=",", skiprows=1,
                      usecols=[header.index(name) for name in names], ndmin=2)
    if data.shape[0] == 0:
        raise ValueError("CSV contains no comparison points.")
    if not np.all(np.isfinite(data[:, :3])):
        raise ValueError("CSV has nonfinite coordinates.")
    phis = np.unique(data[:, 1])
    if list_phi:
        print("Index (one based)       phi [rad]          phi [deg]")
        for index, phi in enumerate(phis, 1):
            print(f"{index:8d}          {phi:17.10g}   {np.degrees(phi):17.10g}")
        return None
    if not 1 <= phi_index <= len(phis):
        raise ValueError(f"--phi-index must be between 1 and {len(phis)}.")
    phi = phis[phi_index - 1]
    plane = data[data[:, 1] == phi]
    # Use coordinates from ALL planes so an absent row or column in the
    # selected plane stays blank rather than being interpolated across.
    r = uniform_axis(data[:, 0])
    z = uniform_axis(data[:, 2])
    ir = np.rint((plane[:, 0] - r[0]) / (r[1] - r[0])).astype(int)
    iz = np.rint((plane[:, 2] - z[0]) / (z[1] - z[0])).astype(int)
    if np.unique(iz * len(r) + ir).size != len(plane):
        raise ValueError("Repeated R/Z coordinates in the selected phi plane.")
    reference = np.full((len(z), len(r), 3), np.nan)
    curl = np.full_like(reference, np.nan)
    reference[iz, ir] = plane[:, 3:6]
    curl[iz, ir] = plane[:, 6:9]
    return r, z, phi, len(phis), reference, curl


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--phi-index", type=int, default=1,
                        help="one-based index among sorted phi planes in CSV (default: 1)")
    parser.add_argument("--list-phi", action="store_true", help="list available plane indices and exit")
    parser.add_argument("--component", choices=["all", "R", "phi", "Z", "magnitude"],
                        default="all", help="rows to plot (default: all components and magnitude)")
    parser.add_argument("-o", "--output", type=Path, help="output image or PDF filename")
    parser.add_argument("--log-error", action="store_true", help="use logarithmic relative-error colors")
    parser.add_argument("--show", action="store_true", help="also open the interactive figure")
    args = parser.parse_args()
    try:
        result = load_plane(args.csv_file, args.phi_index, args.list_phi)
    except (OSError, ValueError, StopIteration) as exc:
        parser.error(str(exc))
    if result is None:
        return
    r, z, phi, nplanes, reference, curl = result

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    refmag = np.linalg.norm(reference, axis=2)
    curlmag = np.linalg.norm(curl, axis=2)
    difference = curl - reference
    denominator = np.where(refmag > 0, refmag, np.nan)
    relative = np.linalg.norm(difference, axis=2) / denominator
    finite = np.isfinite(relative)
    if not np.any(finite):
        parser.error("Selected plane has no finite comparisons with nonzero reference field.")

    components = ["R", "phi", "Z", "magnitude"] if args.component == "all" else [args.component]
    fig, axes = plt.subplots(len(components), 3, figsize=(12, 3.3 * len(components)),
                             squeeze=False, layout="constrained")
    fig.suptitle(f"Toroidal plane {args.phi_index}/{nplanes}: "
                 rf"$\phi={phi:.7g}$ rad ({np.degrees(phi):.4f} deg)", fontsize=13)
    for row, component in enumerate(components):
        if component == "magnitude":
            left, middle, error = refmag, curlmag, relative
            label = r"$|\mathbf{B}|$"
            error_label = r"$|\mathbf{B}_{curl}-\mathbf{B}_{ref}|/|\mathbf{B}_{ref}|$"
            colors = "viridis"
            vmin = min(np.nanmin(left), np.nanmin(middle))
            vmax = max(np.nanmax(left), np.nanmax(middle))
        else:
            index = {"R": 0, "phi": 1, "Z": 2}[component]
            left, middle = reference[:, :, index], curl[:, :, index]
            error = np.abs(difference[:, :, index]) / denominator
            symbol = r"\phi" if component == "phi" else component
            label = rf"$B_{{{symbol}}}$"
            error_label = rf"$|\Delta B_{{{symbol}}}|/|\mathbf{{B}}_{{ref}}|$"
            colors = "RdBu_r"
            vmax = max(np.nanmax(np.abs(left)), np.nanmax(np.abs(middle)))
            vmin = -vmax
        if vmax <= vmin:
            padding = max(abs(vmax) * 0.01, 1e-12)
            vmin -= padding
            vmax += padding
        shared_norm = Normalize(vmin=vmin, vmax=vmax)
        for col, (values, method) in enumerate([(left, "GetBcyl"), (middle, "curl(GetAcyl)")]):
            image = axes[row, col].pcolormesh(r, z, np.ma.masked_invalid(values),
                                              shading="nearest", cmap=colors, norm=shared_norm)
            axes[row, col].set_title(f"{method}: {label}")
            fig.colorbar(image, ax=axes[row, col], label="T", shrink=0.85)
        error_values = np.ma.masked_invalid(error)
        error_norm = None
        if args.log_error:
            positive = error[np.isfinite(error) & (error > 0)]
            if positive.size:
                vmax_error = float(np.max(positive))
                vmin_error = max(float(np.min(positive)), vmax_error * 1e-6)
                if vmin_error >= vmax_error:
                    vmin_error = vmax_error / 10
                error_norm = LogNorm(vmin=vmin_error, vmax=vmax_error)
                # Exact zero errors get the minimum color, not an empty cell.
                error_values = np.ma.maximum(error_values, vmin_error)
        image = axes[row, 2].pcolormesh(r, z, error_values, shading="nearest",
                                        cmap="magma", norm=error_norm,
                                        **({} if error_norm else {"vmin": 0}))
        axes[row, 2].set_title(error_label)
        fig.colorbar(image, ax=axes[row, 2], label="Relative error", shrink=0.85)
        for ax in axes[row]:
            ax.set(xlabel="R [m]", ylabel="Z [m]", aspect="equal")
            ax.set_facecolor("0.92")

    output = args.output or Path(f"{args.csv_file.stem}_phi{args.phi_index:04d}.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    relrms = np.sqrt(np.sum(np.sum(difference[finite]**2, axis=1)) / np.sum(refmag[finite]**2))
    print(f"Selected phi index {args.phi_index}/{nplanes}: {phi:.12g} rad")
    print(f"Compared points: {np.count_nonzero(finite)}")
    print(f"Relative RMS vector error: {relrms:.7e}")
    print(f"Maximum relative vector error: {np.max(relative[finite]):.7e}")
    print(f"Saved: {output}")
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
