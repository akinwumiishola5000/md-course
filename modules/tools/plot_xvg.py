#!/usr/bin/env python
"""
plot_xvg.py - plot GROMACS .xvg files as PNG images.

The Insilico Lab MD Course.

Usage examples:
    python plot_xvg.py rmsd.xvg
    python plot_xvg.py temperature.xvg --avg 10
    python plot_xvg.py rmsd_protein.xvg rmsd_ligand.xvg -o rmsd_compare.png
    python plot_xvg.py rmsd.xvg --unit A          (convert nm to Angstrom)
    python plot_xvg.py potential.xvg --xlabel "Step"

The title, axis labels and legend are read from the .xvg file.
The PNG is saved in the same folder; click it in the file panel to view it.
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")              # no screen needed: save to file
import matplotlib.pyplot as plt
import numpy as np


def clean(label):
    """Turn xmgrace formatting codes (\\S2\\N, \\f{Symbol}...) into plain text."""
    label = re.sub(r"\\S(.*?)\\N", r"^\1", label)
    label = re.sub(r"\\s(.*?)\\N", r"_\1", label)
    label = re.sub(r"\\f\{.*?\}|\\[a-zA-Z]", "", label)
    return label.strip()


def read_xvg(path):
    """Return (data array, title, xlabel, ylabel, legends) from an .xvg file."""
    title, xlabel, ylabel, legends, rows = "", "", "", {}, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("@"):
                quoted = re.findall(r'"(.*)"', line)
                text = clean(quoted[0]) if quoted else ""
                if line.startswith("@    title") or line.startswith("@ title"):
                    title = text
                elif "xaxis" in line and "label" in line:
                    xlabel = text
                elif "yaxis" in line and "label" in line:
                    ylabel = text
                else:
                    match = re.match(r"@\s*s(\d+)\s+legend", line)
                    if match:
                        legends[int(match.group(1))] = text
                continue
            try:
                rows.append([float(value) for value in line.split()])
            except ValueError:
                continue
    if not rows:
        sys.exit(f"Error: no data found in {path}")
    width = min(len(row) for row in rows)
    data = np.array([row[:width] for row in rows])
    return data, title, xlabel, ylabel, legends


def running_average(values, window):
    if window <= 1 or window > len(values):
        return None
    kernel = np.ones(window) / window
    return np.convolve(values, kernel, mode="valid")


def main():
    parser = argparse.ArgumentParser(description="Plot GROMACS .xvg files.")
    parser.add_argument("files", nargs="+", help="one or more .xvg files")
    parser.add_argument("-o", "--output", help="output PNG name")
    parser.add_argument("--avg", type=int, default=0,
                        help="also draw a running average over N points")
    parser.add_argument("--unit", choices=["nm", "A"], default="nm",
                        help="A converts y values from nm to Angstrom")
    parser.add_argument("--title", help="replace the title")
    parser.add_argument("--xlabel", help="replace the x-axis label")
    parser.add_argument("--ylabel", help="replace the y-axis label")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(8, 4.5))
    title, xlabel, ylabel = "", "", ""
    for path in args.files:
        data, title, xlabel, ylabel, legends = read_xvg(path)
        x = data[:, 0]
        for col in range(1, data.shape[1]):
            y = data[:, col]
            if args.unit == "A":
                y = y * 10.0
            name = legends.get(col - 1, "")
            if len(args.files) > 1:
                name = f"{os.path.splitext(os.path.basename(path))[0]} {name}".strip()
            line, = ax.plot(x, y, linewidth=1.2, alpha=0.6 if args.avg else 1.0,
                            label=name or None)
            smooth = running_average(y, args.avg)
            if smooth is not None:
                offset = args.avg // 2
                ax.plot(x[offset:offset + len(smooth)], smooth, color=line.get_color(),
                        linewidth=2.2, label=f"{name} (running avg, {args.avg})".strip())
            print(f"{os.path.basename(path)}  {name or 'column ' + str(col)}: "
                  f"mean = {y.mean():.4g}, std = {y.std():.4g}, "
                  f"min = {y.min():.4g}, max = {y.max():.4g}")

    if args.unit == "A":
        ylabel = ylabel.replace("(nm)", "(Å)").replace("nm", "Å")
    ax.set_title(args.title or title)
    ax.set_xlabel(args.xlabel or xlabel)
    ax.set_ylabel(args.ylabel or ylabel)
    ax.grid(alpha=0.3)
    if ax.get_legend_handles_labels()[1]:
        ax.legend(fontsize=9)
    fig.tight_layout()

    output = args.output or os.path.splitext(args.files[0])[0] + ".png"
    fig.savefig(output, dpi=150)
    print(f"Saved plot: {output}  (click it in the file panel to view)")


if __name__ == "__main__":
    main()
