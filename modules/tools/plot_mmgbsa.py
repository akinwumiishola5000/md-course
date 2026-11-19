#!/usr/bin/env python
"""
plot_mmgbsa.py - plot gmx_MMPBSA results as PNG images.

The Insilico Lab MD Course.

Run it in the folder where gmx_MMPBSA wrote its results:
    python plot_mmgbsa.py
    python plot_mmgbsa.py --top 20

Reads:
    FINAL_RESULTS_MMPBSA.dat   the binding energy and its components
    FINAL_DECOMP_MMPBSA.dat    per-residue contributions (if decomposition was run)
Writes:
    mmgbsa_terms.png           bar chart of the energy components
    mmgbsa_residues.png        the residues contributing most to binding
    mmgbsa_residues.csv        all per-residue contributions as a table
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TERMS = [
    ("ΔVDWAALS", "van der Waals"),
    ("ΔEEL", "Electrostatic"),
    ("ΔEGB", "Polar solvation"),
    ("ΔESURF", "Non-polar solvation"),
    ("ΔTOTAL", "ΔG binding"),
]


def read_terms(path):
    """Average and SD of each term in the 'Delta' block of FINAL_RESULTS."""
    values, in_delta = {}, False
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("Delta (Complex - Receptor - Ligand)"):
                in_delta = True
                continue
            if in_delta:
                parts = line.split()
                if parts and parts[0] in dict(TERMS):
                    values[parts[0]] = (float(parts[1]), float(parts[3]))
                if parts and parts[0] == "ΔTOTAL":
                    break
    return values


def read_residues(path):
    """Residue name, average TOTAL and its SD from the DELTAS block of FINAL_DECOMP."""
    rows, in_delta, header_lines = [], False, 0
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("DELTAS"):
                in_delta = True
                continue
            if not in_delta:
                continue
            if line.startswith(("Total Energy Decomposition", "Residue", ",Avg.")):
                header_lines += 1
                continue
            if not line:
                if rows:
                    break
                continue
            fields = next(csv.reader([line]))
            if len(fields) < 38:
                continue
            # 7 numbers per term: Internal, vdW, Elec, Polar, Non-polar, TOTAL
            rows.append((fields[0], float(fields[36]), float(fields[38])))
    return rows


def short_name(label):
    """'R:A:LEU:84' -> 'LEU84' ; 'L:B:JZ4:165' -> 'JZ4 (ligand)'."""
    parts = label.split(":")
    if len(parts) == 4:
        name = f"{parts[2]}{parts[3]}"
        return name + " (ligand)" if parts[0] == "L" else name
    return label


def main():
    parser = argparse.ArgumentParser(description="Plot gmx_MMPBSA results.")
    parser.add_argument("--top", type=int, default=15,
                        help="number of residues in the per-residue chart (default 15)")
    args = parser.parse_args()

    if not os.path.exists("FINAL_RESULTS_MMPBSA.dat"):
        sys.exit("Error: FINAL_RESULTS_MMPBSA.dat not found. Run this script in the "
                 "folder where gmx_MMPBSA finished.")

    # ---- 1. Energy components
    terms = read_terms("FINAL_RESULTS_MMPBSA.dat")
    labels = [nice for key, nice in TERMS if key in terms]
    means = [terms[key][0] for key, _ in TERMS if key in terms]
    sds = [terms[key][1] for key, _ in TERMS if key in terms]
    colors = ["#3b7dd8" if m < 0 else "#d8573b" for m in means]
    colors[-1] = "#2a2a2a"

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, means, yerr=sds, capsize=4, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("MM/GBSA binding energy and its components")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig("mmgbsa_terms.png", dpi=150)

    print("MM/GBSA energy components (kcal/mol, average ± SD over frames):")
    for (key, nice) in TERMS:
        if key in terms:
            print(f"   {nice:<22} {terms[key][0]:8.2f} ± {terms[key][1]:.2f}")
    print("Saved: mmgbsa_terms.png")

    # ---- 2. Per-residue contributions
    if not os.path.exists("FINAL_DECOMP_MMPBSA.dat"):
        print("No FINAL_DECOMP_MMPBSA.dat: per-residue decomposition was not run.")
        return
    residues = read_residues("FINAL_DECOMP_MMPBSA.dat")
    if not residues:
        print("Could not read per-residue results from FINAL_DECOMP_MMPBSA.dat.")
        return

    with open("mmgbsa_residues.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["residue", "total_kcal_mol", "sd"])
        for name, total, sd in residues:
            writer.writerow([short_name(name), f"{total:.3f}", f"{sd:.3f}"])

    protein = [r for r in residues if not r[0].startswith("L:")]
    shown = sorted(protein, key=lambda r: r[1])[: args.top]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar([short_name(r[0]) for r in shown], [r[1] for r in shown],
           yerr=[r[2] for r in shown], capsize=3,
           color=["#3b7dd8" if r[1] < 0 else "#d8573b" for r in shown])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Contribution to ΔG (kcal/mol)")
    ax.set_xlabel("Residue")
    ax.set_title("Residues contributing most to binding")
    ax.grid(axis="y", alpha=0.3)
    plt.xticks(rotation=60, ha="right")
    fig.tight_layout()
    fig.savefig("mmgbsa_residues.png", dpi=150)

    print(f"\nTop {len(shown)} protein residues (kcal/mol; negative = favours binding):")
    for name, total, sd in shown:
        print(f"   {short_name(name):<10} {total:7.2f} ± {sd:.2f}")
    print("Saved: mmgbsa_residues.png, mmgbsa_residues.csv")


if __name__ == "__main__":
    main()
