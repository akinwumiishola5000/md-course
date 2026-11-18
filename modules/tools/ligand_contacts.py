#!/usr/bin/env python
"""
ligand_contacts.py - which protein residues does the ligand touch, and how often?

The Insilico Lab MD Course.

For every frame of the trajectory, this script finds the protein residues
with any heavy atom within a cutoff distance of any ligand heavy atom, and
the hydrogen bonds between protein and ligand. It then reports how often
(in % of frames) each contact and each hydrogen bond is present.

Usage:
    python ligand_contacts.py complex.tpr md_fit.xtc --ligand JZ4
    python ligand_contacts.py complex.tpr md_fit.xtc --ligand JZ4 --cutoff 4.5 --top 20

Outputs (in the current folder):
    contacts.csv        contact frequency for every residue that touches the ligand
    contacts.png        bar chart of the most frequent contacts
    hbonds.csv          occupancy of each protein-ligand hydrogen bond
"""
import argparse
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


def main():
    parser = argparse.ArgumentParser(description="Protein-ligand contact analysis.")
    parser.add_argument("topology", help="run input or structure file (e.g. complex.tpr)")
    parser.add_argument("trajectory", help="trajectory (e.g. md_fit.xtc)")
    parser.add_argument("--ligand", default="JZ4", help="ligand residue name (default JZ4)")
    parser.add_argument("--cutoff", type=float, default=4.0,
                        help="contact distance in Angstrom (default 4.0)")
    parser.add_argument("--top", type=int, default=15,
                        help="number of residues shown in the bar chart (default 15)")
    args = parser.parse_args()

    universe = mda.Universe(args.topology, args.trajectory)
    ligand = universe.select_atoms(f"resname {args.ligand} and not name H*")
    protein = universe.select_atoms("protein and not name H*")
    if len(ligand) == 0:
        sys.exit(f"Error: no atoms with residue name {args.ligand} were found.")
    n_frames = len(universe.trajectory)
    print(f"Ligand {args.ligand}: {len(ligand)} heavy atoms.  "
          f"Trajectory: {n_frames} frames.  Contact cutoff: {args.cutoff} Å")

    # ---- Contacts: residue counted once per frame if any heavy atom is close
    counts = {}
    labels = {}
    for _ in universe.trajectory:
        distances = distance_array(protein.positions, ligand.positions,
                                   box=universe.dimensions)
        close_atoms = protein[np.any(distances < args.cutoff, axis=1)]
        for residue in close_atoms.residues:
            key = residue.resindex
            counts[key] = counts.get(key, 0) + 1
            labels[key] = f"{residue.resname}{residue.resid}"

    if not counts:
        print("No protein residue came within the cutoff of the ligand.")
        print("The ligand may have left the binding pocket. Check its RMSD and distance.")
        return

    table = pd.DataFrame(
        [(labels[k], 100.0 * v / n_frames) for k, v in counts.items()],
        columns=["residue", "contact_percent"],
    ).sort_values("contact_percent", ascending=False)
    table.to_csv("contacts.csv", index=False, float_format="%.1f")

    shown = table.head(args.top)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(shown["residue"], shown["contact_percent"], color="#3b7dd8")
    ax.set_ylabel("Frames in contact (%)")
    ax.set_xlabel("Residue")
    ax.set_title(f"Protein residues within {args.cutoff} Å of {args.ligand}")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3)
    plt.xticks(rotation=60, ha="right")
    fig.tight_layout()
    fig.savefig("contacts.png", dpi=150)

    print("\nResidue contacts (% of frames):")
    for _, row in shown.iterrows():
        print(f"   {row['residue']:<8} {row['contact_percent']:5.1f} %")
    print("Saved: contacts.csv, contacts.png")

    # ---- Hydrogen bonds (donor-acceptor distance <= 3.0 A, angle >= 150 deg)
    try:
        hbonds = HydrogenBondAnalysis(universe,
                                      between=[f"resname {args.ligand}", "protein"],
                                      d_a_cutoff=3.0, d_h_a_angle_cutoff=150)
        hbonds.hydrogens_sel = hbonds.guess_hydrogens(
            f"resname {args.ligand} or protein")
        hbonds.acceptors_sel = hbonds.guess_acceptors(
            f"resname {args.ligand} or protein")
        hbonds.run()
    except Exception as error:          # e.g. topology without charges or bonds
        print(f"\nHydrogen-bond analysis skipped: {error}")
        return

    # Each row: frame, donor index, hydrogen index, acceptor index, distance, angle.
    # A donor-acceptor pair counts once per frame, whichever hydrogen is involved.
    pairs = {}
    for frame, donor, _, acceptor, _, _ in hbonds.results.hbonds:
        pairs.setdefault((int(donor), int(acceptor)), set()).add(int(frame))
    rows = []
    for (donor, acceptor), frames in pairs.items():
        d, a = universe.atoms[donor], universe.atoms[acceptor]
        rows.append((f"{d.resname}{d.resid}:{d.name}", f"{a.resname}{a.resid}:{a.name}",
                     100.0 * len(frames) / n_frames))
    hb_table = pd.DataFrame(rows, columns=["donor", "acceptor", "occupancy_percent"])
    hb_table = hb_table.sort_values("occupancy_percent", ascending=False)
    hb_table.to_csv("hbonds.csv", index=False, float_format="%.1f")

    print("\nProtein-ligand hydrogen bonds (% of frames):")
    if hb_table.empty:
        print("   none found")
    for _, row in hb_table.iterrows():
        print(f"   {row['donor']:<16} -> {row['acceptor']:<16} {row['occupancy_percent']:5.1f} %")
    print("Saved: hbonds.csv")


if __name__ == "__main__":
    main()
