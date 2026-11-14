# Module 02: Simulating a Protein-Ligand Complex

In drug discovery, we often want to know whether a small molecule stays
bound to its target protein, and how it interacts with the binding
pocket. In this module you simulate **T4 lysozyme (L99A/M102Q mutant)
bound to the ligand JZ4 (2-propylphenol)**, PDB entry **3HTB**.

The new skill in this module is **preparing the ligand**: the protein
force field knows nothing about small molecules, so you will build
the ligand's parameters yourself with GAFF2 and ACPYPE, then join them
to the protein.

**Time needed:** about 3.5 hours, including about 40 minutes of
computing time that runs in the background
**Before you start:** complete Modules 00 and 01

---

## Before the session: install a molecular viewer on your computer

Your codespace has no screen for 3D graphics, so you will **download**
structures from the codespace and look at them on your own computer.
Install at least one of these before the session. All are free for
academic use.

| Program | Best for | Get it from |
|---|---|---|
| **PyMOL** | Looking at structures, binding pockets, making figures | pymol.org (free educational licenses are available; an open-source version also exists) |
| **VMD** | Playing MD trajectories (movies of your simulation) | ks.uiuc.edu/Research/vmd (free registration) |
| **UCSF ChimeraX** | Structures and trajectories, easy to learn | cgl.ucsf.edu/chimerax (free for academic use) |
| **Mol\*** (no installation) | A quick look at a structure in your web browser | molstar.org/viewer |

**Recommended:** PyMOL for structures and VMD for trajectories. If you
cannot install anything, use Mol\* in your browser for structures.

### How to download a file from your codespace

In the file panel on the left, **right-click** the file and choose
**Download**. Your browser saves it to your Downloads folder.

Throughout this module, a **👁 View it** box tells you which file to
download, which program to open it in, and what to look for.

---

## The workflow

```
3HTB.pdb ─┬─ protein ──→ protein topology (pdb2gmx) ─────────────┐
          │                                                       ├→ complex → box → water → ions
          └─ JZ4 ligand → add hydrogens → ligand topology (ACPYPE)┘
                                                                          ↓
   analysis ← production MD ← NPT ← NVT ← index groups ← energy minimization
```

| Part | What you do | Main output |
|---|---|---|
| 1 | Download and inspect 3HTB | `3HTB.pdb` |
| 2 | Separate the protein and the ligand | `protein.pdb`, `JZ4.pdb` |
| 3 | Add hydrogens to the ligand | `JZ4_H.mol2` |
| 4 | Build the ligand topology | `JZ4_GMX.itp`, `JZ4_GMX.gro` |
| 5 | Build the protein topology | `protein.gro`, `topol.top` |
| 6 | Add the ligand to the topology | edited `topol.top` |
| 7 | Join the protein and ligand coordinates | `complex.gro` |
| 8 | Box, water and ions | `solv_ions.gro` |
| 9 | Energy minimization | `em.gro` |
| 10 | Create index groups | `index.ndx` |
| 11–13 | NVT, NPT and production MD | `md.xtc` |
| 14 | Analysis and plots | `.xvg` and `.png` files |
| 15 | Watch the trajectory | `complex_traj.xtc` |

---

## The files you will meet

You will use many file types. Come back to this table whenever you are
unsure what a file is.

| Extension | What it is | Can you read it? |
|---|---|---|
| `.pdb` | Protein Data Bank structure: atoms and coordinates in Å | Yes, text |
| `.mol2` | Small-molecule structure with atom types, **bonds** and charges | Yes, text |
| `.gro` | GROMACS structure: atoms and coordinates in **nm** | Yes, text |
| `.top` | The **topology** of the whole system: which molecules, how many, and their parameters | Yes, text |
| `.itp` | An **included topology**: the parameters of one molecule, or restraints, used inside `.top` | Yes, text |
| `.mdp` | Simulation **settings** (time step, temperature, length...) | Yes, text |
| `.tpr` | **Run input**: structure + topology + settings combined by `grompp` | No, binary |
| `.ndx` | **Index file**: named groups of atom numbers | Yes, text |
| `.cpt` | **Checkpoint**: positions and velocities, to continue or restart a run | No, binary |
| `.edr` | **Energy file**: energies, temperature, pressure over time | No, use `gmx energy` |
| `.xtc` | **Trajectory**: compressed coordinates over time (the "movie") | No, use a viewer |
| `.log` | Log of a run, including performance | Yes, text |
| `.xvg` | Data for plotting (time and values) from GROMACS analysis tools | Yes, text |

---

## Part 0: Get ready

Get the latest course files:

```bash
cd /workspaces/md-course
```

```bash
git pull
```

Make a folder for this module and move into it:

```bash
mkdir -p my-work/02-complex
```

```bash
cd my-work/02-complex
```

Copy the simulation settings and the plotting tool:

```bash
cp /workspaces/md-course/modules/02-protein-ligand/files/*.mdp .
```

```bash
cp /workspaces/md-course/modules/tools/plot_xvg.py .
```

```bash
ls
```

You should see `ions.mdp`, `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp` and
`plot_xvg.py`.

### About the plotting tool

GROMACS writes analysis results as `.xvg` files: plain columns of
numbers. `plot_xvg.py` turns them into pictures and prints the mean,
standard deviation, minimum and maximum. You will use it throughout
this module:

```bash
python plot_xvg.py file.xvg
```

It saves `file.png` in the same folder. **Click the `.png` file in the
file panel to open it.** Useful options:

| Option | Effect |
|---|---|
| `--unit A` | convert nm to Å (most papers report RMSD in Å) |
| `--avg 10` | add a running average over 10 points, to see the trend |
| `-o name.png` | choose the output file name |
| `--xlabel "..."`, `--ylabel "..."`, `--title "..."` | change the labels |
| two `.xvg` files | plot both on the same graph, to compare them |

---

## Part 1: Download and inspect 3HTB

```bash
wget https://files.rcsb.org/download/3HTB.pdb
```

What non-protein molecules does it contain?

```bash
grep "^HETATM" 3HTB.pdb | cut -c18-20 | sort | uniq -c
```

You should see `JZ4` (our ligand), `HOH` (crystal water), and other
small molecules from the crystallization conditions. We keep **only
the protein and JZ4**.

> **👁 View it: the binding pocket**
> Download `3HTB.pdb` and open it in **PyMOL**. Type these commands in
> the PyMOL command line (the `PyMOL>` box), one at a time:
> ```
> hide everything
> show cartoon
> select lig, resn JZ4
> show sticks, lig
> select pocket, byres (polymer within 5 of lig)
> show sticks, pocket
> zoom lig
> ```
> The ligand sits in a buried, mostly hydrophobic pocket. Which pocket
> residues are near the ligand's oxygen atom?
>
> In **Mol\*** (no installation): open molstar.org/viewer, drag
> `3HTB.pdb` into the window, and click on the ligand to zoom to it.

---

## Part 2: Separate the protein and the ligand

The protein and the ligand are prepared with different tools, so they
must be in separate files.

### The protein

Keep only the protein atoms (`ATOM` lines):

```bash
grep "^ATOM" 3HTB.pdb > protein_raw.pdb
```

Some crystal structures contain **alternate locations (altLocs)**: two
possible positions for the same side chain, labeled `A` and `B` in
column 17. GROMACS needs exactly one position per atom. Check:

```bash
cut -c17 protein_raw.pdb | sort | uniq -c
```

If you only see a count next to a blank, there are no altLocs. If you
also see `A` and `B`, keep the blank and `A` positions only. This
command is safe to run either way:

```bash
awk 'substr($0,17,1)==" " || substr($0,17,1)=="A"' protein_raw.pdb > protein.pdb
```

### The ligand

```bash
grep "^HETATM" 3HTB.pdb | grep JZ4 > JZ4_raw.pdb
```

Apply the same altLoc rule to the ligand:

```bash
awk 'substr($0,17,1)==" " || substr($0,17,1)=="A"' JZ4_raw.pdb > JZ4.pdb
```

```bash
cat JZ4.pdb
```

Count the atoms. JZ4 (2-propylphenol, C₉H₁₂O) has 10 heavy atoms: 9
carbons and 1 oxygen. Crystal structures almost never contain hydrogen
atoms, because X-ray diffraction cannot see them clearly. You add them
in the next part.

---

## Part 3: Add hydrogens to the ligand

```bash
obabel JZ4.pdb -O JZ4_H.mol2 -h
```

| Part of the command | Meaning |
|---|---|
| `obabel` | Open Babel, a tool for converting chemical file formats |
| `JZ4.pdb` | input file |
| `-O JZ4_H.mol2` | output file, in MOL2 format |
| `-h` | add hydrogen atoms |

Open Babel works out the bonds from the distances between atoms, then
adds hydrogens to fill each atom's valence. The heavy atoms keep their
exact crystal positions, so the ligand stays in the binding pocket.

Check the chemical formula:

```bash
obabel JZ4_H.mol2 -otxt --append formula
```

It must print **C9H12O**. If it does not, the bonds or hydrogens are
wrong, and the ligand must not be used.

### Understanding the MOL2 file

```bash
cat JZ4_H.mol2
```

A MOL2 file has sections that start with `@<TRIPOS>`:

| Section | Contains |
|---|---|
| `@<TRIPOS>MOLECULE` | The molecule name and the number of atoms and bonds |
| `@<TRIPOS>ATOM` | One line per atom: number, name, x y z (Å), atom type (`C.ar` = aromatic carbon, `C.3` = sp3 carbon, `O.3` = sp3 oxygen, `H`), residue, charge |
| `@<TRIPOS>BOND` | One line per bond: bond number, the two atoms, and the bond type (`1` single, `2` double, `ar` aromatic) |

Unlike a PDB file, a MOL2 file stores **bonds and bond orders**, which
the next tool needs.

> **⚠ In your own research:** adding hydrogens is where most ligand
> errors happen. Drug-like molecules often have groups that are charged
> at pH 7 (amines, carboxylic acids). Always check the formula and the
> protonation state, and set the correct net charge in Part 4.

> **👁 View it: the ligand with hydrogens**
> Download `JZ4_H.mol2` and open it in PyMOL. Check that there are
> hydrogens on every carbon, one on the oxygen (the phenol OH), and that
> the ring is aromatic.

---

## Part 4: Build the ligand topology with ACPYPE

The protein force field (AMBER99SB-ILDN) has parameters only for amino
acids, water and ions. For the ligand we use **GAFF2** (the General
AMBER Force Field 2), which was designed for drug-like molecules and is
compatible with AMBER protein force fields.

**ACPYPE** runs the AmberTools programs that assign GAFF2 parameters and
charges, then writes the results in GROMACS format:

```bash
acpype -i JZ4_H.mol2 -c bcc -n 0 -a gaff2 -b JZ4
```

| Option | Meaning |
|---|---|
| `-i JZ4_H.mol2` | input ligand, with hydrogens |
| `-c bcc` | calculate **AM1-BCC partial charges**, a fast quantum-chemistry method that gives charges compatible with AMBER force fields |
| `-n 0` | the **net charge** of the molecule. JZ4 is neutral. For a charged ligand this must be set correctly (for example `-n 1` or `-n -1`) |
| `-a gaff2` | use the GAFF2 atom types and parameters |
| `-b JZ4` | base name for the output files, and the molecule name |

This takes about a minute. ACPYPE creates a folder called `JZ4.acpype`:

```bash
ls JZ4.acpype
```

You need three of these files. Copy them, one at a time:

```bash
cp JZ4.acpype/JZ4_GMX.itp .
```

```bash
cp JZ4.acpype/JZ4_GMX.gro .
```

```bash
cp JZ4.acpype/posre_JZ4.itp .
```

| File | Contains |
|---|---|
| `JZ4_GMX.itp` | the ligand topology: every parameter GROMACS needs for JZ4 |
| `JZ4_GMX.gro` | the ligand coordinates with hydrogens, in nm, still in the binding-pocket position |
| `posre_JZ4.itp` | position restraints for the ligand, used during equilibration |

### Understanding the ligand topology

Open the topology and read along with the explanations below:

```bash
less JZ4_GMX.itp
```

(Press **q** to quit `less`.) The file is divided into sections, each
starting with a name in square brackets:

| Section | What it describes |
|---|---|
| `[ atomtypes ]` | The GAFF2 atom types used by this molecule, with their van der Waals parameters (sigma and epsilon). For example `ca` = aromatic carbon, `c3` = sp3 carbon, `oh` = hydroxyl oxygen, `ho` = hydroxyl hydrogen, `ha` = aromatic hydrogen, `hc` = aliphatic hydrogen |
| `[ moleculetype ]` | The molecule's name (`JZ4`). This is the name you use in `topol.top` |
| `[ atoms ]` | One line per atom: number, atom type, residue, atom name, **partial charge** and mass. The last column (`qtot`) is the running total of the charge |
| `[ bonds ]` | Each bonded pair, with its equilibrium length and force constant |
| `[ pairs ]` | Atom pairs three bonds apart (1-4 pairs), whose nonbonded interactions are scaled |
| `[ angles ]` | Each three-atom angle, with its equilibrium value and force constant |
| `[ dihedrals ]` | Torsion terms for rotation around bonds, and **impropers** that keep the aromatic ring flat |

Check that the total charge is zero. The last atom's `qtot` value must
match the net charge you gave with `-n`:

```bash
grep "qtot" JZ4_GMX.itp | tail -1
```

Look at the charges of the oxygen and its hydrogen:

```bash
grep -E " oh | ho " JZ4_GMX.itp
```

The oxygen is strongly negative and its hydrogen strongly positive.
This polar OH group is what can form hydrogen bonds with the protein.

---

## Part 5: Build the protein topology

This is the same as Module 01:

```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p -ff amber99sb-ildn -ignh
```

Check the protein's net charge:

```bash
grep "qtot" topol.top | tail -1
```

Write it down; you will neutralize it with ions in Part 8.

---

## Part 6: Add the ligand to the topology

`topol.top` currently describes only the protein. You need to tell it
about JZ4 in **two places**.

First, look at how the file starts:

```bash
head -n 25 topol.top
```

The order of a topology file matters:

1. `#include "amber99sb-ildn.ff/forcefield.itp"` loads the force field,
   including all **atom types** for proteins, water and ions.
2. `[ moleculetype ]` sections define each molecule.
3. `[ system ]` gives the system a name.
4. `[ molecules ]` lists every molecule and how many copies, **in the
   same order as the coordinate file**.

**Rule:** all `[ atomtypes ]` must be defined before the first
`[ moleculetype ]`. `JZ4_GMX.itp` starts with `[ atomtypes ]`, so it must
be included **directly after the force field line**.

Open the topology in nano:

```bash
nano topol.top
```

### Edit 1: include the ligand parameters

Find this line near the top:

```
#include "amber99sb-ildn.ff/forcefield.itp"
```

Put your cursor at the end of that line, press **Enter** to make a new
line, and type:

```
; Include ligand parameters
#include "JZ4_GMX.itp"

; Ligand position restraints
#ifdef POSRES
#include "posre_JZ4.itp"
#endif
```

The `#ifdef POSRES` block switches the ligand restraints on only when
an `.mdp` file contains `define = -DPOSRES`, which is true for NVT and
NPT and false for production.

### Edit 2: add the ligand to the molecule list

Press **Ctrl + End** (or hold the **down arrow**) to go to the end of
the file. You will see:

```
[ molecules ]
; Compound        #mols
Protein             1
```

Add a new line at the end:

```
JZ4                 1
```

Save with **Ctrl + O**, then **Enter**, and exit with **Ctrl + X**.

### Check your edits

```bash
grep -n "JZ4" topol.top
```

You should see three lines: the two `#include` lines near the top, and
`JZ4  1` at the end.

```bash
tail -n 4 topol.top
```

> **Note:** if your protein were split into several chains, the
> molecule names would be `Protein_chain_A`, `Protein_chain_B` and so
> on. For 3HTB there is one chain, so it is `Protein`.

---

## Part 7: Join the protein and ligand coordinates

You now have the protein in `protein.gro` and the ligand in
`JZ4_GMX.gro`. Both use the same coordinate frame (from the crystal
structure), so you can simply join them into one file.

### Understanding the .gro format

```bash
head -n 4 protein.gro
```

```bash
tail -n 2 protein.gro
```

A `.gro` file always has:

| Line | Contains |
|---|---|
| Line 1 | A title |
| Line 2 | The **number of atoms** |
| One line per atom | residue number, residue name, atom name, atom number, then **x, y, z in nm** |
| Last line | The **box size** in nm |

Now look at the ligand file the same way:

```bash
cat JZ4_GMX.gro
```

### Build complex.gro in three steps

Step 1: copy the protein file without its last line (the box):

```bash
head -n -1 protein.gro > complex.gro
```

Step 2: add the ligand atoms, skipping the ligand's title, atom count
and box lines:

```bash
tail -n +3 JZ4_GMX.gro | head -n -1 >> complex.gro
```

Step 3: add the box line back at the end:

```bash
tail -n 1 protein.gro >> complex.gro
```

Note the `>>` in steps 2 and 3: it **adds** to the file instead of
replacing it.

### Fix the atom count

Line 2 still contains the protein's atom count. The correct number is
the total number of lines minus 3 (the title, count and box lines):

```bash
wc -l complex.gro
```

Subtract 3 from that number. Then open the file:

```bash
nano complex.gro
```

Move to **line 2**, delete the old number and type the new one. Save
with **Ctrl + O**, **Enter**, then exit with **Ctrl + X**.

Check the start and end of the file:

```bash
head -n 3 complex.gro
```

```bash
tail -n 3 complex.gro
```

The last atom lines should be the JZ4 atoms, followed by the box.

> **Why must the order match?** GROMACS reads the coordinates in order
> and matches them to `[ molecules ]` in order: first all protein atoms,
> then the JZ4 atoms. If the order or the count is wrong, `grompp` stops
> with an error.

> **👁 View it: the prepared complex**
> Download `complex.gro` and open it in PyMOL (PyMOL reads `.gro`
> files). Type:
> ```
> show sticks, resn JZ4
> zoom resn JZ4
> ```
> The ligand should sit in the pocket, now with hydrogens, exactly
> where it was in the crystal structure.

---

## Part 8: Box, water and ions

These steps are the same as Module 01. Run them one at a time:

```bash
gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -c -d 1.0
```

```bash
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```

```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
```

As in Module 01, `-maxwarn 1` is only for the net-charge warning that
the next command fixes.

```bash
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```

```bash
tail -n 6 topol.top
```

You should see `Protein`, `JZ4`, `SOL`, `NA` and `CL`.

> **👁 View it: the solvated system**
> Download `solv_ions.gro` (about 1.5 MB) and open it in VMD or PyMOL.
> In PyMOL, type `show spheres, resn NA+CL` to see the ions.
> Notice the shape of the dodecahedral box filled with water.

---

## Part 9: Energy minimization

```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
```

```bash
gmx mdrun -deffnm em -nt 2
```

```bash
grep -E "converged|Potential Energy|Maximum force" em.log
```

### Plot the minimization

Extract the potential energy:

```bash
echo "Potential" | gmx energy -f em.edr -o potential.xvg
```

Plot it:

```bash
python plot_xvg.py potential.xvg --xlabel "Minimization step" --ylabel "Potential energy (kJ/mol)" --title "Energy minimization"
```

Click `potential.png` in the file panel.

**What you should see:** the energy falls steeply in the first steps as
the worst contacts are removed, then levels off. A flat end means the
minimization has converged.

> For energy minimization, GROMACS writes the step number in the
> "time" column, which is why we changed the x-axis label.

---

## Part 10: Create index groups

During equilibration, the temperature is controlled separately for the
**solute** (protein + ligand) and the **solvent** (water + ions). This
is more accurate than controlling them together.

The thermostat therefore needs a group containing the protein and the
ligand together. GROMACS does not make this group automatically, so
you create it with `gmx make_ndx`:

```bash
printf '"Protein" | "JZ4"\nq\n' | gmx make_ndx -f em.gro -o index.ndx
```

| Part | Meaning |
|---|---|
| `"Protein" \| "JZ4"` | make a new group containing atoms in Protein **or** JZ4 |
| `q` | quit and save |
| `printf '...\n...'` | sends these answers to `make_ndx`, one per line |

Check that the new group exists:

```bash
grep "^\[" index.ndx
```

The list should include `Protein_JZ4` and `Water_and_ions`. These are
the two groups used by `tc-grps` in the `.mdp` files:

```bash
grep tc-grps nvt.mdp
```

From now on, every `grompp` command needs `-n index.ndx` so it can find
these groups.

---

## Part 11: NVT equilibration (50 ps)

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
```

```bash
nohup gmx mdrun -deffnm nvt -nt 2 > nvt.out 2>&1 &
```

This takes about 5 to 10 minutes. Check progress with:

```bash
tail -n 15 nvt.log
```

When `jobs` shows nothing running, extract and plot the temperature:

```bash
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
```

```bash
python plot_xvg.py temperature.xvg --avg 10
```

**What you should see:** the temperature rises quickly to 300 K and then
fluctuates around it. The thick line (running average) should be flat
at about 300 K. Fluctuations of a few kelvin are normal.

---

## Part 12: NPT equilibration (50 ps)

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
```

```bash
nohup gmx mdrun -deffnm npt -nt 2 > npt.out 2>&1 &
```

When it finishes, extract the pressure and the density **separately**,
because they have very different scales:

```bash
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
```

```bash
echo "Density" | gmx energy -f npt.edr -o density.xvg
```

Plot each one:

```bash
python plot_xvg.py pressure.xvg --avg 10
```

```bash
python plot_xvg.py density.xvg --avg 10
```

**What you should see:**

- **Pressure** jumps up and down by hundreds of bar. That is normal for
  a system this small. Look at the **mean** printed in the terminal and
  the running average: they should be around 1 bar, within the large
  fluctuation.
- **Density** should settle around 1000 to 1020 kg/m³ and stay flat. A
  flat density shows the box size has stabilized.

---

## Part 13: Production MD (100 ps)

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
```

```bash
nohup gmx mdrun -deffnm md -nt 2 > md.out 2>&1 &
```

This takes about 15 to 20 minutes. When it finishes, check the speed:

```bash
grep "Performance:" md.log
```

---

## Part 14: Analysis

### Step 1: remove periodic boundary effects

```bash
echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_noPBC.xtc -pbc mol -center
```

The two answers: center on `Protein`, write out the whole `System`.

### Step 2: protein stability (backbone RMSD)

**RMSD** (root mean square deviation) measures how far a structure has
moved from its starting position, averaged over the selected atoms.

```bash
echo "Backbone Backbone" | gmx rms -s md.tpr -f md_noPBC.xtc -n index.ndx -o rmsd_protein.xvg -tu ps
```

The two answers: **fit** to the backbone (remove overall rotation and
movement), then **measure** the RMSD of the backbone.

```bash
python plot_xvg.py rmsd_protein.xvg --unit A
```

### Step 3: ligand stability (ligand RMSD)

This is the most important check in a protein-ligand simulation: does
the ligand stay in its binding pose?

```bash
echo "Backbone JZ4" | gmx rms -s md.tpr -f md_noPBC.xtc -n index.ndx -o rmsd_ligand.xvg -tu ps
```

Here we **fit** to the protein backbone but **measure** the ligand. So
this RMSD shows how much the ligand moves **relative to the binding
pocket**.

Plot the protein and ligand RMSD together:

```bash
python plot_xvg.py rmsd_protein.xvg rmsd_ligand.xvg --unit A -o rmsd_compare.png --title "Protein and ligand RMSD"
```

**How to read it:**

| Ligand RMSD | Meaning |
|---|---|
| Below about 2 Å and flat | The ligand keeps its binding pose: a stable complex |
| Rises gradually then levels off | The ligand adjusts within the pocket, then settles |
| Rises steadily above 4–5 Å | The ligand is leaving its pose, or leaving the pocket |

### Step 4: does the ligand stay in contact with the protein?

`gmx mindist` calculates the shortest distance between any protein atom
and any ligand atom, at every time point:

```bash
echo "Protein JZ4" | gmx mindist -s md.tpr -f md_noPBC.xtc -n index.ndx -od mindist.xvg -tu ps
```

```bash
python plot_xvg.py mindist.xvg --unit A
```

A bound ligand stays in close contact, at around 2 to 3 Å. If the
distance grows steadily, the ligand is leaving the pocket.

### Step 5: hydrogen bonds between protein and ligand

```bash
echo "Protein JZ4" | gmx hbond-legacy -s md.tpr -f md_noPBC.xtc -n index.ndx -num hbnum.xvg -tu ps
```

```bash
python plot_xvg.py hbnum.xvg
```

The plot shows how many hydrogen bonds exist at each time point. JZ4
has one group that can form hydrogen bonds: its OH.

> **Why `hbond-legacy`?** GROMACS 2024 has a newer `gmx hbond` that
> identifies hydrogen-bond donors and acceptors from element information
> in the topology. Ligand topologies made by ACPYPE do not contain that
> information, so the new tool finds nothing on the ligand. The legacy
> tool identifies them from the atom names and works correctly here.

---

## Part 15: Watch your simulation

The full trajectory includes thousands of water molecules and is too
large to download comfortably. Make a smaller version with only the
protein and ligand, with the protein's rotation removed so it does
not tumble on screen:

```bash
echo "Backbone Protein_JZ4" | gmx trjconv -s md.tpr -f md_noPBC.xtc -n index.ndx -fit rot+trans -o complex_traj.xtc
```

The two answers: **fit** on the backbone, **write out** the protein and
ligand.

A viewer also needs a structure file with the same atoms, to know what
each atom is:

```bash
echo "Protein_JZ4" | gmx trjconv -s md.tpr -f complex_traj.xtc -n index.ndx -o complex_start.pdb -dump 0
```

Download **both** `complex_start.pdb` and `complex_traj.xtc`.

> **👁 View it: the trajectory in VMD**
> 1. Open VMD. Choose **File → New Molecule**, click **Browse**, select
>    `complex_start.pdb`, and click **Load**.
> 2. In the same window, make sure the new molecule is selected under
>    **Load files for**, then browse to `complex_traj.xtc` and click
>    **Load**.
> 3. Choose **Graphics → Representations**. Change **Drawing Method**
>    to **NewCartoon**.
> 4. Click **Create Rep**, type `resname JZ4` in **Selected Atoms**, and
>    set **Drawing Method** to **Licorice**.
> 5. Press **Play** (the ▶ button at the bottom right of the main VMD
>    window).
>
> **In PyMOL instead**, type:
> ```
> load complex_start.pdb
> load_traj complex_traj.xtc, complex_start
> show sticks, resn JZ4
> ```
> then press the **▶** button at the bottom right to play.
>
> **What to look for:** does the ligand stay in the pocket? Does the
> OH group point towards the same residues throughout? Compare what you
> see with your RMSD and hydrogen-bond plots.

---

## Changing the simulation length

As in Module 01, the length of each stage is `nsteps × dt` in its
`.mdp` file. With `dt = 0.002` ps:

| Simulated time | `nsteps` |
|---|---|
| 50 ps | 25000 |
| 100 ps | 50000 |
| 500 ps | 250000 |
| 1 ns | 500000 |

This system is a little larger than lysozyme in Module 01, so expect
each run to take about 20% longer. To change the production length,
edit `md.mdp`:

```bash
nano md.mdp
```

Then rerun the two commands in Part 13. To extend a finished run,
see "Extending a finished run" in Module 01.

---

## If a step fails

Copy only the reference files for the stage where you got stuck, then
continue. First, go to your folder:

```bash
cd /workspaces/md-course/my-work/02-complex
```

**Stuck in Parts 2 to 7 (preparation)?** Copy the prepared complex and
its topology, then continue from Part 8. Note that the reference
topology for this stage is called `topol_complex.top`, and you save it
as `topol.top`:

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/complex.gro .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/topol_complex.top topol.top
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/posre.itp .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/JZ4_GMX.itp .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/posre_JZ4.itp .
```

**Stuck in Parts 8 to 13?** First copy the full topology (with water
and ions) and the parameter files:

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/topol.top .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/posre.itp .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/JZ4_GMX.itp .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/posre_JZ4.itp .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/index.ndx .
```

Then copy the output of the last stage that worked, and continue from
the next part:

| If this part failed | Also copy | Then continue from |
|---|---|---|
| Part 8 (box, water, ions) or Part 9 (minimization) | `em.gro` | Part 11 |
| Part 11 (NVT) | `nvt.gro` and `nvt.cpt` | Part 12 |
| Part 12 (NPT) | `npt.gro` and `npt.cpt` | Part 13 |
| Part 13 (production) | `md.tpr`, `md.xtc` and `md.edr` | Part 14 |

For example, if your NPT run (Part 12) failed, copy the reference NPT
output and continue with production in Part 13:

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/npt.gro .
```

```bash
cp /workspaces/md-course/modules/02-protein-ligand/reference/npt.cpt .
```

Then ask in the help channel about what went wrong, so you learn from it.

### Common errors in protein-ligand setup

| Error message contains | Usual cause and fix |
|---|---|
| `Atomtype ca not found` (or another lowercase type) | `JZ4_GMX.itp` is not included, or is included in the wrong place. It must come directly after the force field `#include` line (Part 6). |
| `Invalid order for directive atomtypes` | `JZ4_GMX.itp` was included after a `[ moleculetype ]`. Move it up, directly after the force field line. |
| `No such moleculetype JZ4` | The ligand `#include` line is missing, or the name in `[ molecules ]` is misspelled. |
| `number of coordinates in coordinate file does not match topology` | The atom count on line 2 of `complex.gro` is wrong, or the order in `[ molecules ]` does not match the coordinates (Protein first, then JZ4). |
| `Group Protein_JZ4 referenced in the .mdp file was not found` | You forgot `-n index.ndx` in the `grompp` command, or skipped Part 10. |
| ACPYPE stops with a charge error | The net charge given with `-n` does not match the molecule. Check the formula and protonation state (Part 3). |

---

## Exercises

1. What other non-protein molecules are in 3HTB, apart from water and
   JZ4? Why did we remove them?
2. How many hydrogen atoms did Open Babel add to JZ4?
3. Using `JZ4_GMX.itp`, list the GAFF2 atom types in JZ4 and what each
   one means. Which atom has the most negative partial charge?
4. What is the net charge of the protein? How many Na⁺ and Cl⁻ ions were
   added?
5. Why must `#include "JZ4_GMX.itp"` come directly after the force field
   line in `topol.top`?
6. What were the average temperature (NVT), pressure and density (NPT)?
   Include your plots.
7. What is the average ligand RMSD during production? Is the binding
   pose stable? Use your `rmsd_compare.png` to explain.
8. What is the average number of protein-ligand hydrogen bonds? Using
   PyMOL or VMD, find which residue the JZ4 oxygen interacts with.
9. **Challenge:** in Part 14, Step 3 we fitted on the protein backbone
   but measured the ligand. Run the command again with
   `echo "JZ4 JZ4"`, which fits **and** measures on the ligand.
   Plot both results together. What is the difference, and what does
   each one tell you?

---

## Files you should have at the end

```
3HTB.pdb  protein.pdb  JZ4.pdb  JZ4_H.mol2  JZ4.acpype/
JZ4_GMX.itp  JZ4_GMX.gro  posre_JZ4.itp
protein.gro  topol.top  posre.itp  complex.gro
newbox.gro  solv.gro  solv_ions.gro  index.ndx
em.*  nvt.*  npt.*  md.*  md_noPBC.xtc
rmsd_protein.xvg  rmsd_ligand.xvg  mindist.xvg  hbnum.xvg
potential.png  temperature.png  pressure.png  density.png
rmsd_compare.png  mindist.png  hbnum.png
complex_start.pdb  complex_traj.xtc
```

Remember to stop your codespace when you finish.
