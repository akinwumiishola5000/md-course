# Module 01: Simulating a Protein in Water

In this module you run a complete molecular dynamics simulation of
**hen egg-white lysozyme (PDB: 1AKI)** in water, from the crystal
structure to a 100 ps production run. Every simulation you run in your
research follows these same steps, and Module 02 builds directly on
them.

**Time needed:** about 3 hours, including about 35 minutes of computing
time that runs in the background
**Before you start:** complete Module 00

---

## Before the session: install a molecular viewer on your computer

Your codespace has no screen for 3D graphics, so you will **download**
structures from the codespace and look at them on your own computer.
Install at least one of these before the session. All are free for
academic use.

| Program | Best for | Get it from |
|---|---|---|
| **PyMOL** | Looking at structures, making figures | pymol.org (free educational licenses are available; an open-source version also exists) |
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
1AKI.pdb  →  clean structure  →  topology  →  box  →  water  →  ions
                                                                  ↓
   production MD  ←  NPT equilibration  ←  NVT equilibration  ←  energy minimization
```

| Part | What you do | GROMACS tool | Main output |
|---|---|---|---|
| 1 | Download and clean the structure | `grep` | `1AKI_clean.pdb` |
| 2 | Build the topology | `gmx pdb2gmx` | `processed.gro`, `topol.top`, `posre.itp` |
| 3 | Define the box | `gmx editconf` | `newbox.gro` |
| 4 | Add water | `gmx solvate` | `solv.gro` |
| 5 | Add ions | `gmx grompp` + `gmx genion` | `solv_ions.gro` |
| 6 | Energy minimization | `gmx grompp` + `gmx mdrun` | `em.gro` |
| 7 | NVT equilibration (50 ps) | `gmx grompp` + `gmx mdrun` | `nvt.gro`, `nvt.cpt` |
| 8 | NPT equilibration (50 ps) | `gmx grompp` + `gmx mdrun` | `npt.gro`, `npt.cpt` |
| 9 | Production MD (100 ps) | `gmx grompp` + `gmx mdrun` | `md.xtc` |
| 10 | First look at the results | `gmx rms`, `gmx gyrate` | plots and a movie |

GROMACS always works in two stages:

1. **`gmx grompp`** (GROMACS pre-processor) combines the **structure**
   (`.gro`), the **topology** (`.top`) and the **settings** (`.mdp`) into
   one **run input file** (`.tpr`), and checks them for errors.
2. **`gmx mdrun`** runs the calculation described in the `.tpr` file.

---

## The files you will meet

You will create many files in this module. Come back to this table
whenever you are unsure what a file is.

| Extension | What it is | Can you read it? |
|---|---|---|
| `.pdb` | Protein Data Bank structure: atoms and coordinates in **Å** | Yes, text |
| `.gro` | GROMACS structure: atoms and coordinates in **nm**, plus the box size | Yes, text |
| `.top` | The **topology** of the whole system: which molecules, how many, and their force field parameters | Yes, text |
| `.itp` | An **included topology**: parameters or restraints for one molecule, used inside `.top` | Yes, text |
| `.mdp` | Simulation **settings** (time step, temperature, length...) | Yes, text |
| `.tpr` | **Run input**: structure + topology + settings combined by `grompp` | No, binary |
| `.cpt` | **Checkpoint**: positions and velocities, to continue or restart a run | No, binary |
| `.edr` | **Energy file**: energies, temperature, pressure over time | No, use `gmx energy` |
| `.xtc` | **Trajectory**: compressed coordinates over time (the "movie") | No, use a viewer |
| `.log` | Log of a run, including performance | Yes, text |
| `.xvg` | Data for plotting, written by GROMACS analysis tools | Yes, text |

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
mkdir -p my-work/01-lysozyme
```

```bash
cd my-work/01-lysozyme
```

Copy the simulation settings and the plotting tool:

```bash
cp /workspaces/md-course/modules/01-protein-in-water/files/*.mdp .
```

```bash
cp /workspaces/md-course/modules/tools/plot_xvg.py .
```

```bash
ls
```

You should see five settings files (`ions.mdp`, `em.mdp`, `nvt.mdp`,
`npt.mdp`, `md.mdp`) and `plot_xvg.py`.

### About the plotting tool

GROMACS writes analysis results as `.xvg` files: plain columns of
numbers. `plot_xvg.py` turns them into pictures, and prints the mean,
standard deviation, minimum and maximum of each line:

```bash
python plot_xvg.py file.xvg
```

It saves `file.png` in the same folder. **Click the `.png` file in the
file panel to open it.** You will use these options in this module:

| Option | Effect |
|---|---|
| `--unit A` | convert nm to Å (most papers report distances in Å) |
| `--avg 10` | add a running average over 10 points, to see the trend through the noise |
| `--cols 1` | plot only the first data column |
| `--xlabel "..."`, `--ylabel "..."`, `--title "..."` | change the labels |

---

## Part 1: Download and clean the structure

```bash
wget https://files.rcsb.org/download/1AKI.pdb
```

Look at what the file contains:

```bash
grep "^TITLE" 1AKI.pdb
```

```bash
grep "REMARK   2 RESOLUTION" 1AKI.pdb
```

```bash
grep "^HETATM" 1AKI.pdb | cut -c18-20 | sort | uniq -c
```

The only non-protein molecules are crystal waters (`HOH`). Remove them,
because in Part 4 you will surround the protein with a complete box of
water instead:

```bash
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
```

Check that no water is left:

```bash
grep -c HOH 1AKI_clean.pdb
```

It should print `0`.

> **In your own research:** always inspect a structure before using it.
> Look for missing residues (`grep "REMARK 465"`), missing atoms
> (`grep "REMARK 470"`), ligands and ions you may need to keep, and
> alternate side-chain positions. 1AKI is a clean, complete structure,
> which is why it is good for learning.

> **👁 View it: the crystal structure**
> Download `1AKI.pdb` and open it in **PyMOL**. Type these commands in
> the PyMOL command line (the `PyMOL>` box), one at a time:
> ```
> hide everything
> show cartoon
> show spheres, resn HOH
> set sphere_scale, 0.3
> show sticks, resn CYS
> ```
> The red dots are the crystal waters you just removed. The sticks are
> the cysteines: lysozyme has four **disulfide bonds** that hold its
> fold together. Notice there are **no hydrogen atoms**: X-ray
> crystallography usually cannot see them.
>
> In **Mol\*** (no installation): open molstar.org/viewer and drag
> `1AKI.pdb` into the window.

---

## Part 2: Build the topology

The **topology** describes every atom in the system: its mass, its
charge, which atoms it is bonded to, and the force field parameters for
every bond, angle and dihedral.

```bash
gmx pdb2gmx -f 1AKI_clean.pdb -o processed.gro -water tip3p -ff amber99sb-ildn -ignh
```

| Option | Meaning |
|---|---|
| `-f 1AKI_clean.pdb` | input structure |
| `-o processed.gro` | output structure, in GROMACS format |
| `-ff amber99sb-ildn` | the **force field**: the equations and parameters describing the forces between atoms |
| `-water tip3p` | the **water model** used with this force field |
| `-ignh` | ignore any hydrogens in the file and add them back with the force field's names |

We use **AMBER99SB-ILDN** with **TIP3P** water because it is a widely
used protein force field, and it is compatible with the GAFF2 ligand
parameters you will use in Module 02.

`pdb2gmx` created three files:

| File | Contains |
|---|---|
| `processed.gro` | the protein, now **with hydrogen atoms** |
| `topol.top` | the topology |
| `posre.itp` | **position restraints**: used during equilibration to hold the protein's heavy atoms in place |

### Understanding the .gro file

```bash
head -n 4 processed.gro
```

```bash
tail -n 1 processed.gro
```

| Line | Contains |
|---|---|
| Line 1 | A title |
| Line 2 | The **number of atoms** |
| One line per atom | residue number, residue name, atom name, atom number, then **x, y, z in nm** |
| Last line | The **box size** in nm |

Compare the number of atoms with the original PDB file. How many
hydrogens did `pdb2gmx` add?

```bash
grep -c "^ATOM" 1AKI_clean.pdb
```

```bash
sed -n 2p processed.gro
```

### Understanding the topology

The topology is long, because it lists every atom and every bond. Look
at its structure instead of reading it all. First, the files it
**includes**:

```bash
grep -n "#include" topol.top
```

| Included file | What it adds |
|---|---|
| `amber99sb-ildn.ff/forcefield.itp` | the force field: every **atom type** and its parameters |
| `posre.itp` | the protein position restraints, used only when an `.mdp` file says `define = -DPOSRES` |
| `amber99sb-ildn.ff/tip3p.itp` | the water model |
| `amber99sb-ildn.ff/ions.itp` | ion parameters (Na⁺, Cl⁻ and others) |

Now the sections inside the topology, which each start with a name in
square brackets:

```bash
grep -n "^\[" topol.top | head -n 12
```

| Section | What it describes |
|---|---|
| `[ moleculetype ]` | the molecule's name. For 1AKI it is `Protein_chain_A`, because the protein is chain A in the PDB file |
| `[ atoms ]` | one line per atom: type, residue, name, **partial charge** and mass |
| `[ bonds ]`, `[ pairs ]`, `[ angles ]`, `[ dihedrals ]` | which atoms are connected, and the force field terms for each connection |
| `[ system ]` | a name for the whole system |
| `[ molecules ]` | every molecule in the system and how many copies, **in the same order as the coordinate file** |

Look at the end of the file:

```bash
tail -n 3 topol.top
```

Right now, the system contains one protein molecule. `gmx solvate` and
`gmx genion` will add lines here for water and ions.

### The protein's charge

The last atom in the `[ atoms ]` section has a running total of the
charge, after `qtot`:

```bash
grep "qtot" topol.top | tail -n 1
```

The number at the end is the protein's **net charge**. Lysozyme is
**+8**: it has more positive (lysine, arginine) than negative
(aspartate, glutamate) residues. Write this down; you neutralize it in
Part 5.

### Position restraints

```bash
head -n 8 posre.itp
```

Each line restrains one heavy atom with a force constant of 1000
kJ/mol/nm² in x, y and z. Hydrogens are not restrained.

> **👁 View it: the protein with hydrogens**
> Download `processed.gro` and open it in PyMOL (PyMOL reads `.gro`
> files). Type:
> ```
> show sticks
> zoom resi 35
> ```
> Every residue now has hydrogens. Residue 35 is glutamate 35, one of
> the two catalytic residues of lysozyme. Compare with `1AKI.pdb`, which
> had none.

---

## Part 3: Define the simulation box

Simulations use **periodic boundary conditions**: the box is surrounded
by copies of itself in every direction, so the protein never meets an
edge or a vacuum.

```bash
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt dodecahedron
```

| Option | Meaning |
|---|---|
| `-c` | centre the protein in the box |
| `-d 1.0` | keep at least 1.0 nm between the protein and the box edge |
| `-bt dodecahedron` | a **rhombic dodecahedron** box, which needs about 30% less water than a cube for the same protein, so the simulation runs faster |

Look at the new box size on the last line:

```bash
tail -n 1 newbox.gro
```

> **Why 1.0 nm?** The protein must never "see" its own copy through
> the periodic boundary. With 1.0 nm on each side, copies are at least
> 2.0 nm apart: more than the 1.0 nm interaction cutoff used in the
> `.mdp` files.

---

## Part 4: Add water

```bash
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
```

| Option | Meaning |
|---|---|
| `-cp newbox.gro` | the protein in its box |
| `-cs spc216.gro` | a small box of pre-equilibrated water that GROMACS copies to fill the space (it works for any three-site water model, including TIP3P) |
| `-p topol.top` | add the water molecules to the topology |

How many water molecules were added?

```bash
tail -n 3 topol.top
```

`SOL` is the name of the water molecules. Expect about 7,300.

---

## Part 5: Add ions

The system has a net charge of +8, which must be neutralized. We also
add salt at a physiological concentration of 0.15 M.

`gmx genion` replaces some water molecules with ions. It needs a `.tpr`
file, so first build one with `grompp`:

```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
```

> **Why `-maxwarn 1`?** `grompp` warns that the system has a net charge.
> That is true and expected: you are about to fix it. This is the
> **only** situation in the course where you override a warning. In
> general, a GROMACS warning means something is wrong and must be
> fixed, not ignored.

Now add the ions:

```bash
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```

| Part | Meaning |
|---|---|
| `echo "SOL" \|` | answers genion's question about which molecules to replace: the water |
| `-pname NA -nname CL` | the names of the positive and negative ions |
| `-neutral` | add enough ions to make the total charge zero |
| `-conc 0.15` | add 0.15 M NaCl on top of that |

Check the result:

```bash
tail -n 4 topol.top
```

Expect about **22 NA and 30 CL**: 8 more chloride than sodium ions, to
cancel the protein's +8 charge.

How many atoms does the whole system have?

```bash
sed -n 2p solv_ions.gro
```

About 24,000 atoms, of which only about 2,000 are protein.

> **👁 View it: the solvated system**
> Download `solv_ions.gro` (about 1 MB) and open it in **PyMOL**. Type:
> ```
> show spheres, resn NA+CL
> show cartoon, polymer
> hide lines, polymer
> ```
> The lines are water molecules. Notice the shape of the dodecahedral box,
> and the ions scattered through the water.

---

## Part 6: Energy minimization

Adding water and ions can place some atoms too close together. Energy
minimization moves atoms slightly to remove these bad contacts, so the
simulation does not become unstable when it starts.

Look at the settings first:

```bash
cat em.mdp
```

| Setting | Meaning |
|---|---|
| `integrator = steep` | **steepest descent**: move each atom downhill in energy |
| `emtol = 1000.0` | stop when the largest force on any atom is below 1000 kJ/mol/nm |
| `nsteps = 50000` | stop anyway after this many steps |
| `define = -DFLEXIBLE` | let water molecules bend during minimization, which avoids instabilities |

Build the run input and run it:

```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
```

```bash
gmx mdrun -deffnm em -nt 2
```

`-deffnm em` uses `em` as the base name for all files: it reads
`em.tpr` and writes `em.gro`, `em.log`, `em.edr` and `em.trr`. `-nt 2`
uses 2 CPU threads, matching the 2 cores of your codespace.

Check the result:

```bash
grep -E "converged|Potential Energy|Maximum force" em.log
```

You should see that steepest descents **converged**, a large negative
**potential energy** (hundreds of thousands of kJ/mol), and a **maximum
force** below 1000.

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

## Part 7: NVT equilibration (50 ps)

Now the system is heated to 300 K at constant **N**umber of atoms,
**V**olume and **T**emperature (**NVT**). The protein is held in place
with position restraints, so the water and ions can arrange themselves
around it without the protein structure changing.

### The settings

```bash
cat nvt.mdp
```

The most important settings, which you will meet in every simulation:

| Setting | Value | Meaning |
|---|---|---|
| `define` | `-DPOSRES` | switch on the position restraints in `posre.itp` |
| `integrator` | `md` | molecular dynamics (the leap-frog algorithm) |
| `dt` | `0.002` | time step: 0.002 ps = 2 fs |
| `nsteps` | `25000` | number of steps: 25000 × 0.002 ps = **50 ps** |
| `nstxout-compressed` | `500` | save coordinates every 500 steps (1 ps) to the trajectory |
| `constraints` | `h-bonds` | keep bonds to hydrogen at fixed length, which allows the 2 fs time step |
| `coulombtype` | `PME` | calculate long-range electrostatics with Particle Mesh Ewald |
| `rcoulomb`, `rvdw` | `1.0` | short-range cutoffs, in nm |
| `tcoupl` | `V-rescale` | the **thermostat** that keeps the temperature at the target |
| `tc-grps` | `Protein Non-Protein` | the protein and the solvent get separate thermostats |
| `ref_t` | `300 300` | the target temperature for each group, in K |
| `pcoupl` | `no` | no pressure control: the volume stays fixed (that is the "V" in NVT) |
| `gen_vel` | `yes` | give the atoms random starting velocities for 300 K |

### Run it

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```

`-c em.gro` is the starting structure, and `-r em.gro` gives the
**reference** positions for the restraints.

```bash
nohup gmx mdrun -deffnm nvt -nt 2 > nvt.out 2>&1 &
```

As in Module 00, `nohup ... &` runs the simulation in the background,
so it keeps going if your browser disconnects. This takes about 5 to 10
minutes. Check on it with:

```bash
tail -n 15 nvt.log
```

To see whether it is still running:

```bash
jobs
```

`jobs` only knows about jobs started **in the same terminal**. If it
prints nothing but you are not sure, check the end of the log: a
finished run ends with `Finished mdrun`.

### Check the temperature

When it has finished, extract the temperature:

```bash
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
```

```bash
python plot_xvg.py temperature.xvg --avg 10
```

**What you should see:** the temperature rises quickly to 300 K and then
fluctuates around it. The thick line (the running average) should be
flat at about 300 K. Fluctuations of a few kelvin are normal: a
system of 24,000 atoms is too small for a perfectly steady temperature.

---

## Part 8: NPT equilibration (50 ps)

Next, the pressure is brought to 1 bar at constant **N**umber,
**P**ressure and **T**emperature (**NPT**). The box size adjusts until
the water reaches the right density. The protein is still restrained.

Compare the settings with NVT:

```bash
grep -E "continuation|pcoupl|tau_p|ref_p|compressibility|gen_vel" npt.mdp
```

| Setting | Meaning |
|---|---|
| `continuation = yes` | continue from NVT instead of starting fresh |
| `pcoupl = C-rescale` | the **barostat** that keeps the pressure at the target |
| `ref_p = 1.0` | the target pressure, in bar |
| `compressibility = 4.5e-5` | how much water's volume changes with pressure |
| `gen_vel = no` | keep the velocities from NVT |

### Run it

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```

`-t nvt.cpt` passes the **checkpoint** from NVT, so the simulation
continues with the same atom velocities instead of starting again.

```bash
nohup gmx mdrun -deffnm npt -nt 2 > npt.out 2>&1 &
```

### Check pressure and density

When it has finished, extract the pressure and the density **separately**,
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

- **Pressure** jumps up and down by hundreds of bar. That is normal for a
  system this small. Look at the **mean** printed in the terminal and
  the running average: they should be near 1 bar, within the large
  fluctuations.
- **Density** should settle around 1000 to 1020 kg/m³ and stay flat. A
  flat density shows the box size has stabilized.

---

## Part 9: Production MD (100 ps)

The system is now equilibrated. The production run has **no
restraints**: the protein moves freely. This is the trajectory you
analyse.

Check what changed from NPT:

```bash
grep -E "define|nsteps" md.mdp
```

There is no `define = -DPOSRES` line, so the restraints are off, and
`nsteps = 50000` gives 100 ps.

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
```

Notice there is no `-r`: without restraints, no reference positions are
needed.

```bash
nohup gmx mdrun -deffnm md -nt 2 > md.out 2>&1 &
```

This takes about 10 to 20 minutes. When it has finished, check the
speed:

```bash
grep "Performance:" md.log
```

The first number is **ns/day**: how many nanoseconds of simulation
your codespace produces per day for this system. Write it down; you use
it below to plan longer simulations.

---

## Part 10: First look at the results

Module 03 teaches trajectory analysis in depth. Here are three quick
checks that every simulation should pass.

### Step 1: fix the periodic boundaries

Because of periodic boundaries, the protein can appear broken across the
edges of the box. Make it whole and centred:

```bash
echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center
```

The two answers: centre on `Protein`, and write out the whole `System`.

### Step 2: RMSD — how much has the structure changed?

**RMSD** (root mean square deviation) is the average distance between
the atoms of each frame and the same atoms in the starting structure,
after superimposing them.

```bash
echo "Backbone Backbone" | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ps
```

The two answers: **fit** on the backbone (remove overall rotation and
movement), then **measure** the RMSD of the backbone.

```bash
python plot_xvg.py rmsd.xvg --unit A --avg 10
```

**What you should see:** the RMSD rises quickly from 0 in the first few
ps, as the structure relaxes from the crystal structure into solution,
then stays roughly level, below about 1.5 Å. A stable, well-folded
protein like lysozyme does not change much in 100 ps.

### Step 3: radius of gyration — is the protein staying compact?

The **radius of gyration (Rg)** measures how spread out the protein's
atoms are around its centre. A folded protein keeps an almost constant
Rg; unfolding makes it increase.

```bash
echo "Protein" | gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg
```

`gyrate.xvg` has four columns of data: the total Rg, and its
components around the x, y and z axes. Plot only the total, the first
column:

```bash
python plot_xvg.py gyrate.xvg --unit A --cols 1
```

**What you should see:** a flat line at about 14 Å, changing by only a
fraction of an Å over 100 ps.

### Step 4: watch the simulation

The full trajectory includes about 22,000 water atoms, which makes it
large and hard to see. Make a smaller version with only the protein,
with its rotation and movement removed so it does not tumble on screen:

```bash
echo "Backbone Protein" | gmx trjconv -s md.tpr -f md_noPBC.xtc -o protein_traj.xtc -fit rot+trans
```

The two answers: **fit** on the backbone, **write out** the protein.

A viewer also needs a structure file with the same atoms, to know what
each atom is. Save the first frame as a PDB file:

```bash
echo "Protein" | gmx trjconv -s md.tpr -f protein_traj.xtc -o protein_start.pdb -dump 0
```

Download **both** `protein_start.pdb` and `protein_traj.xtc`.

> **👁 View it: your simulation in VMD**
> 1. Open VMD. Choose **File → New Molecule**, click **Browse**, select
>    `protein_start.pdb`, and click **Load**.
> 2. In the same window, make sure the new molecule is selected under
>    **Load files for**, then browse to `protein_traj.xtc` and click
>    **Load**.
> 3. Choose **Graphics → Representations**. Change **Drawing Method**
>    to **NewCartoon**.
> 4. Press **Play** (the ▶ button at the bottom right of the main VMD
>    window).
>
> **In PyMOL instead**, type:
> ```
> load protein_start.pdb
> load_traj protein_traj.xtc, protein_start
> hide everything
> show cartoon
> ```
> then press the **▶** button at the bottom right to play.
>
> **What to look for:** the helices and sheets keep their shape, while
> the loops and the two ends of the chain move more. You will measure
> this flexibility properly (RMSF) in Module 03.

---

## Changing the simulation length

The length of every stage is set by two lines in its `.mdp` file:

```
dt      = 0.002    ; time step in ps (2 fs)
nsteps  = 50000    ; number of steps
```

**Simulated time = nsteps × dt.** Keep `dt = 0.002` and change only
`nsteps`:

| Simulated time | `nsteps` (with dt = 0.002) | Approximate run time on 2 cores |
|---|---|---|
| 50 ps | 25000 | 5–10 min |
| 100 ps | 50000 | 10–20 min |
| 500 ps | 250000 | 1–1.5 hours |
| 1 ns | 500000 | 2–3 hours |
| 10 ns | 5000000 | 1–2 days |

Run time depends on your machine and the size of your system. Estimate
it from the ns/day value in Part 9:

```
hours needed = (simulated time in ns) ÷ (ns/day) × 24
```

To change the production run to 500 ps, for example, open `md.mdp`:

```bash
nano md.mdp
```

Change `nsteps` to `250000`, save with **Ctrl + O** and **Enter**, and exit
with **Ctrl + X**. Check the change:

```bash
grep nsteps md.mdp
```

Then run the `grompp` and `mdrun` commands of Part 9 again.

> **Why only 100 ps here?** This course teaches the workflow, not
> research-quality sampling. 100 ps is long enough to see the protein
> move and to practise analysis. Real studies use tens to hundreds of
> nanoseconds, usually on GPUs or HPC clusters (Module 04).

### Extending a finished run

To add 400 ps to a finished 100 ps production run (500 ps in total)
without starting again, first make a new run input:

```bash
gmx convert-tpr -s md.tpr -extend 400 -o md_extended.tpr
```

`-extend` is in ps. Then continue from the checkpoint:

```bash
nohup gmx mdrun -s md_extended.tpr -cpi md.cpt -deffnm md -nt 2 > md_ext.out 2>&1 &
```

New frames are added to the end of `md.xtc`.

### If your codespace stops during a run

If your codespace stops (idle timeout, closed browser, lost connection),
the running simulation stops too. Your files are kept, and GROMACS
writes a checkpoint (`.cpt`) every 15 minutes. Reopen the codespace, go
back to your folder:

```bash
cd /workspaces/md-course/my-work/01-lysozyme
```

and continue from the checkpoint:

```bash
nohup gmx mdrun -deffnm md -cpi md.cpt -nt 2 > md_restart.out 2>&1 &
```

Replace `md` with `nvt` or `npt` (in both places) if that was the stage
that stopped.

If it stopped during an **extended** run, tell GROMACS to use the
extended run input, otherwise it will say the simulation is already
complete:

```bash
nohup gmx mdrun -s md_extended.tpr -cpi md.cpt -deffnm md -nt 2 > md_restart.out 2>&1 &
``` To avoid stops, keep the codespace tab open during long
runs, and set your idle timeout to 240 minutes at
github.com/settings/codespaces.

---

## If a step fails

Each step needs the files from the step before. If a step fails and you
cannot fix it, copy the reference files for that stage and continue.
Copy only the files listed for the stage where you got stuck, so you do
not overwrite your own work.

First, go to your folder:

```bash
cd /workspaces/md-course/my-work/01-lysozyme
```

Set a shortcut to the reference folder, so the copy commands are short:

```bash
R=/workspaces/md-course/modules/01-protein-in-water/reference
```

(The shortcut lasts until you close the terminal. In a new terminal, run
the line again.) Then copy the files for your stage, from the table:

| If this part failed | Copy these files | Then continue from |
|---|---|---|
| Parts 1 to 5 (building the system) | `cp $R/{solv_ions.gro,topol.top,posre.itp} .` | Part 6 |
| Part 6 (energy minimization) | `cp $R/{em.gro,topol.top,posre.itp} .` | Part 7 |
| Part 7 (NVT) | `cp $R/{nvt.gro,nvt.cpt,topol.top,posre.itp} .` | Part 8 |
| Part 8 (NPT) | `cp $R/{npt.gro,npt.cpt,topol.top,posre.itp} .` | Part 9 |
| Part 9 (production) | `cp $R/{md.tpr,md.xtc,md.edr,md.gro} .` | Part 10 |

For example, if your NPT run failed:

```bash
cp $R/{npt.gro,npt.cpt,topol.top,posre.itp} .
```

Then ask in the help channel about what went wrong, so you learn from it.

### Common errors

| Error message contains | Usual cause and fix |
|---|---|
| `No such file or directory` | You are in the wrong folder (check with `pwd`), or a previous step did not produce its output (check with `ls`). |
| `Residue ... not found in residue topology database` | The structure contains a molecule the force field does not know (a ligand or a non-standard residue). Remove it, or parameterize it as in Module 02. |
| `Atom ... in residue ... not found in rtp entry` | Atom names do not match the force field; usually fixed with `-ignh`. |
| `There were N warnings` / `Too many warnings` | Read the warning. Do not add `-maxwarn` to hide it. |
| `number of coordinates in coordinate file does not match topology` | `topol.top` and the `.gro` file are out of step, for example because a step was run twice. Start again from Part 2, or copy the reference files. |
| `LINCS warnings` or `water molecules can not be settled` | The system is unstable. Check that energy minimization converged (Part 6). |
| `The PNG did not appear` | Click the refresh icon at the top of the file panel. |

---

## Exercises

Write down your answers and keep your plots: you will compare them with
the protein-ligand complex in Module 02.

1. What is the resolution of 1AKI? How many hydrogen atoms did
   `pdb2gmx` add?
2. What is the net charge of the protein? How many Na⁺ and Cl⁻ ions were
   added, and why are the two numbers different?
3. How many water molecules are in your system, and how many atoms in
   total?
4. How many steps did energy minimization take, and what was the final
   potential energy? Describe the shape of `potential.png`.
5. What were the average temperature (NVT), pressure and density (NPT)?
   Use the means printed by `plot_xvg.py`.
6. What was the performance of your production run in ns/day? Using the
   formula above, how long would a 100 ns simulation take on your
   codespace?
7. What is the average backbone RMSD over the production run, in Å? Does
   it level off?
8. Look at `nvt.mdp` and `md.mdp`. Which line makes the protein
   restrained during equilibration but free in production? Why do we
   restrain it during equilibration?
9. **Challenge:** extend your production run to 200 ps using
   `gmx convert-tpr`. Repeat Step 1 and Step 2 of Part 10. Does the RMSD
   keep increasing, or level off?

---

## Files you should have at the end

```
1AKI.pdb  1AKI_clean.pdb  processed.gro  topol.top  posre.itp
newbox.gro  solv.gro  solv_ions.gro  ions.tpr
em.tpr  em.gro  em.edr  em.log
nvt.tpr  nvt.gro  nvt.cpt  nvt.edr  nvt.log
npt.tpr  npt.gro  npt.cpt  npt.edr  npt.log
md.tpr  md.gro  md.cpt  md.edr  md.log  md.xtc
md_noPBC.xtc  protein_traj.xtc  protein_start.pdb
potential.png  temperature.png  pressure.png  density.png
rmsd.png  gyrate.png
```

Keep this folder. Module 02 uses the same steps for a protein-ligand
complex, and you will compare the two systems.

Remember to stop your codespace when you finish.
