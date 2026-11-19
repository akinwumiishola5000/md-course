# Module 04: Going Further

You can now set up, run and analyse a protein-ligand simulation. This
final module takes you towards research:

- **Part A (hands-on):** estimate how strongly JZ4 binds to T4 lysozyme
  with **MM/GBSA**, and find out which residues contribute most to
  binding.
- **Part B (reading):** how real studies run simulations that are
  thousands of times longer, on GPUs and HPC clusters.
- **Part C (reading):** a checklist for starting your own
  protein-ligand MD project.

**Time needed:** about 2.5 hours
**Before you start:** complete Modules 02 and 03

---

# Part A: Binding free energy with MM/GBSA

## A1. What does MM/GBSA calculate?

The **binding free energy (ΔG<sub>bind</sub>)** measures how strongly a
ligand binds: the more negative, the stronger the binding.

```
ΔG bind  =  G complex  −  G protein  −  G ligand
```

**MM/GBSA** (Molecular Mechanics / Generalized Born Surface Area)
estimates each **G** for many frames of your trajectory, and averages
the result. For each frame it:

1. takes the complex, and makes two copies: the protein alone and the
   ligand alone,
2. calculates the energy of all three, **removing the explicit water**
   and replacing it with a **continuum** model of water,
3. subtracts: complex − protein − ligand.

The result is split into terms that tell you **why** the ligand binds:

| Term | What it measures | Usually |
|---|---|---|
| **ΔVDWAALS** | van der Waals: close packing, hydrophobic contacts | negative (favours binding) |
| **ΔEEL** | electrostatics: charges, hydrogen bonds | negative |
| **ΔEGB** | **polar solvation**: the cost of stripping water away from polar groups when the ligand enters the pocket | **positive** (opposes binding) |
| **ΔESURF** | **non-polar solvation**: the gain from burying surface area away from water | small, negative |
| **ΔTOTAL** | the sum: the MM/GBSA binding free energy | negative for a binder |

ΔVDWAALS + ΔEEL is called **ΔGGAS** (the interaction in vacuum), and
ΔEGB + ΔESURF is **ΔGSOLV** (the effect of water).

Compare this with the interaction energy in Module 03, Part 9: MM/GBSA
adds the **solvation** terms, which is why it is a better estimate of
binding.

> **⚠ What MM/GBSA leaves out.** It ignores the **entropy** lost when
> the ligand and protein become less flexible on binding, and it
> treats water as a continuum. As a result, MM/GBSA values are usually
> **much more negative than experimental binding free energies**. They
> are best used to **compare similar ligands in the same protein**
> (which binds more strongly?), and to find **key residues**. Do not
> compare them directly with experimental ΔG values.

> **Units:** gmx_MMPBSA reports **kcal/mol**, while GROMACS uses
> kJ/mol. 1 kcal/mol = 4.184 kJ/mol.

## A2. Get ready

```bash
cd /workspaces/md-course
```

```bash
git pull
```

```bash
mkdir -p my-work/04-mmgbsa
```

```bash
cd my-work/04-mmgbsa
```

Copy the input file and the plotting script:

```bash
cp /workspaces/md-course/modules/04-going-further/files/mmgbsa.in .
```

```bash
cp /workspaces/md-course/modules/tools/plot_mmgbsa.py .
```

Copy your results from Module 02:

```bash
cp ../02-complex/md.tpr .
```

```bash
cp ../02-complex/md_noPBC.xtc .
```

```bash
cp ../02-complex/index.ndx .
```

```bash
cp ../02-complex/topol.top .
```

```bash
cp ../02-complex/posre.itp .
```

```bash
cp ../02-complex/JZ4_GMX.itp .
```

```bash
cp ../02-complex/posre_JZ4.itp .
```

```bash
ls
```

> **Why `md_noPBC.xtc`?** MM/GBSA needs every molecule **whole** in every
> frame. A ligand split across the periodic boundary would give
> meaningless energies. `md_noPBC.xtc` was made whole in Module 02.

> **Did not finish Module 02?** Replace `../02-complex/` with
> `/workspaces/md-course/modules/02-protein-ligand/reference/` in the
> commands above. The reference folder has `md.xtc` instead of
> `md_noPBC.xtc`; make it whole first:
> ```bash
> echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_noPBC.xtc -pbc mol -center
> ```

## A3. Find the group numbers

gmx_MMPBSA needs to know which atoms are the **receptor** (the protein)
and which are the **ligand**. It identifies them by their **group
number** in `index.ndx`. List the groups with their numbers:

```bash
awk '/^\[/ {print n++, $2}' index.ndx
```

This prints each group name (the text between `[ ]`) with a number
counting from 0. Find the numbers for **`Protein`** and **`JZ4`**. In the
course system they are **1** and **13**. If yours are different, use
your numbers in A5.

## A4. Understand the input file

```bash
cat mmgbsa.in
```

| Setting | Meaning |
|---|---|
| `interval = 5` | analyse every 5th frame: 21 of your 101 frames. Consecutive frames 1 ps apart are very similar, so skipping some saves time with little loss |
| `PBRadii = 3` | the set of atomic radii used by the GB model (mbondi2, the standard choice with `igb = 5`) |
| `igb = 5` | which **Generalized Born** model of water to use (GB-OBC2, widely used for proteins) |
| `saltcon = 0.150` | salt concentration (M), the same as in your simulation |
| `idecomp = 2` | also split the energy **per residue** |
| `print_res = "within 4"` | report only residues within 4 Å of the ligand |

Lines starting with `#` are comments, which gmx_MMPBSA ignores.

## A5. Run MM/GBSA

```bash
nohup gmx_MMPBSA -O -i mmgbsa.in -cs md.tpr -ci index.ndx -cg 1 13 -ct md_noPBC.xtc -cp topol.top -nogui > mmgbsa.log 2>&1 &
```

| Option | Meaning |
|---|---|
| `-O` | overwrite results from an earlier run |
| `-i mmgbsa.in` | the input file from A4 |
| `-cs md.tpr` | the complex **structure** (the run input from Module 02) |
| `-ci index.ndx` | the complex **index** file |
| `-cg 1 13` | the receptor group and the ligand group (from A3) |
| `-ct md_noPBC.xtc` | the complex **trajectory** |
| `-cp topol.top` | the complex **topology**, which includes the JZ4 parameters |
| `-nogui` | do not try to open the results viewer (there is no screen in a codespace) |

This takes about 3 to 5 minutes. Follow its progress:

```bash
tail -f mmgbsa.log
```

Press **Ctrl + C** to stop watching (the calculation continues). It is
finished when `jobs` shows nothing and the log ends with timing
information, including `Total time taken`.

You may see this warning:

```
[WARNING] ... invalid DIHEDRAL_PERIODICITY = 0 found in Complex topology
```

It is harmless: gmx_MMPBSA corrects a formatting difference between
GROMACS and Amber topologies automatically.

gmx_MMPBSA creates many files. The ones starting with `_GMXMMPBSA_` are
temporary working files. The results are in:

| File | Contains |
|---|---|
| `FINAL_RESULTS_MMPBSA.dat` | the binding free energy and its components |
| `FINAL_DECOMP_MMPBSA.dat` | the per-residue contributions |
| `.csv` versions of both | the same numbers, per frame, for Excel or Google Sheets |

## A6. Read the binding free energy

```bash
sed -n '/Delta (Complex - Receptor - Ligand)/,/ΔTOTAL/p' FINAL_RESULTS_MMPBSA.dat
```

This prints the **Delta** section: the difference complex − receptor −
ligand, which is the binding energy. The columns that matter:

| Column | Meaning |
|---|---|
| `Average` | the average over the analysed frames |
| `SD` | the standard deviation: how much the value varies between frames |
| `SEM` | the standard error of the mean: how uncertain the average is |

Now plot it:

```bash
python plot_mmgbsa.py
```

Click `mmgbsa_terms.png`.

**What to expect for JZ4** (check this against your own results):

- **ΔVDWAALS** is the largest favourable term, because JZ4 sits in a
  hydrophobic pocket.
- **ΔEEL** is favourable but smaller: JZ4 has only one polar group, its OH.
- **ΔEGB** is **positive**: removing water from around the OH group
  costs energy.
- **ΔTOTAL** is negative: JZ4 binds.

This matches what you saw in Module 03: mostly hydrophobic contacts, and
the Lennard-Jones interaction energy larger than the Coulomb one.

## A7. Which residues contribute most?

`plot_mmgbsa.py` also made `mmgbsa_residues.png`. Click it.

Each bar is one residue's contribution to ΔG<sub>bind</sub>.
**Negative bars help binding**; positive bars oppose it. The terminal
lists the same residues with their values, and `mmgbsa_residues.csv`
contains all of them.

**Compare with Module 03.** Open your Module 03 contacts chart
(`../03-analysis/contacts.png`). Residues that are both **frequently in
contact** and **strongly favourable** in MM/GBSA are the **hot spots** of
the binding site: the residues a medicinal chemist would design new
ligands around, and the residues whose mutation would most weaken
binding.

> **👁 View it: the hot spots on the structure**
> Download `complex_start.pdb` from your Module 02 folder
> (`my-work/02-complex`) and open it in **PyMOL**. Replace the residue
> numbers below with your top 5 residues from the terminal output:
> ```
> hide everything
> show cartoon
> color grey80
> show sticks, resn JZ4
> select hotspots, resi 84+87+99+111+118
> show sticks, hotspots and not name N+C+O
> color orange, hotspots and elem C
> zoom resn JZ4, 6
> ```
> Where are the hot-spot residues relative to the ligand's ring and its
> OH group?

## A8. How reliable is the number?

Look at the `SEM` of ΔTOTAL. It tells you how much the average could
change if you analysed more frames from the same process. But it
**underestimates** the real uncertainty, because 21 frames from 100 ps
are not independent samples of all the ways the complex can move.

In a research study, the result is made more reliable by:

1. **Longer simulations** (tens to hundreds of ns), analysing frames
   spread over the equilibrated part.
2. **Replicas:** three or more independent simulations started with
   different random velocities (`gen_seed` in `nvt.mdp`). Report the
   average and spread across replicas.
3. **Checking the settings:** results depend on the GB model (`igb`),
   the radii and the dielectric constants. Use the same settings for
   every ligand you compare.

---

# Part B: Running research-scale simulations

How fast was your Module 02 production run? Check:

```bash
grep "Performance:" /workspaces/md-course/my-work/02-complex/md.log
```

On a 2-core codespace, the lysozyme-JZ4 system runs at roughly 10 ns
per day. A typical study needs 100 ns or more per replica, three replicas
per system, often for many ligands. That needs faster computers.

## B1. GPUs

A modern graphics card (GPU) runs GROMACS **10 to 50 times faster** than
a few CPU cores. On a GPU machine, the same `mdrun` command takes extra
options that move the work to the GPU:

```
gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu
```

| Option | Moves this work to the GPU |
|---|---|
| `-nb gpu` | short-range nonbonded forces (the most expensive part) |
| `-pme gpu` | long-range electrostatics (PME) |
| `-bonded gpu` | bonds, angles and dihedrals |
| `-update gpu` | updating positions and velocities each step |

Your `.mdp` files do not need to change. The GROMACS installed in your
codespace is a CPU-only version; GPU machines need a GPU-enabled build,
which HPC clusters usually provide.

## B2. HPC clusters and SLURM

A **high-performance computing (HPC) cluster** is a group of powerful
computers shared by many users. You do not run programs on it directly:
you write a **job script** describing what to run and what resources
you need, and submit it to a **scheduler**, which runs it when the
resources are free. The most common scheduler is **SLURM**.

Look at the example job script:

```bash
cat /workspaces/md-course/modules/04-going-further/files/slurm_md.sh
```

It has three parts:

1. **`#SBATCH` lines:** the resources requested: number of CPU cores,
   GPUs, memory and maximum run time. They look like comments, but
   SLURM reads them.
2. **`module load`:** loads the cluster's installed GROMACS.
3. **`gmx mdrun`:** the same command you have used all course, with GPU
   options and two important additions:
   - `-maxh` stops the run cleanly, with a checkpoint, just before the
     time limit;
   - `-cpi md.cpt` continues from that checkpoint, so submitting the same
     script again extends the run. This is how simulations of hundreds of
     ns run on clusters with 24- or 48-hour time limits.

The commands you will use on a cluster:

| Command | What it does |
|---|---|
| `sbatch slurm_md.sh` | submit the job; prints a job number |
| `squeue -u $USER` | list your jobs and whether they are waiting or running |
| `scancel 123456` | cancel job number 123456 |
| `sinfo` | list the partitions (queues) and their state |
| `module avail gromacs` | list the GROMACS versions installed |

Everything you prepared in this course (`.tpr` files, topologies,
index files) works unchanged on a cluster. A common workflow is to
**prepare and equilibrate on a small machine, then copy `md.tpr` to the
cluster for the long production run**, and copy the trajectory back for
analysis.

## B3. Where to find computing time

| Option | Notes |
|---|---|
| **Your university's HPC centre** | Usually free for students and staff; ask your supervisor or IT department how to get an account |
| **National research computing centres** | Many countries run shared HPC facilities for researchers; access is often through a short application |
| **Collaborators** | A research group with cluster access may run simulations as part of a collaboration |
| **Free notebook services with GPUs** (such as Google Colab or Kaggle) | Limited session lengths and weekly quotas, but usable for short production runs, using `-maxh` and checkpoints to continue across sessions |

Whatever you use, keep the same settings you used in this course:
`.mdp` files, force fields and the analysis workflow do not depend on
the computer.

---

# Part C: Starting your own project

Use this checklist when you set up a protein-ligand simulation for your
own research.

### 1. The protein structure

- Choose a high-resolution structure (ideally below 2.5 Å), preferably
  with your ligand or a similar one bound.
- Check for **missing residues** (`grep "REMARK 465"`) and **missing
  atoms** (`grep "REMARK 470"`). Gaps in the protein must be modelled
  (for example with PDBFixer, which is installed in your codespace,
  MODELLER or SWISS-MODEL) or carefully justified.
- Decide which **waters, ions and cofactors** to keep. Some are
  essential for binding.
- Check **protonation states** of histidine, aspartate, glutamate and
  lysine, especially near the binding site. Tools such as PROPKA or the
  H++ server predict them.

### 2. The ligand

- Check the ligand's **formula** after adding hydrogens (Module 02,
  Part 3).
- Decide its **protonation state and net charge at pH 7.4**. Amines,
  carboxylic acids and many heterocycles are charged. Give the correct
  charge to ACPYPE with `-n`.
- If the ligand comes from **docking** rather than a crystal structure,
  the starting pose is a prediction: watch the ligand RMSD closely, and
  consider simulating several docking poses.

### 3. The simulation

- Use a protein force field compatible with your ligand parameters
  (AMBER + GAFF2, as in this course, or CHARMM + CGenFF).
- Always minimize, then equilibrate with NVT and NPT, and **check** the
  temperature, pressure and density plots.
- Run **production long enough to answer your question** (typically
  100 ns or more), with **at least three replicas**.

### 4. The analysis

- Decide the equilibrated part from the RMSD first (Module 03, Part 3).
- Report averages **with their uncertainty** across replicas.
- Combine analyses: RMSD + contacts + hydrogen bonds + MM/GBSA together
  tell a more convincing story than any one of them.

### 5. Reporting

In the Methods section of a report or paper, include: the software and
versions, the force fields and water model, how the ligand was
parameterized (GAFF2, AM1-BCC charges, ACPYPE), the box type and
padding, the ion concentration, the equilibration protocol (length,
restraints, thermostat, barostat), the production length and number of
replicas, and the MM/GBSA settings. Cite the software you used. The
main tools in this course were GROMACS, AmberTools, ACPYPE, Open Babel,
MDAnalysis and gmx_MMPBSA.

---

## Common errors (Part A)

| Error message contains | Usual cause and fix |
|---|---|
| `The maximum number of iterations to select interaction residues was reached` | No residue is within 4 Å of the ligand: it has probably left the pocket. Check its RMSD and contacts (Module 03). To run anyway, change `print_res` to `"all"`. |
| `group ... not found` or wrong atoms selected | The group numbers after `-cg` are wrong. Check them with the `awk` command in A3. |
| `ERROR` mentioning the topology or an unknown atom type | `JZ4_GMX.itp` or `posre_JZ4.itp` is missing from this folder. Copy them (A2). |
| Energies of thousands of kcal/mol | The trajectory was not made whole. Use `md_noPBC.xtc`, not `md.xtc`. |
| `FINAL_RESULTS_MMPBSA.dat not found` (from `plot_mmgbsa.py`) | gmx_MMPBSA has not finished or stopped with an error. Read the end of `mmgbsa.log`. |

---

## Exercises

1. What is the MM/GBSA ΔTOTAL for JZ4, with its SD and SEM?
2. Which energy term favours binding most? Which term opposes it, and
   why?
3. List the five residues contributing most to binding. How many of them
   were also among your top five contacts in Module 03?
4. Change `interval = 5` to `interval = 1` in `mmgbsa.in` and run
   gmx_MMPBSA again. Does ΔTOTAL change? Does the SEM change? What does
   that tell you about using consecutive frames?
5. Change `igb = 5` to `igb = 2` and run again. How much does ΔTOTAL
   change? Why must you use the same settings when comparing ligands?
6. Looking at `slurm_md.sh`: what would you change to request 48 hours
   and 16 CPU cores? Why is `-maxh` set slightly lower than the time
   limit?
7. **Project plan:** choose a protein-ligand complex related to your own
   research from the PDB. Using the Part C checklist, write half a page
   describing how you would prepare and simulate it: which structure,
   what needs fixing, the ligand's charge, the simulation length and
   replicas, and which analyses would answer your question.

---

## Files you should have at the end

```
mmgbsa.in  mmgbsa.log
FINAL_RESULTS_MMPBSA.dat  FINAL_RESULTS_MMPBSA.csv
FINAL_DECOMP_MMPBSA.dat   FINAL_DECOMP_MMPBSA.csv
mmgbsa_terms.png  mmgbsa_residues.png  mmgbsa_residues.csv
```

Congratulations on completing the course. Remember to stop your
codespace, and download any results you want to keep before you delete
it.
