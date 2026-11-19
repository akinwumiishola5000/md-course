# Module 05: Simulating Your Own System

In Modules 01 to 04 you ran every step by hand, so that you understand
what each one does. In your research you will repeat the same steps for
many proteins and ligands, and typing about 50 commands each time is
slow and easy to get wrong.

This module gives you a **template**: one settings file and one script
that run the whole workflow of the course for **any protein, or any
protein-ligand complex**, from a PDB file to plots and a binding free
energy. You will:

1. test the template on the systems you already know (3HTB and 1AKI),
   and check that it reproduces your earlier results,
2. learn how to read its checks and fix the problems real structures have,
3. run it on **your own** protein or complex,
4. learn what to change to run a **longer, research-quality** simulation.

**Time needed:** about 2 hours, plus computing time
**Before you start:** complete Modules 01 to 04. The template does
nothing you have not already done by hand; the comments in it tell you
which module each step comes from.

---

## What is in the template

```bash
ls /workspaces/md-course/modules/05-your-own-system/template
```

| File | What it is | Do you edit it? |
|---|---|---|
| `config.sh` | All the settings: your PDB file, your ligand, simulation lengths, temperature... Every setting is explained inside the file. | **Yes**: this is the only file you edit |
| `run_md.sh` | The script that runs the workflow, in stages. It is heavily commented, so you can read exactly what it does. | No |

The script runs in **stages**, which match the modules:

| Stage | What it does | Same as |
|---|---|---|
| `check` | Inspects your structure and ligand, reports problems. **Changes nothing.** | Module 01 Part 1, Module 02 Parts 1–3 |
| `prepare` | Cleans the structure, builds the protein and ligand topologies, joins them, adds a box, water and ions | Module 01 Parts 1–5, Module 02 Parts 2–10 |
| `minimize` | Energy minimization, with a check that it converged | Module 01 Part 6 |
| `equilibrate` | NVT then NPT, with checks of temperature and density | Module 01 Parts 7–8 |
| `production` | Production MD. Prints how long it will take before starting | Module 01 Part 9 |
| `analyse` | RMSD, RMSF, radius of gyration; for a ligand also ligand RMSD, distance and contacts | Module 03 |
| `mmgbsa` | MM/GBSA binding free energy and per-residue contributions (ligand only) | Module 04 |
| `all` | All of the above, in order | |

It creates three folders:

| Folder | Contains |
|---|---|
| `run/` | every GROMACS file: structures, topology, `.mdp` files, trajectories. The same files you made by hand in the earlier modules. |
| `results/` | plots (`.png`), tables (`.csv`), structures to view, and `summary.txt` |
| `logs/` | the complete output of every command. When something fails, the script tells you which log to read. |

### What the script does for you that you did by hand

- It **writes the `.mdp` files** from your settings, with the same
  parameters as Modules 01–02, so lengths and temperature always match.
- It **edits the topology** to include the ligand, for one or many
  protein chains.
- It finds the **index group numbers** itself (for the thermostat and for
  MM/GBSA), so you never type the wrong number.
- It keeps **all protein chains and the ligand together** when it removes
  the periodic boundary effects (`gmx trjconv -pbc cluster`). In Modules
  02–03 you used `-pbc mol`, which keeps each molecule whole and works well
  for a buried ligand like JZ4. With a ligand in a shallow surface pocket,
  or a protein with several chains, `-pbc mol` can leave one molecule in a
  neighbouring periodic image, and the analysis would then see it far away.
- It **checks the results** of each stage and stops with an explanation
  in plain words if something is wrong.
- It can **continue** after an interruption: run the same command again
  and finished stages are skipped, and production continues from its
  checkpoint.

---

## Part 1: Set up a project folder

Each system gets its own folder with its own copy of the template.
Start with the system you know best, the T4 lysozyme-JZ4 complex:

```bash
cd /workspaces/md-course
```

```bash
git pull
```

```bash
mkdir -p my-work/05-3HTB
```

```bash
cd my-work/05-3HTB
```

```bash
cp /workspaces/md-course/modules/05-your-own-system/template/* .
```

```bash
ls
```

You should see `config.sh` and `run_md.sh`.

The settings in `config.sh` are already filled in for 3HTB with JZ4,
using the same short lengths as Module 02 (50 ps NVT, 50 ps NPT, 100 ps
production). You do not need a PDB file in the folder: when
`INPUT_PDB="3HTB.pdb"` is not found, the script downloads 3HTB from the
PDB automatically.

---

## Part 2: Read the settings file

```bash
less config.sh
```

(Use the arrow keys to scroll and **q** to quit.) The file has seven
sections. You will change these settings most often:

| Setting | What it controls | Example values |
|---|---|---|
| `INPUT_PDB` | your structure file, or a PDB ID to download | `"3HTB.pdb"`, `"my_protein.pdb"` |
| `LIGAND_RESNAME` | the ligand's 3-letter residue name; **empty for protein only** | `"JZ4"`, `"LIG"`, `""` |
| `LIGAND_FILE` | a separate ligand file, such as a docking pose; empty if the ligand is in `INPUT_PDB` | `""`, `"pose.sdf"` |
| `LIGAND_CHARGE` | the ligand's net charge; `auto` works it out at `LIGAND_PH` | `"auto"`, `"0"`, `"-1"` |
| `MD_PS` | production length in ps | `"100"`, `"100000"` (100 ns) |
| `SAVE_EVERY_PS` | how often a frame is saved | `"1"`, `"50"` |
| `TEMPERATURE` | temperature in K | `"300"`, `"310"` |

The rules for editing the file:

- Write `NAME="value"`, with **no spaces around `=`**.
- Keep the **quotes**.
- Lines starting with `#` are comments and are ignored.

---

## Part 3: Check the structure

Always run `check` first. It reads your structure and ligand, reports
what it finds, and **changes nothing** (except downloading the PDB file
if it is not in the folder yet):

```bash
bash run_md.sh check
```

Read the output from top to bottom. For 3HTB you should see:

- the **title** and **resolution** of the structure,
- the **protein chains** and how many residues each has,
- the **other molecules** (HETATM residues), with a reminder that all of
  them except your ligand will be removed,
- whether there are **gaps** in the chains,
- the **ligand**: its heavy atoms, its formula after adding hydrogens,
  and its **net charge**.

For JZ4 the formula must be **C9H12O** and the charge **0**, as in
Module 02.

### What to do with each kind of warning

| `check` reports... | What it means | What to do |
|---|---|---|
| `Alternate locations found` | some atoms have two possible positions (A and B) | Nothing: position A is kept automatically |
| `The file has N models` | an NMR structure, which contains many models | Nothing: the first model is used. Consider whether a crystal structure exists |
| A modified amino acid in the HETATM list (`MSE`, `SEP`, `TPO`, `PTR`...) | non-standard residues the force field does not know | Keep `FIX_PROTEIN="yes"`: PDBFixer converts them to the standard amino acid |
| `Gaps in the protein chain` | residues are missing in the middle of a chain | `prepare` stops until you choose what to do: set `MISSING_RESIDUES`, or use a complete structure. See "Gaps in the chain" below |
| `cis peptide bond(s) not present in the input structure` (after `minimize` or `equilibrate`) | a rare and usually wrong geometry, typically in rebuilt residues | See "Gaps in the chain" below |
| `N copies of ... found` | the ligand appears more than once (often one per protein chain) | Only the first copy is simulated. Consider keeping only the chain you need (see "Keeping only one chain") |
| `No HETATM residue called ...` | the ligand name is wrong | Look at the HETATM list and correct `LIGAND_RESNAME` |
| A charge you did not expect | the ligand's protonation at pH 7.4 | Check the chemistry (below), and set `LIGAND_CHARGE` if needed |
| `LIGAND_CHARGE=... but the protonated structure has charge ...` | your setting disagrees with the structure | Decide which is right; a wrong charge gives a wrong simulation |

### Is the ligand charge right?

The net charge is the most common source of errors with real ligands.
At pH 7.4:

| Group in the ligand | Usually |
|---|---|
| carboxylic acid (–COOH) | deprotonated, **−1** |
| aliphatic amine (–NH₂, –NHR, –NR₂) | protonated, **+1** |
| amidine, guanidine | protonated, **+1** |
| phosphate, sulfonate | negative (**−1** or **−2**) |
| aniline, amide, pyridine, phenol, alcohol | **neutral** |

Add up the groups in your ligand and compare with the charge `check`
reports. If they disagree, set `LIGAND_CHARGE` yourself, or provide the
ligand in the correct protonation state as an SDF or MOL2 file
(`LIGAND_FILE`). Tools such as the free web server of Chemicalize,
MarvinSketch or the `dimorphite_dl` Python package predict protonation
states.

---

## Part 4: Run the simulation

The first time, run the stages **one at a time**, and look at the plots
after each one. Then you know exactly where a problem starts.

### Stage 1: prepare

```bash
bash run_md.sh prepare
```

This takes 1 to 5 minutes, mostly for the ligand charges (ACPYPE).
It reports the number of atoms, waters and ions at the end.

### Stage 2: minimize

```bash
bash run_md.sh minimize
```

Click `results/1_potential.png`. The energy should fall and level off,
as in Module 01.

### Stage 3: equilibrate

The equilibration runs take several minutes, so run them in the
background, as in Module 01:

```bash
nohup bash run_md.sh equilibrate > equilibrate.log 2>&1 &
```

Watch the progress (press **Ctrl + C** to stop watching; the run
continues):

```bash
tail -f equilibrate.log
```

When it has finished, the log shows the average temperature, pressure
and density, and the plots are `results/2_temperature.png`,
`3_pressure.png` and `4_density.png`. Compare them with Module 02.

### Stage 4: production

```bash
nohup bash run_md.sh production > production.log 2>&1 &
```

```bash
tail -f production.log
```

Before starting, the script prints the **speed** of your system on your
computer and the **estimated time** for the production run.

### Stages 5 and 6: analysis and MM/GBSA

```bash
bash run_md.sh analyse
```

```bash
bash run_md.sh mmgbsa
```

### Or everything at once

Once you trust your settings, one command runs every stage:

```bash
nohup bash run_md.sh all > run.log 2>&1 &
```

```bash
tail -f run.log
```

If anything stops it (an error you then fix, or your codespace going to
sleep), run the same command again: finished stages are skipped, and
production continues from its last checkpoint.

---

## Part 5: The results

```bash
ls results
```

| File | Contents | From |
|---|---|---|
| `summary.txt` | the main numbers from every stage, in one place | |
| `1_potential.png` | energy minimization | Module 01 |
| `2_temperature.png`, `3_pressure.png`, `4_density.png` | equilibration checks | Module 01 |
| `5_rmsd.png` | backbone RMSD (and ligand RMSD) | Module 03, Part 3 |
| `6_rmsf.png` | flexibility of each residue | Module 03, Part 4 |
| `7_gyrate.png` | radius of gyration | Module 03, Part 5 |
| `8_mindist.png` | shortest protein-ligand distance | Module 02 |
| `9_contacts.png`, `contacts.csv`, `hbonds.csv` | residues in contact with the ligand, hydrogen bonds | Module 03, Part 8 |
| `mmgbsa/mmgbsa_terms.png`, `mmgbsa/mmgbsa_residues.png` | binding free energy and hot-spot residues | Module 04 |
| `start.pdb`, `trajectory.xtc` | protein (and ligand) without water, fitted, for viewing | Module 03, Part 2 |
| `rmsf_bfactors.pdb` | structure coloured by flexibility | Module 03, Part 4 |
| `*.xvg` | the data behind every plot | |

Read the summary:

```bash
cat results/summary.txt
```

**Validation exercise:** compare the numbers for 3HTB with your own
results from Modules 02, 03 and 04. They will not be identical, because
every simulation starts from different random velocities, but they
should be similar. This is how you check that a new workflow is
correct before trusting it with a new system.

> **👁 View it: the trajectory and the flexibility**
> Download `results/start.pdb` and `results/trajectory.xtc`, and play
> them in VMD or PyMOL exactly as in Module 02, Part 15.
> Download `results/rmsf_bfactors.pdb` and colour it by flexibility in
> PyMOL as in Module 03, Part 4:
> ```
> spectrum b, blue_white_red
> cartoon putty
> ```

---

## Part 6: Your own system

### Step 1: a new folder and your structure

Make a folder for your system (change the name to suit it) and copy the
template:

```bash
mkdir -p /workspaces/md-course/my-work/05-mysystem
```

```bash
cd /workspaces/md-course/my-work/05-mysystem
```

```bash
cp /workspaces/md-course/modules/05-your-own-system/template/* .
```

Put your structure in this folder:

- **From the PDB:** you do not need to download it. Set
  `INPUT_PDB="XXXX.pdb"` with its 4-character ID, and `check` downloads it.
- **From your computer** (a model, a prepared structure or a docking
  result): in the file panel, open `my-work`, then your folder, and
  **drag the file from your computer onto the folder**. It uploads into
  your codespace.

### Step 2: edit the settings

```bash
nano config.sh
```

Choose the case that matches your system:

**A protein on its own:**
```
INPUT_PDB="my_protein.pdb"
LIGAND_RESNAME=""
```

**A protein-ligand complex from the PDB (ligand in the same file):**
```
INPUT_PDB="XXXX.pdb"
LIGAND_RESNAME="ABC"     <- the ligand's name from the HETATM list
LIGAND_FILE=""
```

**A docking result (protein and ligand in separate files):**
```
INPUT_PDB="receptor.pdb"
LIGAND_RESNAME="LIG"
LIGAND_FILE="pose.sdf"   <- the docked pose: .sdf, .mol2 or .pdb
```
The ligand file must contain **only the pose you want** (or have it
first), already **in the binding site**: docking programs write poses in
the receptor's coordinates, so this is normally true. If your docking
program wrote `.pdbqt` files, convert the chosen pose first:

```bash
obabel docked.pdbqt -l 1 -O pose.sdf
```

(`-l 1` keeps only the first pose. Change the file names to yours.)

Save with **Ctrl + O**, **Enter**, and exit with **Ctrl + X**.

### Step 3: check, fix, run

```bash
bash run_md.sh check
```

Fix anything it reports (Part 3 and the sections below), run `check`
again until it looks right, then run the stages as in Part 4.

---

## Part 7: Running longer, research-quality simulations

The default lengths **teach the workflow but are too short for research**.
For a study you will publish, change these settings in `config.sh`:

| Setting | Course default | Research | Why |
|---|---|---|---|
| `NVT_PS` | 50 | 100–500 | larger systems need longer to reach the temperature |
| `NPT_PS` | 50 | 100–1000 | the density must be stable before production |
| `MD_PS` | 100 | **100000 (100 ns) or more** | slow motions and ligand stability need nanoseconds |
| `SAVE_EVERY_PS` | 1 | **10–100** | keeps the trajectory a manageable size: aim for 1000–5000 frames |
| `MMGBSA_FRAMES` | 20 | 100–500 | more frames give a more reliable average |

**Always change `SAVE_EVERY_PS` together with `MD_PS`.** 100 ns saved every
1 ps would be 100,000 frames: a trajectory of many gigabytes that is slow
to analyse. With `SAVE_EVERY_PS="50"`, 100 ns gives 2,000 frames.

### How long will it take?

The `production` stage measures your system's speed and prints an
estimate before starting. You can also calculate it yourself from any
earlier run of the same system:

```
hours = (MD_PS ÷ 1000) ÷ (ns/day) × 24
```

On a 2-core codespace, a system like 3HTB runs at roughly 10 ns/day, so
100 ns takes about 10 days. That is why research simulations run on GPUs
or HPC clusters (Module 04, Part B), where the same system runs at 100 to
500 ns/day.

### A codespace for preparation, a GPU for production

A practical approach is to do everything **except production** in your
codespace, and only production on a faster computer:

1. In your codespace, set the research lengths in `config.sh`, and run:
   ```bash
   bash run_md.sh prepare
   ```
   ```bash
   bash run_md.sh minimize
   ```
   ```bash
   nohup bash run_md.sh equilibrate > equilibrate.log 2>&1 &
   ```
2. Pack the whole project folder into one file (run this in the
   folder **above** it, `my-work`):
   ```bash
   cd /workspaces/md-course/my-work
   ```
   ```bash
   tar czf 05-mysystem.tar.gz 05-mysystem
   ```
   Download `05-mysystem.tar.gz` (right-click → **Download**) and copy it
   to the GPU computer or cluster. Unpack it there with
   `tar xzf 05-mysystem.tar.gz`.
3. On the GPU computer or cluster, the production stage only needs
   GROMACS and Python 3. In `config.sh`, set:
   ```
   THREADS="auto"
   MDRUN_OPTIONS="-nb gpu -pme gpu -bonded gpu -update gpu"
   ```
   If the cluster's GROMACS command is `gmx_mpi` rather than `gmx`
   (check with `module avail gromacs` and your cluster's documentation),
   also set `GMX_COMMAND="gmx_mpi"`. Then run `bash run_md.sh production`,
   for example inside the SLURM script from Module 04: replace its
   `gmx mdrun ...` line with `bash run_md.sh production`.
4. Pack the folder again on the cluster, upload the `.tar.gz` file into
   `my-work` in your codespace (drag it onto the folder in the file
   panel), unpack it with `tar xzf`, and run `bash run_md.sh analyse`
   and `bash run_md.sh mmgbsa` there.

If a cluster job reaches its time limit, submit it again: production
continues from the checkpoint.

### Replicas

Research results need **at least three independent simulations**
(replicas), so you can report an average and its spread. Each run starts
from different random velocities automatically, so a replica is simply
another copy of the folder:

```bash
cd /workspaces/md-course/my-work
```

```bash
mkdir 05-mysystem-rep2
```

```bash
cp 05-mysystem/config.sh 05-mysystem/run_md.sh 05-mysystem-rep2/
```

Copy your structure (and ligand file) into the new folder too, then run
`bash run_md.sh all` in it. Compare the `summary.txt` files of all
replicas.

---

## Special situations

### Gaps in the chain

Crystal structures often lack residues that were too flexible to be
seen, usually in loops. If `check` reports **gaps**, residues are missing
in the **middle** of a chain.

**Why gaps matter.** GROMACS does not know residues are missing. It
would join the two ends of a gap with one very long "peptide bond", and
energy minimization would then drag the ends together and distort the
protein around the gap, **with no error message**. That is why, by
default (`MISSING_RESIDUES="stop"`), `prepare` stops when it finds a gap,
and asks you to choose what to do.

Missing residues at the **start or end** of a chain are a different case:
they are usually disordered tails, the chain simply starts or ends a
little later, and nothing is added there.

**Your options, from best to simplest:**

| Option | How | When to use it |
|---|---|---|
| **1. A complete structure** | Use another PDB entry of the same protein without the gap, or a complete model: an AlphaFold model (alphafold.ebi.ac.uk), or the missing residues modelled with SWISS-MODEL or MODELLER. Use it as `INPUT_PDB`. | **Always the best choice**, and the only good one if the gap is **near your binding site** |
| **2. Rebuild short gaps:** `MISSING_RESIDUES="fill"` | PDBFixer rebuilds gaps of up to `MAX_FILL_LENGTH` residues (default 5), using the sequence in the file's `SEQRES` records | Short gaps (a few residues) **far from the binding site** |
| **3. Leave the gap open:** `MISSING_RESIDUES="split"` | Nothing is added. Each side of the gap becomes a free chain end, so no false bond is made | Gaps **far from the binding site**, when you have no complete structure; also for long gaps |

**About option 2.** PDBFixer gets the sequence right, but builds the loop
quickly and roughly. In testing on lysozyme, a rebuilt 3-residue loop
contained a **cis peptide bond**, a geometry that is very rare in real
proteins. The script therefore checks for new cis peptide bonds after
minimization and after equilibration, and warns you if it finds any:

```
[WARNING] cis peptide bond(s) not present in the input structure:
  SER72-ARG73
```

If you see this warning in a region that matters for your question, use
option 1 instead. If the region is far from the binding site, you can
usually continue, and mention it in your methods.

`fill` needs `SEQRES` records, which every file downloaded from the PDB
has. Models from other sources (docking receptors, homology models)
often do not have them: then `prepare` stops and tells you so, and you
can use option 1 or 3.

**About option 3.** The two new chain ends are charged (–NH₃⁺ and –COO⁻),
as real chain ends are. That is a small change to the chemistry near
the gap, which is why this option is for gaps far from the binding site.
GROMACS lists the pieces as separate molecules (for example
`Protein_chain_A` and `Protein_chain_A2`); the analysis keeps them
together.

### Keeping only one chain

Many structures contain several copies of the protein (chains A, B, C...)
because of how the crystal packs, not because the protein works that
way. Simulating all of them makes the system several times larger and
slower. To keep only chain A (and anything else labelled chain A, such
as its ligand):

```bash
awk '/^(ATOM|HETATM)/ && substr($0,22,1)=="A"' XXXX.pdb > XXXX_chainA.pdb
```

The chain letter is in column 22 of every `ATOM` and `HETATM` line;
`awk` keeps the lines where it is `A`. Change `XXXX` to your file name,
and `"A"` to another chain if you need it.

Then set `INPUT_PDB="XXXX_chainA.pdb"`. Check the biological assembly
on the PDB website first: if the protein really works as a dimer, keep
both chains.

### Ions, metals and cofactors

`prepare` **removes every HETATM molecule except your ligand**. That is
right for crystallization additives (sulfate, glycerol, PEG...), but
not always for:

- **Structural metal ions** (zinc in zinc fingers, calcium in
  calcium-binding sites): without them, the site can fall apart. They
  need special parameters, beyond this template.
- **Cofactors** (heme, NAD, FAD, ATP): they are part of the active site.
  A cofactor can be simulated as the ligand instead (set
  `LIGAND_RESNAME` to its name), but only one HETATM molecule can be
  kept by this template.

If your protein depends on a metal or cofactor, discuss it with your
instructor before trusting the results.

### What this template cannot do

| System | Why | Where to look |
|---|---|---|
| Membrane proteins | need a lipid bilayer | CHARMM-GUI Membrane Builder |
| More than one ligand, or ligand + cofactor | the template keeps one HETATM molecule | build the topology by hand, as in Module 02, once per molecule |
| Covalently bound ligands | the ligand is part of the protein | special parameterization |
| Metal-containing ligands, or elements GAFF2 does not cover (B, Si, metals) | no GAFF2 parameters | MCPB.py (AmberTools), specialist force fields |
| DNA, RNA | not tested with this template | AMBER nucleic acid force fields |
| Very large ligands (more than ~100 atoms) | AM1-BCC charges become very slow | patience, or other charge methods |

---

## Troubleshooting

When a stage fails, the script stops and shows:

1. **`[FAILED]`** and the name of the step,
2. the **last lines of its output**,
3. **"What this usually means"**: an explanation and what to do,
4. the path of the **full log** in `logs/`.

To read a full log:

```bash
less logs/pdb2gmx.log
```

Common problems:

| Problem | Fix |
|---|---|
| `residue ... not found in residue topology database` | A non-protein molecule is still in the protein. Check the HETATM list (`check`), set `LIGAND_RESNAME`, keep `FIX_PROTEIN="yes"`. |
| Atoms missing in a residue | Set `FIX_PROTEIN="yes"`. |
| ACPYPE fails | Usually a wrong ligand charge, or a ligand file without proper bonds. Check the formula and charge in `check`. |
| `Maximum force ... is far too high` after minimize | Atoms overlap: a gap in the chain, a ligand pose that clashes with the protein, or a wrong ligand structure. |
| Temperature or density warnings | Equilibration was too short for this system. Increase `NVT_PS` or `NPT_PS` and run `equilibrate` again. |
| `Run the 'X' stage first` | Stages must run in order. Run the stage named. |
| MM/GBSA: `no protein residue within 4 A` | The ligand left the binding site during the simulation. Look at `5_rmsd.png` and `8_mindist.png`. |
| You changed `config.sh` and want to start again | Run `bash run_md.sh prepare`. An earlier finished run is kept in a `run_backup_...` folder. |
| You changed only `MD_PS` or `SAVE_EVERY_PS` | Run `bash run_md.sh prepare` again too: the `.mdp` files are written during `prepare`. |
| You want to make a **finished** production run longer, without starting again | Extend it as in Module 01 ("Extending a finished run"), inside the `run/` folder: `gmx convert-tpr -s md.tpr -extend 400 -o md.tpr` (400 ps more), then `bash run_md.sh production`. Then run `analyse` again. |

If you are stuck, post in the help channel: the command you ran, the
`[FAILED]` message, and your `config.sh`.

---

## Exercises

1. Run the template on 3HTB with the default settings. Compare the
   ligand RMSD, the top five contacts and the MM/GBSA ΔG with your results
   from Modules 03 and 04. Which numbers are close, and which differ
   more? Why can they differ at all?
2. Make a second folder for **1AKI, protein only** (`INPUT_PDB="1AKI.pdb"`,
   `LIGAND_RESNAME=""`). Run it and compare the backbone RMSD and radius of
   gyration with Module 01.
3. Look inside `run/md.mdp` for the 3HTB run. Which settings came from
   `config.sh`? Which are the same as in Module 02?
4. Change `TEMPERATURE` to `310` in a copy of the 3HTB folder. Which `.mdp`
   files and which lines change? (Run `prepare`, then compare
   `run/nvt.mdp` in the two folders with `diff`.)
5. You want to simulate a complex for 200 ns with about 4,000 frames. What
   values do you set for `MD_PS` and `SAVE_EVERY_PS`? Using the speed your
   3HTB run reported, how long would it take in your codespace?
6. **Your project:** run the template on a protein or protein-ligand
   complex from your own research, with the default (short) lengths. Write
   half a page: what `check` reported and how you fixed it, the main
   numbers from `summary.txt`, and which settings you would change for a
   research-quality study (Part 7).

---

Remember to stop your codespace when you finish.
