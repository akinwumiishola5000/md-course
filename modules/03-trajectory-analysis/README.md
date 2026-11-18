# Module 03: Analysing MD Trajectories

A trajectory is a movie of thousands of atoms. On its own it answers
nothing: **analysis** turns it into answers. In this module you analyse
the protein-ligand simulation from Module 02 and ask the questions
drug-discovery scientists ask of every complex:

1. Has the simulation **equilibrated**, so that the results mean something?
2. Is the **protein stable**, and which parts of it are **flexible**?
3. Does the **ligand stay bound**, and in what **pose**?
4. **Which residues** hold the ligand in place, through which **interactions**?
5. How strong are those **interactions**?
6. What are the **main motions** of the protein?

**Time needed:** about 3 hours
**Before you start:** complete Module 02. You also need PyMOL or VMD
on your computer (see the start of Module 02).

---

## Your trajectory

You analyse the 100 ps production run you made in Module 02: 100 ps of
the T4 lysozyme-JZ4 complex, saved as one frame every 1 ps.

> **An honest warning about 100 ps.** Proteins move on time scales from
> femtoseconds to seconds. 100 ps captures fast motions only: side-chain
> wiggles, loop breathing, the first adjustments of the ligand. Every
> analysis in this module works the same way on 100 ns as on 100 ps,
> which is why we learn it here, but **do not draw scientific
> conclusions from a 100 ps simulation**. Published studies typically
> use 100 ns or more, often repeated three times.

---

## The analyses in this module

| Part | Analysis | GROMACS tool | Question it answers |
|---|---|---|---|
| 2 | Prepare the trajectory | `trjconv`, `convert-tpr` | (makes everything else faster and correct) |
| 3 | RMSD | `gmx rms` | Has the system equilibrated? Is the ligand pose stable? |
| 4 | RMSF | `gmx rmsf` | Which residues are flexible? |
| 5 | Radius of gyration | `gmx gyrate` | Does the protein stay compact? |
| 6 | Solvent-accessible surface area | `gmx sasa` | How buried is the ligand? |
| 7 | Ligand-pocket distance | `gmx distance` | Does the ligand drift away from the pocket? |
| 8 | Contacts and hydrogen bonds | `ligand_contacts.py` | Which residues hold the ligand, and how? |
| 9 | Interaction energy | `gmx mdrun -rerun` | How strongly does the ligand interact with the protein? |
| 10 | Clustering | `gmx cluster` | What are the main ligand poses? |
| 11 | Principal component analysis | `gmx covar`, `gmx anaeig` | What are the main motions of the protein? |

---

## Part 0: Get ready

```bash
cd /workspaces/md-course
```

```bash
git pull
```

```bash
mkdir -p my-work/03-analysis
```

```bash
cd my-work/03-analysis
```

Copy the analysis tools:

```bash
cp /workspaces/md-course/modules/tools/plot_xvg.py .
```

```bash
cp /workspaces/md-course/modules/tools/ligand_contacts.py .
```

`plot_xvg.py` has new options in this module, in addition to those
from Module 02:

| Option | Effect |
|---|---|
| `--bar` | bar chart instead of a line |
| `--scatter` | one point per frame, coloured from start to end of the simulation |
| `--cols 2` | plot only the chosen data column(s) |
| `--xmax 20` | plot only points up to x = 20 |

Copy your results from Module 02. The `..` means "the folder above", so
`../02-complex` is your Module 02 folder:

```bash
cp ../02-complex/md.tpr .
```

```bash
cp ../02-complex/md.xtc .
```

```bash
cp ../02-complex/md_noPBC.xtc .
```

```bash
cp ../02-complex/index.ndx .
```

For the interaction energy in Part 9 you also need the topology files:

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
cp ../02-complex/md.mdp .
```

```bash
cp ../02-complex/md.gro .
```

```bash
ls
```

> **Did not finish Module 02?** Use the reference files instead. Replace
> `../02-complex/` in each command above with
> `/workspaces/md-course/modules/02-protein-ligand/reference/`, and
> create `md_noPBC.xtc` with the first command in Part 2.

---

## Part 1: What is in a trajectory?

Before analysing anything, check what the trajectory contains:

```bash
gmx check -f md.xtc
```

Read the summary at the end:

| Line | Meaning |
|---|---|
| `Coords` | number of frames that contain coordinates |
| `Last frame ... time` | the time of the last frame, in ps |
| the number of atoms | how many atoms each frame has |

Your trajectory should have about 101 frames (time 0 to 100 ps), one
every 1 ps. That spacing was set by `nstxout-compressed = 500` in
`md.mdp`: 500 steps × 0.002 ps = 1 ps.

**Why it matters:** every analysis result is an average over frames.
More frames, and especially more **independent** frames spread over a
longer time, give more reliable averages.

---

## Part 2: Prepare the trajectory for analysis

Raw trajectories have three problems for analysis:

1. **Broken molecules.** Because of periodic boundaries, part of the
   protein can appear on the other side of the box.
2. **Tumbling.** The whole protein rotates and drifts in the box, which
   hides its internal motions.
3. **Size.** 95% of the atoms are water, which you do not need for most
   analyses and which makes everything slow.

### Step 1: make molecules whole and centre the protein

You already did this in Module 02, producing `md_noPBC.xtc`. If you do
not have that file, run:

```bash
echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_noPBC.xtc -pbc mol -center
```

### Step 2: remove tumbling and remove water

```bash
echo "Backbone Protein_JZ4" | gmx trjconv -s md.tpr -f md_noPBC.xtc -n index.ndx -o md_fit.xtc -fit rot+trans
```

The two answers:

- **`Backbone`**: superimpose every frame on the first one, using the
  backbone atoms. This removes overall rotation and translation (`-fit
  rot+trans`), leaving only internal motion.
- **`Protein_JZ4`**: write out only the protein and the ligand.

### Step 3: make a matching run input file

GROMACS tools need a `.tpr` file that matches the atoms in the
trajectory. `md.tpr` contains the water too, so it no longer matches
`md_fit.xtc`. Make a smaller one:

```bash
echo "Protein_JZ4" | gmx convert-tpr -s md.tpr -n index.ndx -o complex.tpr
```

### Step 4: make a matching index file

```bash
echo "q" | gmx make_ndx -f complex.tpr -o complex.ndx
```

`q` quits immediately, saving the default groups for the smaller
system: `Protein`, `Backbone`, `C-alpha`, `JZ4`, and more.

Compare the sizes:

```bash
ls -lh md_noPBC.xtc md_fit.xtc
```

The new trajectory is about 10 times smaller. **From now on, every
analysis uses `complex.tpr`, `md_fit.xtc` and `complex.ndx`**, except
the interaction energy in Part 9.

### Step 5: a structure file for viewers

Molecular viewers cannot read `.tpr` files. To watch `md_fit.xtc` in
PyMOL or VMD, they need a structure file with the same atoms. Save the
first frame as a PDB file:

```bash
echo "System" | gmx trjconv -s complex.tpr -f md_fit.xtc -o start.pdb -dump 0
```

You will use `start.pdb` several times in this module.

> **👁 View it: the prepared trajectory**
> Download `start.pdb` and `md_fit.xtc`, load them in VMD or PyMOL as in
> Module 02, Part 15, and play the movie. Compared with Module 02, the
> protein no longer drifts or rotates, and there is no water: only the
> internal motions of the protein and the ligand remain.

---

## Part 3: RMSD — has the system equilibrated?

**RMSD** (root mean square deviation) is the average distance between
the atoms of each frame and the same atoms in a reference structure,
after superimposing them:

```
RMSD = √( Σ (position now − reference position)² / number of atoms )
```

You met RMSD in Module 02. Here you use it for its most important job:
deciding **which part of the trajectory to analyse**.

```bash
echo "Backbone Backbone" | gmx rms -s complex.tpr -f md_fit.xtc -n complex.ndx -o rmsd_backbone.xvg -tu ps
```

```bash
echo "Backbone JZ4" | gmx rms -s complex.tpr -f md_fit.xtc -n complex.ndx -o rmsd_ligand.xvg -tu ps
```

```bash
python plot_xvg.py rmsd_backbone.xvg rmsd_ligand.xvg --unit A --avg 10 -o rmsd.png --title "Backbone and ligand RMSD"
```

Click `rmsd.png`.

**What you should see:** the RMSD rises from zero at the start (the
structure relaxes away from the crystal structure), then levels off
and fluctuates around a stable value. The point where it levels off is
where the system has **equilibrated**.

**How to use it:** analyses of averages (RMSF, contacts, energies,
clusters) should use only the equilibrated part. Every GROMACS analysis
tool accepts `-b` (begin) to skip the start. For example, if the RMSD
levels off at 20 ps, add `-b 20` to the analysis commands. In a 100 ns
study, it is common to skip the first 10 to 20 ns.

> In this module we analyse the whole 100 ps, because there is so little
> of it. In your research, always decide the equilibrated part from the
> RMSD first.

---

## Part 4: RMSF — which residues are flexible?

**RMSF** (root mean square fluctuation) is calculated **per atom** or
**per residue**: how much each one fluctuates around its own **average**
position over the whole trajectory.

| | RMSD | RMSF |
|---|---|---|
| Calculated for | each **frame** | each **residue** |
| Plotted against | time | residue number |
| Tells you | how the whole structure changes over time | which parts move most |

```bash
echo "C-alpha" | gmx rmsf -s complex.tpr -f md_fit.xtc -n complex.ndx -o rmsf.xvg -res -oq bfactors.pdb
```

| Option | Meaning |
|---|---|
| `C-alpha` | use one atom per residue, the alpha carbon |
| `-res` | report one value per residue |
| `-oq bfactors.pdb` | also save a structure with each residue's RMSF stored in the B-factor column, for colouring in PyMOL |

```bash
python plot_xvg.py rmsf.xvg --unit A --xlabel "Residue number"
```

**What you should see:** low values (below about 1 Å) for most of the
protein, and peaks at a few places. Peaks are usually:

- the **N- and C-terminal** residues (free ends),
- **loops** on the protein surface,

while **α-helices and β-sheets**, especially buried ones, are rigid.

**Why it matters in drug discovery:** if pocket residues are very
flexible, the pocket may change shape, which affects how ligands bind.
Comparing the RMSF of the protein with and without a ligand shows
whether the ligand **stabilises** parts of the protein.

> **👁 View it: flexibility on the structure**
> Download `bfactors.pdb` and open it in **PyMOL**. Type:
> ```
> hide everything
> show cartoon
> spectrum b, blue_white_red
> cartoon putty
> ```
> Blue and thin = rigid; red and thick = flexible. Where are the flexible
> regions? Is the binding pocket rigid or flexible?
>
> In **ChimeraX**: open the file and type `color bfactor palette bluered`.

---

## Part 5: Radius of gyration — does the protein stay compact?

The **radius of gyration (Rg)** measures how spread out the protein's
atoms are around its centre of mass. A folded protein keeps an almost
constant Rg; unfolding makes it increase.

```bash
echo "Protein" | gmx gyrate -s complex.tpr -f md_fit.xtc -n complex.ndx -o gyrate.xvg
```

```bash
python plot_xvg.py gyrate.xvg --unit A
```

The plot has four lines: the total Rg and its components around the x,
y and z axes. Focus on the **total**.

**What you should see:** a flat line, changing by less than about 0.2 Å
over 100 ps. A steady increase would suggest the protein is unfolding,
usually a sign of a problem in the setup.

---

## Part 6: SASA — how buried is the ligand?

The **solvent-accessible surface area (SASA)** is the area of a
molecule's surface that a water molecule could touch. A ligand deep in
a pocket has a small SASA; the same ligand floating in water has a
large one.

Calculate the ligand's SASA **in the complex**, where the protein
shields part of it:

```bash
gmx sasa -s complex.tpr -f md_fit.xtc -n complex.ndx -surface 'group "System"' -output 'group "JZ4"' -o sasa_bound.xvg
```

Then the SASA of the ligand **on its own**, as if the protein were not there:

```bash
gmx sasa -s complex.tpr -f md_fit.xtc -n complex.ndx -surface 'group "JZ4"' -output 'group "JZ4"' -o sasa_free.xvg
```

| Option | Meaning |
|---|---|
| `-surface` | the atoms that make up the surface being measured |
| `-output` | the part of that surface to report |

Each file has two data columns: the **total** area of the `-surface`
atoms, and the area of the `-output` atoms (`JZ4`). You only want the
ligand, which is the **second** data column. `--cols 2` plots only that
column:

```bash
python plot_xvg.py sasa_bound.xvg sasa_free.xvg --cols 2 -o ligand_sasa.png --title "Ligand SASA: bound vs free"
```

Click `ligand_sasa.png`. The terminal prints the mean of each line, in
nm². Calculate the **fraction buried**:

```
fraction buried = 1 − (SASA bound ÷ SASA free)
```

**What it means:** above 0.8 (80%), the ligand is deeply buried, which
is typical of tight binders in enclosed pockets. Close to 0, the
ligand is barely touching the protein.

---

## Part 7: Does the ligand drift away from the pocket?

RMSD tells you the ligand moved, but not whether it moved **within**
the pocket or **out of** it. The distance between the ligand's centre
and the pocket's centre answers that.

### Step 1: define the pocket

The pocket is every protein residue with an atom within 0.5 nm (5 Å) of
the ligand, in the first frame:

```bash
gmx select -s complex.tpr -f md_fit.xtc -n complex.ndx -select 'same residue as (group "Protein" and within 0.5 of resname JZ4)' -on pocket.ndx -e 0
```

Reading the selection from the inside out:

| Part | Meaning |
|---|---|
| `resname JZ4` | the ligand atoms |
| `within 0.5 of ...` | atoms closer than 0.5 nm to any ligand atom |
| `group "Protein" and ...` | ...that belong to the protein |
| `same residue as (...)` | expand to whole residues, not just the touching atoms |
| `-e 0` | use only the first frame (time 0), so the pocket stays fixed |
| `-on pocket.ndx` | save the selected atoms as an index group |

```bash
head -n 3 pocket.ndx
```

The file has one group (group 0) containing the pocket's atom numbers.

### Step 2: measure the distance

```bash
gmx distance -s complex.tpr -f md_fit.xtc -n pocket.ndx -select 'com of resname JZ4 plus com of group 0' -oall pocket_distance.xvg -tu ps
```

- `com of resname JZ4`: the ligand's centre of mass
- `com of group 0`: the centre of mass of the pocket (group 0 in `pocket.ndx`)

```bash
python plot_xvg.py pocket_distance.xvg --unit A --avg 10
```

**What you should see:** a roughly constant distance. A ligand that is
leaving its pocket shows a distance that grows steadily; typically by
several Å before it loses contact completely.

---

## Part 8: Which residues hold the ligand?

This is where MD adds most to docking: instead of one snapshot, you see
how **often** each interaction exists.

`ligand_contacts.py` is a Python script using the **MDAnalysis**
library. For every frame, it finds:

- **contacts**: protein residues with any heavy atom within 4 Å of any
  ligand heavy atom, and
- **hydrogen bonds**: donor-acceptor distance ≤ 3.0 Å and
  donor-hydrogen-acceptor angle ≥ 150°.

It then reports each one as a **percentage of frames** (its
**occupancy**).

```bash
python ligand_contacts.py complex.tpr md_fit.xtc --ligand JZ4
```

It prints the results and saves three files:

| File | Contents |
|---|---|
| `contacts.png` | bar chart of the residues most often in contact with the ligand |
| `contacts.csv` | the same as a table (opens in Excel or Google Sheets) |
| `hbonds.csv` | each hydrogen bond, with its occupancy |

Click `contacts.png`.

**How to read it:**

| Occupancy | Meaning |
|---|---|
| Above 70% | a persistent interaction: part of the binding mode |
| 30–70% | an interaction that forms and breaks |
| Below 30% | an occasional contact |

JZ4 is mostly hydrophobic, so expect many hydrophobic contacts
(leucine, isoleucine, valine, phenylalanine, alanine) and at most one
hydrogen bond, through its OH group.

Want a different contact distance? Add `--cutoff 4.5`. Want more bars?
Add `--top 25`.

You can look inside the script to see how it works:

```bash
less ligand_contacts.py
```

> **👁 View it: the interacting residues**
> Download `start.pdb` (from Part 2) and open it in **PyMOL**. Type the
> residue numbers of your top contacts in place of the example numbers:
> ```
> hide everything
> show cartoon
> show sticks, resn JZ4
> select contacts, resi 84+87+91+99+111+118+121
> show sticks, contacts and not name N+C+O
> distance hbond, resn JZ4, contacts, mode=2
> zoom resn JZ4, 6
> ```
> `distance ... mode=2` draws dashed lines for hydrogen bonds PyMOL
> detects. This view is a good figure for a report.

---

## Part 9: Interaction energy

How strongly do the protein and ligand attract each other? GROMACS can
split the energy into contributions between groups of atoms. The
production run did not do this, so you **rerun** the trajectory: GROMACS
reads each saved frame and recalculates the energies, without moving
any atom.

This part uses the **full** system (with water), because the energies
depend on the environment.

### Step 1: a settings file with energy groups

```bash
cp md.mdp rerun.mdp
```

```bash
echo "energygrps = Protein JZ4" >> rerun.mdp
```

```bash
tail -n 2 rerun.mdp
```

`energygrps` asks GROMACS to report the energy between the protein and
JZ4 separately.

### Step 2: build the run input and rerun

```bash
gmx grompp -f rerun.mdp -c md.gro -p topol.top -n index.ndx -o rerun.tpr
```

```bash
gmx mdrun -s rerun.tpr -rerun md.xtc -deffnm rerun -nt 2
```

This takes a minute or two.

### Step 3: extract the protein-ligand energy

```bash
printf "Coul-SR:Protein-JZ4\nLJ-SR:Protein-JZ4\n\n" | gmx energy -f rerun.edr -o interaction_energy.xvg
```

| Term | Meaning |
|---|---|
| `Coul-SR:Protein-JZ4` | **electrostatic** energy between protein and ligand (charges; includes hydrogen bonds) |
| `LJ-SR:Protein-JZ4` | **Lennard-Jones** (van der Waals) energy: close packing and hydrophobic contacts |

```bash
python plot_xvg.py interaction_energy.xvg --avg 10 --ylabel "Energy (kJ/mol)"
```

**What you should see:** both terms are negative (attractive). For a
mostly hydrophobic ligand like JZ4, the **Lennard-Jones term is larger**
than the electrostatic term.

> **⚠ This is not a binding free energy.** The interaction energy
> ignores the cost of removing the ligand from water, the loss of
> flexibility on binding, and entropy. It is useful for comparing
> **which interactions** dominate, or comparing similar ligands in the
> same pocket, but it cannot be compared with experimental binding
> affinities. Methods that estimate binding free energy (MM/GBSA,
> MM/PBSA, free-energy perturbation) are introduced in Module 04.

> `-SR` means **short-range**: within the 1.0 nm cutoff. The long-range
> part of the electrostatics (calculated by PME) cannot be split between
> groups, so it is not included.

---

## Part 10: Clustering — what are the main ligand poses?

Clustering groups similar frames together, so you can describe a
trajectory with a few **representative structures** instead of
hundreds of frames.

The **GROMOS** method used here:

1. calculates the RMSD between every pair of frames,
2. finds the frame with the most neighbours within a cutoff (here
   0.1 nm = 1 Å), and calls it the centre of cluster 1,
3. removes that frame and its neighbours, and repeats for cluster 2,
   and so on.

We cluster on the **ligand position**. `md_fit.xtc` is already
superimposed on the protein backbone, so we switch off further fitting
(`-nofit`), and the ligand RMSD then measures how the ligand moves
**in the pocket**.

```bash
echo "JZ4 System" | gmx cluster -s complex.tpr -f md_fit.xtc -n complex.ndx -nofit -method gromos -cutoff 0.1 -g cluster.log -cl clusters.pdb -sz cluster_sizes.xvg -dist rmsd_distribution.xvg
```

The two answers: calculate the RMSD on **`JZ4`**, and write out the
whole complex (**`System`**) for each representative.

How many clusters, and how big is each?

```bash
grep -A 12 "cl. | #st" cluster.log
```

```bash
python plot_xvg.py cluster_sizes.xvg --bar --xlabel "Cluster" --ylabel "Number of frames"
```

**How to read it:** cluster 1 is always the largest. If it contains most
of the frames, the ligand has **one dominant pose**. Several clusters of
similar size mean the ligand **switches between poses**.

The cutoff controls how fine the clustering is. Try `-cutoff 0.05`
(more, smaller clusters) and `-cutoff 0.2` (fewer, larger ones).

> **👁 View it: the representative poses**
> Download `clusters.pdb`. It contains one **model** per cluster, the
> central structure of each.
> In **PyMOL**: open it, type `split_states clusters`, then
> `show sticks, resn JZ4`. Each cluster appears as a separate object in
> the list on the right; click them on and off to compare the poses.
> Is cluster 1 close to the crystal pose in `3HTB.pdb`?

---

## Part 11: PCA — the main motions of the protein

A protein with 160 residues moves in hundreds of directions at once.
**Principal component analysis (PCA)** finds the few directions
(**principal components**, PCs) that account for most of the motion.
This is also called **essential dynamics**.

### Step 1: calculate the motions

```bash
echo "C-alpha C-alpha" | gmx covar -s complex.tpr -f md_fit.xtc -n complex.ndx -o eigenvalues.xvg -v eigenvectors.trr
```

The answers: fit on **`C-alpha`**, and analyse the motion of the
**`C-alpha`** atoms. `gmx covar` builds a **covariance matrix** (how
each atom's motion is correlated with every other atom's) and finds its
**eigenvectors** (the directions of motion) and **eigenvalues** (how
much motion each direction carries).

### Step 2: how much motion do the first PCs explain?

```bash
python plot_xvg.py eigenvalues.xvg --bar --xmax 20 --xlabel "Principal component" --ylabel "Eigenvalue (nm²)" -o eigenvalues.png
```

`--xmax 20` shows only the first 20 components; there are as many as
three times the number of C-alpha atoms, but the rest are tiny.

**What you should see:** the first few bars are much larger than the
rest. Typically the first 2 or 3 PCs explain a large part of all motion.

### Step 3: project the trajectory onto PC1 and PC2

```bash
echo "C-alpha C-alpha" | gmx anaeig -s complex.tpr -f md_fit.xtc -n complex.ndx -v eigenvectors.trr -first 1 -last 2 -2d pc1_pc2.xvg
```

```bash
python plot_xvg.py pc1_pc2.xvg --scatter --xlabel "PC1 (nm)" --ylabel "PC2 (nm)" --title "Conformational space"
```

Each point is one frame, coloured from the start (purple) to the end
(yellow) of the simulation.

**How to read it:** points that form one cloud mean the protein explored
one region of conformations. Separate clouds mean it moved between
distinct states. A path that keeps moving in one direction means the
simulation had not finished exploring: **it was too short**, which is
expected at 100 ps.

### Step 4: watch the main motion

```bash
echo "C-alpha C-alpha" | gmx anaeig -s complex.tpr -f md_fit.xtc -n complex.ndx -v eigenvectors.trr -first 1 -last 1 -extr pc1_motion.pdb -nframes 20
```

This writes a short movie that moves the protein along PC1, from one
extreme to the other.

> **👁 View it: the dominant motion**
> Download `pc1_motion.pdb` and open it in **VMD** (File → New Molecule).
> In **Graphics → Representations**, set Drawing Method to **Tube**.
> Press **Play** and set the loop mode (the button next to the speed
> slider) to **Rock**. You are watching the single largest motion of
> the protein. Which regions move most? Compare with your RMSF plot.
>
> In **PyMOL**: open the file, type `show ribbon`, and press ▶.

---

## Summary: what each analysis told you

Fill in this table with your results. It is the core of a "Results"
section for any protein-ligand MD study.

| Question | Analysis | Your result |
|---|---|---|
| Equilibrated after... | RMSD (Part 3) | ___ ps |
| Protein stable? | Backbone RMSD, Rg | average ___ Å, Rg ___ Å |
| Most flexible regions | RMSF (Part 4) | residues ___ |
| Ligand pose stable? | Ligand RMSD, clusters | average ___ Å, cluster 1 = ___% |
| Ligand stays in pocket? | Pocket distance (Part 7) | ___ Å, steady / increasing |
| Ligand buried? | SASA (Part 6) | ___% buried |
| Key residues | Contacts (Part 8) | ___ |
| Hydrogen bonds | `hbonds.csv` | ___ with ___% occupancy |
| Dominant interaction | Interaction energy (Part 9) | LJ ___, Coulomb ___ kJ/mol |

---

## Common errors

| Error message contains | Usual cause and fix |
|---|---|
| `Atom index ... larger than the number of atoms` | The `.tpr` and the trajectory do not match. Use `complex.tpr` and `complex.ndx` with `md_fit.xtc`; use `md.tpr` and `index.ndx` with `md.xtc` or `md_noPBC.xtc`. |
| `Group ... not found` or the wrong group selected | Check the group names in your index file: `grep "^\[" complex.ndx` |
| `No such file or directory` | Check you are in `my-work/03-analysis` (`pwd`) and the file exists (`ls`). |
| `energy group ... not found` (Part 9) | Use `-n index.ndx` (the full-system index) in the `grompp` command, not `complex.ndx`. |
| `Invalid selection ... syntax error` | A typing error in a selection. Copy the command again exactly, including both kinds of quotes: `'...'` around the selection and `"..."` around group names. |
| The PNG did not appear | Click the refresh icon at the top of the file panel. |

---

## Exercises

1. How many frames does your trajectory contain, and how much time
   separates them? How would you change `md.mdp` to save a frame every
   0.5 ps?
2. Using `rmsd.png`, at what time does the backbone RMSD level off?
   Rerun the RMSF in Part 4 with `-b` set to that time. Does the RMSF
   change?
3. Which three residues have the highest RMSF? Are they in loops,
   helices or at the ends of the chain? (Use PyMOL.)
4. What fraction of JZ4's surface is buried by the protein?
5. List the five residues most often in contact with JZ4, with their
   occupancies. Are they hydrophobic or polar?
6. Does JZ4 form a hydrogen bond with the protein? With which residue
   and atom, and for what percentage of the time?
7. What are the average Coulomb and Lennard-Jones interaction energies?
   Which interaction dominates, and why does that fit the chemistry of JZ4?
8. How many clusters did you find with `-cutoff 0.1`? How many with
   `-cutoff 0.05`?
9. **Challenge:** in Part 3 you decided when the RMSD levels off. Repeat
   the contacts analysis (Part 8) and the clustering (Part 10) using only
   the frames after that time. For `gmx cluster`, add `-b` followed by
   the time in ps. For `ligand_contacts.py`, make a shorter trajectory
   first:
   ```bash
   echo "System" | gmx trjconv -s complex.tpr -f md_fit.xtc -o md_fit_eq.xtc -b 20
   ```
   (replace `20` with your time), then run the script on `md_fit_eq.xtc`.
   Do the results change when the start of the trajectory is left out?

---

## Files you should have at the end

```
complex.tpr  complex.ndx  md_fit.xtc  start.pdb
rmsd_backbone.xvg  rmsd_ligand.xvg  rmsd.png
rmsf.xvg  rmsf.png  bfactors.pdb
gyrate.xvg  gyrate.png
sasa_bound.xvg  sasa_free.xvg  ligand_sasa.png
pocket.ndx  pocket_distance.xvg  pocket_distance.png
contacts.png  contacts.csv  hbonds.csv
rerun.tpr  rerun.edr  interaction_energy.xvg  interaction_energy.png
cluster.log  clusters.pdb  cluster_sizes.xvg  cluster_sizes.png
eigenvalues.xvg  eigenvalues.png  pc1_pc2.xvg  pc1_pc2.png  pc1_motion.pdb
```

Remember to stop your codespace when you finish.
