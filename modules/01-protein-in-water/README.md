# Module 01: Simulating a Protein in Water

In this module you run a complete molecular dynamics simulation of
**hen egg-white lysozyme (PDB: 1AKI)** in water, from the crystal
structure to a 100 ps production run. Every simulation you run in your
research will follow these same steps.

**Time needed:** about 3 hours, including about 35 minutes of computing
time that runs in the background
**Before you start:** complete Module 00

---

## The workflow

```
1AKI.pdb  →  clean structure  →  topology  →  box  →  water  →  ions
                                                                  ↓
   production MD  ←  NPT equilibration  ←  NVT equilibration  ←  energy minimization
```

| Step | GROMACS tool | Output |
|---|---|---|
| 1. Clean the structure | `grep` | `1AKI_clean.pdb` |
| 2. Build the topology | `gmx pdb2gmx` | `processed.gro`, `topol.top`, `posre.itp` |
| 3. Define the box | `gmx editconf` | `newbox.gro` |
| 4. Add water | `gmx solvate` | `solv.gro` |
| 5. Add ions | `gmx grompp` + `gmx genion` | `solv_ions.gro` |
| 6. Energy minimization | `gmx grompp` + `gmx mdrun` | `em.gro` |
| 7. NVT equilibration (50 ps) | `gmx grompp` + `gmx mdrun` | `nvt.gro`, `nvt.cpt` |
| 8. NPT equilibration (50 ps) | `gmx grompp` + `gmx mdrun` | `npt.gro`, `npt.cpt` |
| 9. Production MD (100 ps) | `gmx grompp` + `gmx mdrun` | `md.xtc`, `md.edr`, `md.gro` |

GROMACS always works in two stages. **`gmx grompp`** combines the
structure (`.gro`), the topology (`.top`) and the settings (`.mdp`) into
one run input file (`.tpr`). **`gmx mdrun`** then runs the calculation
described in the `.tpr` file.

---

## Part 0: Get ready

Get the latest course files and make a folder for this module:

```bash
cd /workspaces/md-course && git pull
mkdir -p my-work/01-lysozyme
cd my-work/01-lysozyme
cp /workspaces/md-course/modules/01-protein-in-water/files/*.mdp .
ls
```

You should see five settings files: `ions.mdp`, `em.mdp`, `nvt.mdp`,
`npt.mdp` and `md.mdp`.

---

## Part 1: Get and clean the structure

```bash
wget https://files.rcsb.org/download/1AKI.pdb
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
```

Crystal waters are removed because we will surround the protein with
a complete box of water molecules in Part 4.

> **In your own research:** always inspect a structure before using it.
> Look for missing residues (`grep "REMARK 465"`), missing atoms
> (`grep "REMARK 470"`), ligands and ions you may need to keep, and
> multiple conformations of a side chain. 1AKI is a clean structure,
> which is why it is good for learning.

---

## Part 2: Build the topology

The **topology** describes every atom in the system: its mass, charge,
which atoms it is bonded to, and the force field parameters for every
bond, angle and dihedral.

```bash
gmx pdb2gmx -f 1AKI_clean.pdb -o processed.gro -water tip3p -ff amber99sb-ildn -ignh
```

| Option | Meaning |
|---|---|
| `-f 1AKI_clean.pdb` | input structure |
| `-o processed.gro` | output structure, in GROMACS format |
| `-ff amber99sb-ildn` | the **force field**: the set of equations and parameters describing the forces between atoms |
| `-water tip3p` | the **water model** to use with this force field |
| `-ignh` | ignore any hydrogens in the file and add them back correctly |

We use AMBER99SB-ILDN with TIP3P water because it is a widely used
protein force field, and it is compatible with the GAFF2 ligand
parameters you will use in Module 02.

Three files were created:

- `processed.gro`: the protein, now with hydrogen atoms
- `topol.top`: the topology
- `posre.itp`: **position restraints**, used later to hold the protein in place during equilibration

What is the total charge of the protein? Look at the last atom in the
topology:

```bash
grep "qtot" topol.top | tail -1
```

The number after `qtot` is the protein's net charge. Lysozyme should be
**+8**. Write this down; you will neutralize it in Part 5.

---

## Part 3: Define the simulation box

Simulations use **periodic boundary conditions**: the box is surrounded
by copies of itself in every direction, so the protein never sees an
edge or a vacuum.

```bash
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt dodecahedron
```

| Option | Meaning |
|---|---|
| `-c` | center the protein in the box |
| `-d 1.0` | keep at least 1.0 nm between the protein and the box edge |
| `-bt dodecahedron` | a rhombic dodecahedron box, which needs about 30% less water than a cube for the same protein |

---

## Part 4: Add water

```bash
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
```

`-cp` is the protein in its box, `-cs spc216.gro` is a small box of
pre-equilibrated water that GROMACS copies to fill the space, and
`-p topol.top` tells GROMACS to add the water molecules to the topology.

How many water molecules were added?

```bash
tail -3 topol.top
```

---

## Part 5: Add ions

The system has a net charge of +8, which must be neutralized. We also
add salt at a physiological concentration of 0.15 M.

`gmx genion` replaces water molecules with ions. It needs a `.tpr` file,
so first build one with `grompp`:

```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
```

> **Why `-maxwarn 1`?** `grompp` warns that the system has a net charge.
> That is true and expected: we are about to fix it. This is the **only**
> place in the course where you override a warning. In general, a
> GROMACS warning means something is wrong and must be fixed, not
> ignored.

Now add the ions:

```bash
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```

`echo "SOL"` answers genion's question about which molecules to replace
(water, called SOL). `-neutral` adds enough ions to make the total
charge zero, and `-conc 0.15` adds 0.15 M NaCl on top of that.

Check the result:

```bash
tail -4 topol.top
```

There should be 8 more CL than NA ions. Why?

---

## Part 6: Energy minimization

Adding water and ions can place some atoms too close together. Energy
minimization moves atoms slightly to remove these bad contacts, so the
simulation does not become unstable when it starts.

```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -nt 2
```

`-deffnm em` uses `em` as the base name for all files (`em.tpr` in;
`em.gro`, `em.log`, `em.edr`, `em.trr` out).

When it finishes, look for the line:

```
Steepest Descents converged to Fmax < 1000 in ... steps
```

The **potential energy** should be large and negative, and the
**maximum force** below 1000 kJ/mol/nm:

```bash
grep -E "converged|Potential Energy|Maximum force" em.log
```

Plot how the energy dropped during minimization:

```bash
echo "Potential" | gmx energy -f em.edr -o potential.xvg
```

---

## Part 7: NVT equilibration (50 ps)

Now the system is heated to 300 K at constant **N**umber of atoms,
**V**olume and **T**emperature (NVT). The protein is held in place with
position restraints (`define = -DPOSRES` in `nvt.mdp`), so the water
and ions can settle around it without the protein structure changing.

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
nohup gmx mdrun -deffnm nvt -nt 2 > nvt.out 2>&1 &
```

`-r em.gro` gives the reference positions for the restraints.

This takes about 5 to 10 minutes. While it runs, read through
`nvt.mdp` (`less nvt.mdp`) and find the settings for the thermostat,
the time step and the number of steps. Check progress with:

```bash
tail -n 15 nvt.log
```

When `jobs` shows nothing running, check the temperature:

```bash
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
```

The average should be close to 300 K.

---

## Part 8: NPT equilibration (50 ps)

Next, the pressure is brought to 1 bar at constant **N**umber, **P**ressure
and **T**emperature (NPT). The box size adjusts until the water
reaches the right density.

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
nohup gmx mdrun -deffnm npt -nt 2 > npt.out 2>&1 &
```

`-t nvt.cpt` passes the **checkpoint** from NVT, so the simulation
continues with the same atom velocities instead of starting again.

When it finishes, check pressure and density:

```bash
printf "Pressure\nDensity\n" | gmx energy -f npt.edr -o pressure_density.xvg
```

- **Density** should be close to 1000 kg/m³ (water with a protein and salt).
- **Pressure** will fluctuate by hundreds of bar. That is normal for a
  small system: look at the **average**, which should be near 1 bar
  within its large error.

---

## Part 9: Production MD (100 ps)

The system is now equilibrated. The production run has **no
restraints**: the protein moves freely. This is the trajectory you
analyze.

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
nohup gmx mdrun -deffnm md -nt 2 > md.out 2>&1 &
```

This takes about 10 to 20 minutes. Check the speed once it has started:

```bash
grep "Performance:" md.log    # appears when the run finishes
```

---

## Part 10: First look at the results

Detailed analysis is in Module 03. Here are three quick checks.

### Fix the periodic boundaries

Because of periodic boundaries, the protein can appear broken across
the edges of the box. Make it whole and centered:

```bash
echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center
```

The two answers are: center on `Protein`, write out the whole `System`.

### RMSD: how much has the structure changed?

```bash
echo "Backbone Backbone" | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ps
grep -v "^[#@]" rmsd.xvg | tail
```

The values are in nm. For a stable protein over 100 ps, expect small
values, around 0.05 to 0.15 nm.

### Radius of gyration: is the protein staying compact?

```bash
echo "Protein" | gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg
grep -v "^[#@]" gyrate.xvg | tail -3
```

The second column is the radius of gyration in nm. It should stay
nearly constant.

### Save a structure to look at

```bash
echo "Protein" | gmx trjconv -s md.tpr -f md_noPBC.xtc -o final_protein.pdb -dump 100
```

To download it, right-click `final_protein.pdb` in the file panel and
choose **Download**. Open it in PyMOL, ChimeraX or VMD on your computer,
or upload it to an online viewer such as molstar.org/viewer.

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

Run time depends on your machine. To estimate it for your own system:

```
hours needed = (simulated time in ns) ÷ (ns/day from md.log) × 24
```

To change the production run to 500 ps, for example:

```bash
nano md.mdp        # change nsteps to 250000, save with Ctrl+O, exit with Ctrl+X
grep nsteps md.mdp # check the change
```

Then run `grompp` and `mdrun` again as in Part 9.

> **Why only 100 ps here?** This course teaches the workflow, not
> research-quality sampling. 100 ps is long enough to see the protein
> move and to practice analysis. Real studies use tens to hundreds of
> nanoseconds, usually on GPUs or HPC clusters (Module 04).

### Extending a finished run

To add 400 ps to a finished 100 ps production run (500 ps total),
without starting again:

```bash
gmx convert-tpr -s md.tpr -extend 400 -o md_extended.tpr
nohup gmx mdrun -s md_extended.tpr -cpi md.cpt -deffnm md -nt 2 > md_ext.out 2>&1 &
```

`-extend` is in ps. New frames are added to the end of `md.xtc`.

### If your codespace stops during a run

If your codespace stops (idle timeout, closed browser, lost connection),
the running simulation stops too. Your files are kept, and GROMACS
writes a checkpoint (`.cpt`) every 15 minutes. Reopen the codespace,
go back to your folder, and continue from the checkpoint:

```bash
cd /workspaces/md-course/my-work/01-lysozyme
nohup gmx mdrun -deffnm md -cpi md.cpt -nt 2 > md_restart.out 2>&1 &
```

Replace `md` with `nvt` or `npt` if that was the stage that stopped.
To avoid stops, keep the codespace tab open during long runs, and set
your idle timeout to 240 minutes at github.com/settings/codespaces.

---

## If a step fails

Each step needs the files from the step before. If a step fails and you
cannot fix it, copy the reference files for that stage and continue.
Copy only the files listed for the stage where you got stuck, so you
do not overwrite your own work.

First, go to your folder and set a shortcut to the reference files:

```bash
cd /workspaces/md-course/my-work/01-lysozyme
R=/workspaces/md-course/modules/01-protein-in-water/reference
```

| If this part failed | Copy these files | Then continue from |
|---|---|---|
| Parts 1 to 5 (building the system) | `cp $R/{solv_ions.gro,topol.top,posre.itp} .` | Part 6 |
| Part 6 (energy minimization) | `cp $R/{em.gro,topol.top,posre.itp} .` | Part 7 |
| Part 7 (NVT) | `cp $R/{nvt.gro,nvt.cpt,topol.top,posre.itp} .` | Part 8 |
| Part 8 (NPT) | `cp $R/{npt.gro,npt.cpt,topol.top,posre.itp} .` | Part 9 |
| Part 9 (production) | `cp $R/{md.tpr,md.xtc,md.edr,md.gro} .` | Part 10 |

Then ask in the help channel about what went wrong, so you learn from it.

Common errors:

| Error message contains | Usual cause |
|---|---|
| `No such file or directory` | You are in the wrong folder (`pwd`), or a previous step did not produce its output (`ls`). |
| `Residue ... not found in residue topology database` | The structure contains a molecule the force field does not know (a ligand or non-standard residue). Remove it or parameterize it (Module 02). |
| `There were N warnings` / `Too many warnings` | Read the warning. Do not add `-maxwarn` to hide it. |
| `Atom ... in residue ... not found in rtp entry` | Atom names do not match the force field; often fixed with `-ignh`. |
| `LINCS warnings` or `water molecules can not be settled` | The system is unstable. Check the energy minimization converged. |

---

## Exercises

1. What is the net charge of 1AKI before adding ions? How many Na⁺ and
   Cl⁻ ions were added, and why are the numbers different?
2. How many water molecules are in your system? Approximately how many
   atoms in total? (Hint: `head -2 solv_ions.gro`)
3. How many steps did energy minimization take, and what was the final
   potential energy?
4. What were the average temperature (NVT) and average density (NPT)?
5. What was the performance of your production run in ns/day? Using
   the formula above, how long would a 100 ns simulation take on your
   codespace?
6. What is the final backbone RMSD after 100 ps?
7. Look at `nvt.mdp` and `md.mdp`. Which line makes the protein
   restrained in NVT but free in production? Why do we restrain it
   during equilibration?
8. **Challenge:** extend your production run to 200 ps using
   `gmx convert-tpr`. Does the RMSD keep increasing, or level off?

---

## Files you should have at the end

```
1AKI.pdb  1AKI_clean.pdb  processed.gro  topol.top  posre.itp
newbox.gro  solv.gro  solv_ions.gro
em.gro  em.edr  em.log
nvt.gro  nvt.cpt  nvt.edr  nvt.log
npt.gro  npt.cpt  npt.edr  npt.log
md.gro   md.cpt   md.edr   md.log   md.xtc   md.tpr
md_noPBC.xtc  rmsd.xvg  gyrate.xvg  final_protein.pdb
```

Keep this folder. Module 03 analyzes these files in detail.

Remember to stop your codespace when you finish.
