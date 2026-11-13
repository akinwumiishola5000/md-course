# Module 00: Working in the Linux Terminal

Every molecular dynamics simulation in this course is run by typing
commands in a terminal. In this module you learn those commands using
the same files you will use later: protein structures from the Protein
Data Bank and GROMACS input files. By the end, you will run your first
(very small) GROMACS simulation.

**Time needed:** about 2 hours
**What you need:** your codespace open, with the terminal at the bottom

---

## Part 1: The terminal

The terminal is the panel at the bottom of the screen showing:

```
(md-course) /workspaces/md-course $
```

This is the **prompt**. It tells you where you are (`/workspaces/md-course`)
and that it is waiting for a command. Type a command after the `$` and
press **Enter**.

Four keys that will save you a lot of time:

| Key | What it does |
|---|---|
| **Tab** | Completes file and folder names. Type the first few letters and press Tab. |
| **Up arrow** | Shows your previous commands, so you can run them again. |
| **Ctrl + C** | Stops the command that is running. |
| `clear` | Clears the screen (type it and press Enter). |

To copy from the terminal, select the text with your mouse. To paste,
use **Ctrl + Shift + V** (on a Mac, **Cmd + V**).

---

## Part 2: Moving around

Your files live in **folders** (also called **directories**). Try each
command below and look at what it prints.

```bash
pwd          # print working directory: where am I?
ls           # list what is in this folder
ls modules   # list what is inside the modules folder
```

Now create your own working folder and move into it:

```bash
mkdir my-work      # make a new folder called my-work
cd my-work         # change directory: move into my-work
pwd                # check where you are now
```

> **Important rule for the whole course:** do all your work inside
> `my-work`. Never edit files inside the `modules` folder. That way,
> new course materials can always be added to your codespace without
> problems.

Useful ways to move:

```bash
cd ..          # go up one folder
cd my-work     # go back into my-work
cd             # go to your home folder
cd /workspaces/md-course/my-work   # go to an exact location
```

Go back into `my-work` before continuing.

---

## Part 3: Downloading protein structures

Protein structures are stored in the **Protein Data Bank (PDB)** at
rcsb.org. Each structure has a four-character ID. You will download two:

- **1AKI**: hen egg-white lysozyme. You will simulate this protein in Module 01.
- **3HTB**: T4 lysozyme with a bound small molecule (ligand) called **JZ4**.
  You will simulate this protein-ligand complex in Module 02.

```bash
wget https://files.rcsb.org/download/1AKI.pdb
wget https://files.rcsb.org/download/3HTB.pdb
ls
```

You should now see `1AKI.pdb` and `3HTB.pdb`.

To see more detail about your files, including their sizes:

```bash
ls -lh
```

---

## Part 4: Looking inside a file

A PDB file is plain text. Each line starts with a **record type**
that tells you what the line contains.

```bash
head 1AKI.pdb          # show the first 10 lines
head -n 25 1AKI.pdb    # show the first 25 lines
tail 1AKI.pdb          # show the last 10 lines
wc -l 1AKI.pdb         # count the lines in the file
```

To scroll through the whole file, use `less`:

```bash
less 1AKI.pdb
```

Inside `less`: use the **arrow keys** or **Space** to scroll, type
`/ATOM` and press Enter to jump to the first atom, and press **q** to quit.

The record types you will meet most often:

| Record | Contains |
|---|---|
| `TITLE` | The name of the structure |
| `REMARK   2` | The resolution of the crystal structure |
| `ATOM` | One atom of the protein (standard amino acids) |
| `HETATM` | One atom of anything else: water, ions, ligands |
| `TER` | The end of a protein chain |

---

## Part 5: Searching with grep

`grep` finds every line that contains a piece of text.

```bash
grep TITLE 1AKI.pdb
grep "REMARK   2 RESOLUTION" 1AKI.pdb
```

The `^` symbol means "at the start of the line". This lets you select
lines by record type:

```bash
grep "^ATOM" 1AKI.pdb | head      # the first protein atoms
grep "^HETATM" 1AKI.pdb | head    # the first non-protein atoms
```

Add `-c` to **count** matching lines instead of printing them:

```bash
grep -c "^ATOM" 1AKI.pdb      # how many protein atoms?
grep -c "^HETATM" 1AKI.pdb    # how many non-protein atoms?
```

### Pipes: joining commands together

The `|` symbol (called a **pipe**) sends the output of one command into
the next. This is one of the most useful ideas in Linux.

Count the water molecules (water is called `HOH` in PDB files; each
crystal water is one oxygen atom):

```bash
grep "^HETATM" 1AKI.pdb | grep -c HOH
```

Count the amino acid residues. Every amino acid has exactly one alpha
carbon, named `CA`, so counting `CA` atoms counts residues:

```bash
grep "^ATOM" 1AKI.pdb | grep " CA " | wc -l
```

The spaces around `CA` matter: without them, `grep` would also match
other text containing "CA".

---

## Part 6: Cutting out columns

PDB files use fixed columns. The residue name is always in characters
18 to 20, and the chain ID is always character 22. The `cut` command
extracts those columns.

Which chains does the protein have?

```bash
grep "^ATOM" 1AKI.pdb | cut -c22 | sort -u
```

Here `sort -u` sorts the results and keeps only unique values.

Which amino acids does lysozyme contain, and how many of each?

```bash
grep "^ATOM" 1AKI.pdb | grep " CA " | cut -c18-20 | sort | uniq -c
```

To see the most common amino acids first, add one more sort:

```bash
grep "^ATOM" 1AKI.pdb | grep " CA " | cut -c18-20 | sort | uniq -c | sort -rn
```

Read this command from left to right: take the protein atoms, keep the
alpha carbons, cut out the residue names, sort them, count each one,
then sort by count from largest to smallest.

---

## Part 7: Saving output to a file

The `>` symbol saves output to a file instead of printing it.

This is a real preparation step for molecular dynamics: removing the
crystal waters before setting up the simulation. The `-v` option tells
`grep` to keep every line that does **not** match:

```bash
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
grep -c HOH 1AKI_clean.pdb      # should print 0
```

> **Careful:** `>` replaces the file if it already exists. Use `>>`
> to add to the end of a file instead.

### Extracting a ligand from 3HTB

First, find out what non-protein molecules 3HTB contains:

```bash
grep "^HETATM" 3HTB.pdb | cut -c18-20 | sort | uniq -c
```

Now save the ligand JZ4 into its own file. You will do exactly this in
Module 02 before parameterizing the ligand.

```bash
grep "^HETATM" 3HTB.pdb | grep JZ4 > JZ4.pdb
cat JZ4.pdb
```

`cat` prints a whole file at once, which is fine for short files.

Why `grep "^HETATM"` first? The name JZ4 also appears in header lines
that describe the ligand. Selecting `HETATM` lines first makes sure you
only get the ligand's atoms.

---

## Part 8: Copying, moving and deleting

```bash
mkdir structures                      # make a folder
cp 1AKI.pdb structures/               # copy a file into it
mv 3HTB.pdb structures/               # move a file into it
ls structures
mv structures/1AKI.pdb structures/lysozyme.pdb   # rename a file
rm structures/lysozyme.pdb            # delete a file
```

> **Warning:** `rm` deletes immediately. There is no recycle bin.
> Check the file name carefully before pressing Enter.

---

## Part 9: Editing a file with nano

GROMACS reads its simulation settings from `.mdp` files. You will
often need to change a setting, so you need a text editor that works
in the terminal.

Copy the practice files from the module folder into a new folder:

```bash
mkdir water-test
cp /workspaces/md-course/modules/00-linux-basics/files/* water-test/
cd water-test
ls
```

Look at the MD settings file:

```bash
cat md.mdp
```

Each line has the form `setting = value`. Everything after `;` is a
comment that GROMACS ignores. Find `nsteps` and `dt`. The total
simulated time is `nsteps × dt`.

Open the file in nano:

```bash
nano md.mdp
```

Use the **arrow keys** to move to the `nsteps` line and change `5000`
to `10000`. Then:

- **Ctrl + O**, then **Enter** to save
- **Ctrl + X** to exit

Check your change:

```bash
grep nsteps md.mdp
```

How many picoseconds will this simulation now run? (`dt` is in ps.)

---

## Part 10: Your first GROMACS simulation

You will now run a short simulation of a small box of water. In Module
01 you will learn what each step means in detail; for now, focus on
how commands are run.

Make sure you are inside `water-test`, then run these one at a time:

```bash
# 1. Fill a 3 x 3 x 3 nm box with water molecules
gmx solvate -cs spc216.gro -box 3 3 3 -o water.gro -p water.top

# 2. Prepare the energy minimization
gmx grompp -f em.mdp -c water.gro -p water.top -o em.tpr

# 3. Run the energy minimization
gmx mdrun -deffnm em -nt 2
```

GROMACS prints a lot of text. That is normal. Look for the line saying
that steepest descents **converged**.

How many water molecules were added? `gmx solvate` updated the topology:

```bash
tail -3 water.top
```

### Running a job in the background

Real simulations run for hours. You do not want to keep the terminal
busy the whole time, and you want the job to keep running if your
browser disconnects. For that, use `nohup` and `&`:

```bash
# 4. Prepare the MD run
gmx grompp -f md.mdp -c em.gro -p water.top -o md.tpr

# 5. Run it in the background
nohup gmx mdrun -deffnm md -nt 2 > mdrun.out 2>&1 &
```

- `nohup` keeps the job running even if you close the terminal.
- `> mdrun.out 2>&1` saves all of the screen output to `mdrun.out`.
- `&` at the end runs the job in the background and gives you your prompt back.
- `-nt 2` uses 2 CPU threads, matching the 2 cores of your codespace.

Check that it is running, and watch the log file grow:

```bash
jobs
tail -f md.log
```

Press **Ctrl + C** to stop watching. This stops `tail`, **not** the
simulation. Run `jobs` again to check whether the simulation has
finished.

---

## Part 11: Reading the results

How fast was the simulation?

```bash
grep "Performance:" md.log
```

The first number is **ns/day**: how many nanoseconds of simulation this
machine can produce in one day. This number tells you how long a
simulation will take before you start it.

What was the average temperature? `gmx energy` normally asks you to
choose a term interactively. You can send it the answer with `echo`
and a pipe instead:

```bash
echo "Temperature" | gmx energy -f md.edr -o temperature.xvg
```

Look at the summary at the end of the output. The average should be
close to 300 K, the temperature set by `ref_t` in `md.mdp`.

The temperature over time is saved in `temperature.xvg`. Lines starting
with `#` or `@` are formatting information, so remove them to see the data:

```bash
grep -v "^[#@]" temperature.xvg | head
```

The first column is time (ps) and the second is temperature (K).

---

## Exercises

Write down your answers. Your instructor will go through them at the
start of the next session.

1. What is the resolution of 1AKI, in Ångströms?
2. How many amino acid residues does 1AKI have?
3. How many crystal water molecules are in 1AKI?
4. Which amino acid occurs most often in 1AKI, and how many times?
5. How many chains does 3HTB have?
6. Apart from water and JZ4, what other non-protein molecules are in 3HTB?
7. How many atoms are in `JZ4.pdb`? Why are there no hydrogen atoms?
8. Copy `md.mdp` to `md_310.mdp`, change the temperature to 310 K
   (two settings must change: find them), and run the simulation again
   using `-deffnm md310`. What average temperature do you get?
9. Using the ns/day value from Part 11, how long would 1 ns of this
   water box take on your codespace?

**Bonus:** write a single command that counts how many `HOH` molecules
are in 3HTB and saves the number to a file called `water_count.txt`.

---

## Command summary

| Command | What it does |
|---|---|
| `pwd` | Show the current folder |
| `ls`, `ls -lh` | List files (with sizes) |
| `cd folder`, `cd ..` | Move into a folder, move up one |
| `mkdir name` | Make a folder |
| `wget URL` | Download a file |
| `head`, `tail`, `cat`, `less` | Look at a file |
| `wc -l` | Count lines |
| `grep text file` | Find lines containing text |
| `grep -c`, `grep -v` | Count matches, keep non-matching lines |
| `cut -c18-20` | Extract columns 18 to 20 |
| `sort`, `sort -u`, `uniq -c` | Sort, keep unique, count repeats |
| `a \| b` | Send the output of `a` into `b` |
| `> file`, `>> file` | Save output (replace / add to the end) |
| `cp`, `mv`, `rm` | Copy, move or rename, delete |
| `nano file` | Edit a file |
| `nohup command &` | Run a job in the background |
| `jobs` | List background jobs |
| `tail -f file` | Watch a file as it grows |

---

## Before you leave

Stop your codespace so it does not use your free hours: go to
https://github.com/codespaces, click the three dots next to your
codespace, and choose **Stop codespace**.

Your files in `my-work` are kept. Next session, reopen the same codespace.
