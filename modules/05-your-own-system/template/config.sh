# =====================================================================
#  config.sh - settings for run_md.sh
#  The Insilico Lab MD Course, Module 05
#
#  This is the ONLY file you need to edit. Each setting is explained.
#  Lines starting with # are comments. Do not put spaces around "=".
# =====================================================================


# ---------------------------------------------------------------------
# 1. YOUR STRUCTURE
# ---------------------------------------------------------------------

# The PDB file of your protein (or protein-ligand complex), in this folder.
# If the file is not here and you give a 4-character PDB ID instead
# (for example INPUT_PDB="3HTB"), it is downloaded from the RCSB PDB.
INPUT_PDB="3HTB.pdb"


# ---------------------------------------------------------------------
# 2. YOUR LIGAND (leave LIGAND_RESNAME empty for a protein on its own)
# ---------------------------------------------------------------------

# The 3-character residue name of the ligand.
#  - Ligand inside the PDB file: use its name in the HETATM lines
#    (find it with:  bash run_md.sh check).
#  - Ligand in a separate file (LIGAND_FILE below): any 3 capital
#    letters, for example LIG.
#  - Protein only: leave empty:  LIGAND_RESNAME=""
LIGAND_RESNAME="JZ4"

# A separate ligand file, for example a docking pose (.sdf, .mol2 or
# .pdb). It must already be positioned in the binding site of
# INPUT_PDB. Leave empty to take the ligand from INPUT_PDB instead.
LIGAND_FILE=""

# The ligand's net charge at the pH below. "auto" works it out from the
# protonated structure; CHECK it with  bash run_md.sh check.
# Write a number (for example 0, 1 or -1) to set it yourself.
LIGAND_CHARGE="auto"

# The pH used to add hydrogens to the ligand (charged amines,
# carboxylates...). 7.4 is physiological pH.
LIGAND_PH="7.4"


# ---------------------------------------------------------------------
# 3. PREPARING THE PROTEIN
# ---------------------------------------------------------------------

# Repair the protein with PDBFixer before building the topology:
# adds missing side-chain atoms and converts modified amino acids
# (for example selenomethionine MSE -> MET). Recommended: yes
FIX_PROTEIN="yes"

# What to do with GAPS: residues missing in the MIDDLE of a chain (run
# "bash run_md.sh check" to see if your structure has any).
#   "stop"  stop before building anything, so you can decide (default).
#   "fill"  rebuild the missing residues with PDBFixer. Needs SEQRES
#           records in the PDB file (files from the PDB have them), and
#           only gaps up to MAX_FILL_LENGTH residues. The rebuilt loops are
#           ROUGH models: use only for short gaps far from your binding site.
#   "split" add nothing: each gap becomes a chain break, with two free
#           chain ends. For gaps far from the binding site.
# For gaps near the binding site, model the missing residues with a
# dedicated tool first (see 'Gaps in the chain' in the README).
# Missing residues at the START or END of a chain are never added.
MISSING_RESIDUES="stop"
MAX_FILL_LENGTH="5"

# ---------------------------------------------------------------------
# 4. FORCE FIELD AND SOLVENT
# ---------------------------------------------------------------------

# Protein force field and water model. Keep the defaults if you have a
# ligand: the ligand is parameterized with GAFF2, which is compatible
# with AMBER force fields.
FORCEFIELD="amber99sb-ildn"
WATER_MODEL="tip3p"

# Minimum distance between the protein and the box edge (nm).
# 1.0 is standard. Larger values add more water (slower).
BOX_DISTANCE="1.0"

# Box shape: dodecahedron (least water, fastest), cubic or triclinic.
BOX_TYPE="dodecahedron"

# Salt concentration in mol/L. 0.15 is physiological.
SALT_CONC="0.15"


# ---------------------------------------------------------------------
# 5. SIMULATION LENGTH  (all times in picoseconds, 1 ns = 1000 ps)
# ---------------------------------------------------------------------
#
# The defaults are SHORT, so a course codespace finishes in under an
# hour for a protein of ~150-300 residues. They teach the workflow but
# are NOT long enough for research. For a research study, typical
# values are:
#
#   Setting           Course (default)   Research
#   NVT_PS                    50           100 - 500
#   NPT_PS                    50           100 - 1000
#   MD_PS                    100           100000 (100 ns) or more
#   SAVE_EVERY_PS              1           10 - 100
#
# Before production starts, run_md.sh measures the speed of YOUR system
# on YOUR computer and prints how long MD_PS will take.

NVT_PS="50"           # NVT equilibration: heat to TEMPERATURE
NPT_PS="50"           # NPT equilibration: adjust the density
MD_PS="100"           # production run: the trajectory you analyse

# How often a frame is saved to the trajectory (ps). Aim for about
# 1000-5000 frames in total: MD_PS / SAVE_EVERY_PS.
# Example: 100 ns (100000 ps) with SAVE_EVERY_PS=50 gives 2000 frames.
SAVE_EVERY_PS="1"

# Temperature in kelvin. 300 K is standard; 310 K is body temperature.
TEMPERATURE="300"


# ---------------------------------------------------------------------
# 6. COMPUTER
# ---------------------------------------------------------------------

# Number of CPU threads. "auto" uses all CPU cores of this computer.
THREADS="auto"

# The GROMACS command. "gmx" in the codespace and on most computers.
# Some HPC clusters only provide an MPI version called "gmx_mpi".
GMX_COMMAND="gmx"

# Extra options for gmx mdrun. Leave empty in a codespace.
# On a computer with an NVIDIA GPU and a GPU build of GROMACS, use:
#   MDRUN_OPTIONS="-nb gpu -pme gpu -bonded gpu -update gpu"
MDRUN_OPTIONS=""

# Folder containing plot_xvg.py, ligand_contacts.py and plot_mmgbsa.py.
# The default works in the course codespace. Elsewhere, copy those
# three files into this folder and set TOOLS_DIR="."
TOOLS_DIR="/workspaces/md-course/modules/tools"


# ---------------------------------------------------------------------
# 7. BINDING FREE ENERGY (only used when there is a ligand)
# ---------------------------------------------------------------------

# Run MM/GBSA (gmx_MMPBSA, Module 04) after the analysis: yes or no.
RUN_MMGBSA="yes"

# Number of frames to use for MM/GBSA, spread evenly over production.
# More frames = more reliable but slower (~10 s per frame per 2 cores).
MMGBSA_FRAMES="20"
