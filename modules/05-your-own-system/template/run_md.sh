#!/usr/bin/env bash
# =====================================================================
#  run_md.sh - molecular dynamics of a protein or protein-ligand
#  complex, from a PDB file to analysed results.
#  The Insilico Lab MD Course, Module 05
#
#  Edit config.sh, then run one of:
#
#    bash run_md.sh check        inspect your structure, change nothing
#    bash run_md.sh all          run every stage below, in order
#
#  or run the stages one at a time:
#
#    bash run_md.sh prepare      clean, topology, box, water, ions  (Modules 01-02)
#    bash run_md.sh minimize     energy minimization                (Module 01, Part 6)
#    bash run_md.sh equilibrate  NVT and NPT                        (Module 01, Parts 7-8)
#    bash run_md.sh production   production MD                      (Module 01, Part 9)
#    bash run_md.sh analyse      RMSD, RMSF, Rg, contacts, plots    (Module 03)
#    bash run_md.sh mmgbsa       binding free energy                (Module 04)
#
#  Long runs: start it in the background so it continues if your
#  browser disconnects, and watch the progress:
#
#    nohup bash run_md.sh all > run.log 2>&1 &
#    tail -f run.log
#
#  "all" skips stages that have already finished, and production
#  continues from its checkpoint, so if anything stops, simply run
#  "bash run_md.sh all" again.
#
#  Folders created:
#    run/      all GROMACS files (structures, topology, trajectories)
#    results/  plots (.png), tables (.csv) and summary.txt
#    logs/     the full output of every command, for troubleshooting
# =====================================================================

set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE" || exit 1
if [ ! -f config.sh ]; then
  echo "ERROR: config.sh not found in $HERE"; exit 1
fi
# shellcheck disable=SC1091
source config.sh

# Defaults for settings added in later versions of config.sh, so an older
# config.sh still works.
MISSING_RESIDUES=${MISSING_RESIDUES:-stop}
MAX_FILL_LENGTH=${MAX_FILL_LENGTH:-5}
GMX_COMMAND=${GMX_COMMAND:-gmx}

RUN="$HERE/run"
# Overwrite files when a stage is rerun, instead of keeping GROMACS
# backups (#file.1#). A finished earlier run is kept by "prepare" instead.
export GMX_MAXBACKUP=-1
RES="$HERE/results"
LOGS="$HERE/logs"
DT=0.002                                   # time step in ps (2 fs), as in every module
STAGES=(prepare minimize equilibrate production analyse mmgbsa)

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
say()  { echo; echo "=== $(date +%H:%M:%S)  $*"; }
info() { echo "    $*"; }
ok()   { echo "    [OK]   $*"; }
warn() { echo "    [WARNING] $*"; }
die()  { echo; echo "    [STOPPED] $*"; echo; exit 1; }

# Explain the most common GROMACS and AmberTools errors in plain words.
explain_error() {
  local log=$1
  echo "    ---------------------------------------------------------"
  echo "    What this usually means:"
  if grep -q "not found in residue topology database" "$log"; then
    echo "    The protein contains a residue the force field does not know"
    echo "    (a ligand, ion or modified amino acid). Run 'bash run_md.sh check'"
    echo "    and look at the HETATM list. Set LIGAND_RESNAME for your ligand,"
    echo "    and keep FIX_PROTEIN=\"yes\" to convert modified amino acids."
  elif grep -qiE "is missing in residue|not found in rtp entry|atom .* missing" "$log"; then
    echo "    Atoms are missing from some residues. Set FIX_PROTEIN=\"yes\" in"
    echo "    config.sh so PDBFixer can add them, then run prepare again."
  elif grep -q "Atomtype .* not found" "$log"; then
    echo "    The ligand parameters are not included correctly in the topology."
    echo "    Run 'bash run_md.sh prepare' again from the start."
  elif grep -q "does not match topology" "$log"; then
    echo "    The structure and the topology have different numbers of atoms."
    echo "    Run 'bash run_md.sh prepare' again from the start."
  elif grep -qE "LINCS|can not be settled|Segmentation fault|blowing up" "$log"; then
    echo "    The simulation became unstable. Usually: atoms too close"
    echo "    together (check that minimization converged), a wrong ligand"
    echo "    structure or charge, or a gap in the protein chain (see 'check')."
  elif grep -qiE "sqm|antechamber|charge" "$log" && [[ "$log" == *acpype* ]]; then
    echo "    Ligand parameterization failed. The most common cause is a wrong"
    echo "    net charge or missing hydrogens. Check the formula and charge"
    echo "    printed by 'bash run_md.sh check', and set LIGAND_CHARGE."
  elif grep -q "interaction residues" "$log"; then
    echo "    gmx_MMPBSA found no protein residue within 4 A of the ligand: the"
    echo "    ligand has probably left the binding pocket. Look at the ligand RMSD"
    echo "    and distance plots in results/ before trusting any binding energy."
  elif grep -q "Too many warnings" "$log"; then
    echo "    grompp found problems. Read the WARNING lines in the log above."
  else
    echo "    Read the error message above. The full output is in:"
  fi
  echo "    Full log: $log"
  echo "    ---------------------------------------------------------"
}

# Run a command, saving its output to logs/NAME.log. Stop with an
# explanation if it fails.
run() {
  local name=$1; shift
  local log="$LOGS/$name.log"
  if ! "$@" > "$log" 2>&1; then
    echo
    echo "    [FAILED] $name"
    echo "    Last lines of the output:"
    grep -v "^\s*$" "$log" | tail -n 12 | sed 's/^/      | /'
    explain_error "$log"
    exit 1
  fi
}

# Same, for GROMACS tools that ask a question: $1 = answer(s), $2 = name
run_answer() {
  local answer=$1 name=$2; shift 2
  local log="$LOGS/$name.log"
  if ! printf "%b\n" "$answer" | "$@" > "$log" 2>&1; then
    echo; echo "    [FAILED] $name"
    grep -v "^\s*$" "$log" | tail -n 12 | sed 's/^/      | /'
    explain_error "$log"; exit 1
  fi
}

nsteps()  { python3 -c "print(int(round($1 / $DT)))"; }

# Each stage writes its own part of results/summary.txt, so rerunning a
# stage replaces its lines instead of adding them again.
summary() {           # summary STAGE  (reads the text from stdin)
  cat > "$RES/.summary_$1"
  : > "$RES/summary.txt"
  for s in "${STAGES[@]}"; do
    [ -f "$RES/.summary_$s" ] && cat "$RES/.summary_$s" >> "$RES/summary.txt"
  done
}
mark_done()    { touch "$RUN/.done_$1"; }
is_done()      { [ -f "$RUN/.done_$1" ]; }
need()         { is_done "$1" || die "Run the '$1' stage first:  bash run_md.sh $1"; }
clear_after()  {                      # rerunning a stage invalidates the later ones
  local found=0
  for s in "${STAGES[@]}"; do
    [ $found = 1 ] && rm -f "$RUN/.done_$s"
    [ "$s" = "$1" ] && found=1
  done
  return 0
}
has_ligand() { [ -n "${LIGAND_RESNAME}" ]; }
avg_energy() { printf "%s\n" "$1" | gmx energy -f "$2" -o "$RUN/tmp_energy.xvg" 2>/dev/null | awk -v t="$1" '$1==t {print $2}'; rm -f "$RUN/tmp_energy.xvg"; }

plot() {          # plot_xvg.py if available
  if [ -n "$PLOT" ]; then python "$PLOT" "$@" >> "$LOGS/plots.log" 2>&1 || warn "plot failed: $*"; fi
}

[ "$THREADS" = "auto" ] && THREADS=$(nproc)

# GROMACS command: "gmx" normally. If GMX_COMMAND is different (for
# example gmx_mpi on a cluster), every "gmx" below runs that instead.
GMX_COMMAND=${GMX_COMMAND:-gmx}
if [ "$GMX_COMMAND" != "gmx" ]; then gmx() { "$GMX_COMMAND" "$@"; }; fi
command -v "$GMX_COMMAND" > /dev/null || die "GROMACS command '$GMX_COMMAND' not found. Check GMX_COMMAND in config.sh, or load GROMACS (module load ...)."
NT_FLAG="-nt"                                   # threads option of mdrun
[[ "$GMX_COMMAND" == *mpi* ]] && NT_FLAG="-ntomp" # MPI builds use -ntomp
PLOT=""; CONTACTS=""; MMPLOT=""
for d in "$TOOLS_DIR" "$HERE"; do
  [ -z "$PLOT" ]     && [ -f "$d/plot_xvg.py" ]       && PLOT="$d/plot_xvg.py"
  [ -z "$CONTACTS" ] && [ -f "$d/ligand_contacts.py" ] && CONTACTS="$d/ligand_contacts.py"
  [ -z "$MMPLOT" ]   && [ -f "$d/plot_mmgbsa.py" ]     && MMPLOT="$d/plot_mmgbsa.py"
done

if has_ligand && ! [[ "$LIGAND_RESNAME" =~ ^[A-Z0-9]{1,3}$ ]]; then
  die "LIGAND_RESNAME must be 1 to 3 capital letters or digits (you wrote \"$LIGAND_RESNAME\")."
fi

# ---------------------------------------------------------------------
# Inspect the input structure and the ligand.
# Used by "check", and at the start of "prepare".
# Writes clean files into the folder given as $1.
# ---------------------------------------------------------------------
inspect() {
  local dir=$1
  mkdir -p "$dir"

  # --- Find the input file (download it if a PDB ID was given)
  local pdb="$INPUT_PDB"
  if [ ! -f "$pdb" ]; then
    local id=${INPUT_PDB%.pdb}
    if [[ "$id" =~ ^[0-9A-Za-z]{4}$ ]]; then
      info "Downloading $id from the RCSB PDB..."
      if ! wget -q "https://files.rcsb.org/download/${id^^}.pdb" -O "${id^^}.pdb.part"; then
        rm -f "${id^^}.pdb.part"     # never leave an empty or partial file behind
        die "Could not download $id. Check the ID and your internet connection, or put the file in this folder yourself."
      fi
      mv "${id^^}.pdb.part" "${id^^}.pdb"
      pdb="${id^^}.pdb"
    else
      die "INPUT_PDB file \"$INPUT_PDB\" not found in $HERE"
    fi
  fi
  info "Input structure: $pdb"
  grep "^TITLE" "$pdb" | head -n 2 | sed 's/^TITLE */    title: /'
  grep "REMARK   2 RESOLUTION" "$pdb" | sed 's/.*RESOLUTION\. */    resolution: /'

  # --- Keep only the first model (NMR files contain many) and the
  #     first alternate location (altLoc A) of each atom. Remove water.
  awk '/^ENDMDL/ {exit}
       /^(ATOM|HETATM)/ { alt=substr($0,17,1); if (alt!=" " && alt!="A") next;
                          if (substr($0,18,3)=="HOH" || substr($0,18,3)=="WAT") next;
                          print substr($0,1,16) " " substr($0,18) ; next }
       /^(TER|SEQRES)/ {print}' "$pdb" > "$dir/clean.pdb"
  local nalt; nalt=$(awk '/^(ATOM|HETATM)/ {a=substr($0,17,1); if (a!=" ") n++} END {print n+0}' "$pdb")
  local nmodel; nmodel=$(grep -c "^MODEL" "$pdb" || true)
  [ "$nmodel" -gt 1 ] && info "The file has $nmodel models (NMR?); only the first is used."
  [ "$nalt" -gt 0 ] && info "Alternate locations found ($nalt atoms): position A is kept."

  # --- Chains, residues and non-protein molecules
  info "Protein chains (residues with a C-alpha atom):"
  grep "^ATOM" "$dir/clean.pdb" | awk 'substr($0,13,4)==" CA " {print substr($0,22,1)}' \
     | sort | uniq -c | awk '{printf "      chain %s: %d residues\n", ($2==""?"(blank)":$2), $1}'
  info "Other molecules (HETATM residue names, with atom counts):"
  if grep -q "^HETATM" "$dir/clean.pdb"; then
    grep "^HETATM" "$dir/clean.pdb" | cut -c18-20 | sort | uniq -c | sed 's/^/      /'
    info "Everything in this list is REMOVED, except your ligand (LIGAND_RESNAME)"
    info "and modified amino acids such as MSE, which FIX_PROTEIN=yes converts."
    info "Metal ions or cofactors that matter for your protein need extra care:"
    info "see 'Ions, metals and cofactors' in the Module 05 README."
  else
    info "  none (water is removed automatically)"
  fi

  # --- Gaps in the chains (missing residues)
  local gaps; gaps=$(find_gaps "$dir/clean.pdb")
  if [ -z "$gaps" ]; then
    ok "No gaps found in the protein chains."
  else
    warn "Gaps in the protein chain (residues missing from the structure):"
    echo "$gaps" | awk -F'\t' '{printf "      chain %s: between %s and %s (%s A apart; bonded residues are 3.8 A)\n", $1, $2, $3, $4}'
    local seqres="no"; grep -q "^SEQRES" "$dir/clean.pdb" && seqres="yes"
    info "The file has SEQRES records (the full sequence): $seqres"
    case "$MISSING_RESIDUES" in
      fill)  info "MISSING_RESIDUES=\"fill\": gaps up to $MAX_FILL_LENGTH residues will be rebuilt with PDBFixer."
             [ "$seqres" = "no" ] && warn "Without SEQRES records the missing residues are unknown, so they cannot be rebuilt." ;;
      split) info "MISSING_RESIDUES=\"split\": each gap will become a chain break (no residues are added)." ;;
      *)     info "MISSING_RESIDUES=\"stop\": 'prepare' will stop. Choose how to handle the gaps:"
             info "see 'Gaps in the chain' in the Module 05 README." ;;
    esac
  fi

  # --- The ligand
  if ! has_ligand; then
    info "No ligand (LIGAND_RESNAME is empty): protein-only simulation."
    return 0
  fi

  if [ -n "$LIGAND_FILE" ]; then
    [ -f "$LIGAND_FILE" ] || die "LIGAND_FILE \"$LIGAND_FILE\" not found in $HERE"
    info "Ligand from file: $LIGAND_FILE"
    local fmt=${LIGAND_FILE##*.}
    obabel -i"$fmt" "$LIGAND_FILE" -l 1 -opdb -O "$dir/ligand_raw.pdb" > /dev/null 2>&1 \
      || die "Open Babel could not read $LIGAND_FILE"
    # remove any hydrogens: they are added again at LIGAND_PH below
    awk '/^(ATOM|HETATM)/ && $NF!="H" {print}' "$dir/ligand_raw.pdb" > "$dir/ligand.pdb"
  else
    grep "^HETATM" "$dir/clean.pdb" | awk -v r="$LIGAND_RESNAME" 'substr($0,18,3)==sprintf("%3s",r) || substr($0,18,3)==sprintf("%-3s",r)' > "$dir/ligand_all.pdb"
    [ -s "$dir/ligand_all.pdb" ] || die "No HETATM residue called $LIGAND_RESNAME in $INPUT_PDB. Check the list of other molecules above, and LIGAND_RESNAME in config.sh."
    # if there are several copies (one per chain, for example), keep the first
    local first; first=$(head -n 1 "$dir/ligand_all.pdb" | cut -c22-27)
    awk -v f="$first" 'substr($0,22,6)==f' "$dir/ligand_all.pdb" > "$dir/ligand.pdb"
    local ncopies; ncopies=$(cut -c22-27 "$dir/ligand_all.pdb" | sort -u | wc -l)
    [ "$ncopies" -gt 1 ] && warn "$ncopies copies of $LIGAND_RESNAME found; only the first (chain/residue $first) is simulated."
    # remove hydrogens if the file has any: they are added again at LIGAND_PH
    awk '{el=substr($0,77,2); gsub(/ /,"",el); if (el!="H") print}' "$dir/ligand.pdb" > "$dir/ligand_noH.pdb"
    mv "$dir/ligand_noH.pdb" "$dir/ligand.pdb"
  fi
  info "Ligand heavy atoms: $(grep -c . "$dir/ligand.pdb")"

  # Add hydrogens for LIGAND_PH (charged groups are protonated or
  # deprotonated as at that pH), and save as MOL2 with bonds.
  obabel "$dir/ligand.pdb" -O "$dir/ligand_H.mol2" -p "$LIGAND_PH" > /dev/null 2>&1 \
    || die "Open Babel could not add hydrogens to the ligand."
  # Give the ligand the residue name LIGAND_RESNAME. Ligands read from SDF
  # files are otherwise called "UNL", and ACPYPE keeps that name.
  python3 - "$dir/ligand_H.mol2" "$LIGAND_RESNAME" << 'PY'
import sys
path, name = sys.argv[1], sys.argv[2]
out, section = [], ""
for line in open(path):
    if line.startswith("@<TRIPOS>"):
        section = line.strip()
    elif section == "@<TRIPOS>ATOM" and line.strip():
        f = line.split()
        if len(f) >= 9:
            line = (f"{int(f[0]):7d} {f[1]:<8s} {float(f[2]):10.4f} {float(f[3]):10.4f} "
                    f"{float(f[4]):10.4f} {f[5]:<6s} {1:3d}  {name:<8s} {float(f[8]):9.4f}\n")
    elif section == "@<TRIPOS>SUBSTRUCTURE" and line.strip():
        f = line.split()
        f[1] = name
        line = " ".join(f) + "\n"
    out.append(line)
open(path, "w").writelines(out)
PY
  local props
  props=$(python3 - "$dir/ligand_H.mol2" << 'PY'
import sys
from openbabel import pybel
m = next(pybel.readfile("mol2", sys.argv[1]))
print(m.formula, m.charge, sum(1 for a in m.atoms if a.atomicnum == 1))
PY
)
  local formula=${props%% *}; props=${props#* }
  local auto_charge=${props%% *}; local nh=${props#* }
  info "Ligand formula with hydrogens at pH $LIGAND_PH: $formula  ($nh hydrogens)"
  info "Ligand net charge at pH $LIGAND_PH (from the structure): $auto_charge"
  if [ "$LIGAND_CHARGE" = "auto" ]; then
    LIG_Q=$auto_charge
  else
    LIG_Q=$LIGAND_CHARGE
    [ "$LIG_Q" != "$auto_charge" ] && warn "LIGAND_CHARGE=$LIG_Q in config.sh, but the protonated structure has charge $auto_charge. Make sure $LIG_Q is correct."
  fi
  info "Charge that will be used: $LIG_Q"
  [ "$nh" -eq 0 ] && warn "No hydrogens were added. Check that the ligand file has correct bonds."
  echo "$LIG_Q" > "$dir/ligand_charge.txt"
}

# ---------------------------------------------------------------------
# Gaps: consecutive residues of a chain whose C-alpha atoms are more than
# 4.2 A apart (bonded residues are always ~3.8 A apart). Modified amino
# acids stored as HETATM (MSE...) are included, so they are not gaps.
# Prints one line per gap:  chain <TAB> residue before <TAB> residue after
# <TAB> distance <TAB> residue number of the residue after (for splitting)
# ---------------------------------------------------------------------
find_gaps() {
  python3 - "$1" "${LIGAND_RESNAME}" << 'PY'
import sys, math
path, ligand = sys.argv[1], sys.argv[2]
ca = {}
for l in open(path):
    if l.startswith(("ATOM", "HETATM")) and l[12:16] == " CA " and l[17:20].strip() != ligand:
        ca.setdefault(l[21], []).append((l[22:27].strip(), l[17:20].strip(),
            (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
for ch, res in ca.items():
    for a, b in zip(res, res[1:]):
        d = math.dist(a[2], b[2])
        if d > 4.2:
            print(f"{ch}\t{a[1]}{a[0]}\t{b[1]}{b[0]}\t{d:.1f}\t{b[0]}")
PY
}

# Insert TER lines at the gaps, so pdb2gmx treats each piece of the chain
# as a separate chain with its own two ends (MISSING_RESIDUES="split").
split_at_gaps() {
  python3 - "$1" "$2" << 'PY'
import sys
pdb, gaps = sys.argv[1], sys.argv[2]
starts = {(g.split("\t")[0], g.split("\t")[4]) for g in gaps.splitlines() if g.strip()}
out, prev = [], None
for l in open(pdb):
    if l.startswith(("ATOM", "HETATM")):
        key = (l[21], l[22:27].strip())
        if key != prev and key in starts:
            out.append("TER\n")
        prev = key
    out.append(l)
open(pdb, "w").writelines(out)
PY
}

# cis peptide bonds (omega dihedral within 30 degrees of 0) that are NOT
# before a proline. These are rare in real proteins; new ones after
# minimization usually come from badly rebuilt residues.
find_cis() {   # find_cis FILE (.pdb or .gro)
  python3 - "$1" << 'PY'
import sys, math
path = sys.argv[1]
atoms, order = {}, []
gro = path.endswith(".gro")
lines = open(path).read().splitlines()
if gro:
    n = int(lines[1]); recs = lines[2:2 + n]
    for a in recs:
        key = (int(a[0:5]), a[5:10].strip()); name = a[10:15].strip()
        if key not in atoms: atoms[key] = {}; order.append(key)
        atoms[key][name] = [float(a[k:k + 8]) * 10 for k in (20, 28, 36)]
else:
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            key = (l[21] + l[22:27].strip(), l[17:20].strip()); name = l[12:16].strip()
            if key not in atoms: atoms[key] = {}; order.append(key)
            atoms[key][name] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
def sub(a, b): return [x - y for x, y in zip(a, b)]
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = sub(p0, p1), sub(p2, p1), sub(p3, p2)
    n1, n2 = cross(b0, b1), cross(b1, b2)
    m = cross(n1, [x / math.sqrt(dot(b1, b1)) for x in b1])
    return math.degrees(math.atan2(dot(m, n2), dot(n1, n2)))
for r1, r2 in zip(order, order[1:]):
    a, b = atoms[r1], atoms[r2]
    if not all(k in a for k in ("CA", "C")) or not all(k in b for k in ("N", "CA")):
        continue
    if math.dist(a["C"], b["N"]) > 2.0:     # not bonded: a chain end
        continue
    if r2[1].endswith("PRO"):                 # cis before proline is normal
        continue
    if abs(dihedral(a["CA"], a["C"], b["N"], b["CA"])) < 30:
        num = lambda r: str(r[0])[1:] if not gro else str(r[0])
        print(f"{num(r1)}-{num(r2)} {r1[1][-3:]}{num(r1)}-{r2[1][-3:]}{num(r2)}")
PY
}

# Quality check after minimization and equilibration: report non-proline
# cis peptide bonds that were not in the input structure.
cis_check() {   # cis_check FILE.gro   (run inside run/)
  local cis_in cis_now new_cis
  cis_in=$(find_cis input/clean.pdb)
  cis_now=$(find_cis "$1")
  # compare by residue numbers (GROMACS renames some residues, e.g. HIS -> HIE)
  new_cis=$(awk 'NR==FNR {seen[$1]; next} NF && !($1 in seen) {print $2}' <(echo "$cis_in") <(echo "$cis_now"))
  if [ -n "$new_cis" ]; then
    warn "cis peptide bond(s) not present in the input structure:"
    echo "$new_cis" | sed 's/^/      /'
    info "Non-proline cis peptides are rare in real proteins. They usually come from"
    info "rebuilt residues (MISSING_RESIDUES=\"fill\"): that region is unreliable. If it"
    info "matters for your question, model the missing residues properly (README, 'Gaps')."
  else
    ok "no new cis peptide bonds"
  fi
}

# =====================================================================
# Stage: check
# =====================================================================
stage_check() {
  say "CHECK: inspecting $INPUT_PDB (no files are changed)"
  local tmp; tmp=$(mktemp -d)
  inspect "$tmp"
  rm -rf "$tmp"
  echo
  info "Settings: force field $FORCEFIELD/$WATER_MODEL, $BOX_TYPE box (${BOX_DISTANCE} nm),"
  info "${SALT_CONC} M NaCl, ${TEMPERATURE} K, NVT ${NVT_PS} ps, NPT ${NPT_PS} ps, production ${MD_PS} ps,"
  info "a frame every ${SAVE_EVERY_PS} ps ($(python3 -c "print(int($MD_PS/$SAVE_EVERY_PS))") frames), $THREADS CPU threads."
  echo
  info "If everything above looks right, run:  bash run_md.sh all"
}

# =====================================================================
# Stage: prepare  (Module 01 Parts 1-5, Module 02 Parts 2-8)
# =====================================================================
stage_prepare() {
  say "PREPARE: cleaning the structure and building the system"
  # Starting again from scratch: keep a backup of an earlier run, so a
  # finished simulation is never lost, then start with an empty folder.
  if [ -f "$RUN/md.xtc" ]; then
    local backup="run_backup_$(date +%Y%m%d_%H%M%S)"
    mv "$RUN" "$HERE/$backup"; info "Earlier run moved to $backup/"
  else
    rm -rf "$RUN"
  fi
  mkdir -p "$RUN" "$RES" "$LOGS"
  inspect "$RUN/input"
  cd "$RUN" || exit 1
  cp "$HERE/config.sh" "$RUN/config_used.sh"     # a record of the settings

  # --- 1. The protein
  say "Protein: building the clean protein structure"
  if [ "$FIX_PROTEIN" = "yes" ]; then
    # PDBFixer: replace modified amino acids, remove everything that is
    # not protein, optionally rebuild missing residues, add missing atoms.
    run pdbfixer python3 - input/clean.pdb protein.pdb "$MISSING_RESIDUES" "$MAX_FILL_LENGTH" << 'PY'
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fixer = PDBFixer(filename=sys.argv[1])
fixer.findNonstandardResidues()
for res, std in fixer.nonstandardResidues:
    print(f"replacing modified residue {res.name}{res.id} -> {std}")
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingResidues()
mode, max_len = sys.argv[3], int(sys.argv[4])
chains = list(fixer.topology.chains())
keep = {}
for (c, i), names in fixer.missingResidues.items():
    n_res = len(list(chains[c].residues()))
    if i == 0 or i == n_res:
        # Missing residues at the start or end of a chain are usually
        # disordered tails: they are never added.
        print(f"not adding {len(names)} missing residues at an end of chain {chains[c].id}")
    elif mode != "fill":
        continue
    elif len(names) > max_len:
        print(f"NOT rebuilding a gap of {len(names)} residues in chain {chains[c].id}: longer than MAX_FILL_LENGTH={max_len}")
    else:
        keep[(c, i)] = names
        print(f"rebuilding missing residues {' '.join(names)} in chain {chains[c].id}")
fixer.missingResidues = keep
fixer.findMissingAtoms()
n_missing = sum(len(a) for a in fixer.missingAtoms.values())
print(f"adding {n_missing} missing heavy atoms")
fixer.addMissingAtoms()
with open(sys.argv[2], "w") as out:
    PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
PY
    grep -E "replacing|adding|rebuilding|NOT rebuilding" "$LOGS/pdbfixer.log" | sed 's/^/    PDBFixer: /'
  else
    grep "^ATOM" input/clean.pdb > protein.pdb
    echo "END" >> protein.pdb
  fi
  grep -q "^ATOM" protein.pdb || die "No protein atoms found in $INPUT_PDB."

  # Gaps left in the protein? GROMACS would join the two ends of every gap
  # with one long bond, then minimization would drag them together and
  # distort the protein, with no error. So never continue silently.
  local gaps; gaps=$(find_gaps protein.pdb)
  if [ -n "$gaps" ]; then
    local list; list=$(echo "$gaps" | awk -F'\t' '{printf "      chain %s: between %s and %s (%s A)\n", $1, $2, $3, $4}')
    if [ "$MISSING_RESIDUES" = "split" ]; then
      split_at_gaps protein.pdb "$gaps"
      warn "Gaps kept as chain breaks (MISSING_RESIDUES=\"split\"). Each gap now has two free chain ends:"
      echo "$list"
    else
      echo; echo "    [STOPPED] The protein has gaps (missing residues):"; echo "$list"
      if [ "$MISSING_RESIDUES" = "fill" ]; then
        if [ "$FIX_PROTEIN" != "yes" ]; then
          info "MISSING_RESIDUES=\"fill\" needs FIX_PROTEIN=\"yes\"."
        elif ! grep -q "^SEQRES" input/clean.pdb; then
          info "They could not be rebuilt: $INPUT_PDB has no SEQRES records, so the missing"
          info "residues are unknown. Model them with another tool, or use MISSING_RESIDUES=\"split\"."
        else
          info "They were not rebuilt: longer than MAX_FILL_LENGTH=$MAX_FILL_LENGTH residues (see above)."
          info "Model long gaps with a dedicated tool, or use MISSING_RESIDUES=\"split\"."
        fi
      else
        info "Set MISSING_RESIDUES in config.sh to \"fill\" or \"split\", or model the missing"
        info "residues with another tool first. See 'Gaps in the chain' in the Module 05 README."
      fi
      echo; exit 1
    fi
  fi

  # --- 2. Protein topology (Module 01, Part 2)
  run pdb2gmx gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top \
      -ff "$FORCEFIELD" -water "$WATER_MODEL" -ignh
  ok "protein topology built ($(sed -n 2p protein.gro | tr -d ' ') atoms with hydrogens)"

  # --- 3. Ligand topology (Module 02, Parts 3-7)
  if has_ligand; then
    LIG_Q=$(cat input/ligand_charge.txt)
    say "Ligand: GAFF2 parameters with AM1-BCC charges (ACPYPE), net charge $LIG_Q"
    info "This takes from under a minute (small ligands) to 10+ minutes (large ones)."
    cp input/ligand_H.mol2 ligand_H.mol2
    rm -rf "${LIGAND_RESNAME}.acpype"
    run acpype acpype -i ligand_H.mol2 -c bcc -n "$LIG_Q" -a gaff2 -b "$LIGAND_RESNAME"
    local acp="${LIGAND_RESNAME}.acpype"
    cp "$acp/${LIGAND_RESNAME}_GMX.itp" "$acp/${LIGAND_RESNAME}_GMX.gro" "$acp/posre_${LIGAND_RESNAME}.itp" .
    local qtot; qtot=$(grep "qtot" "${LIGAND_RESNAME}_GMX.itp" | tail -n 1 | awk '{print $NF}')
    ok "ligand topology built (total charge $qtot)"

    # Add the ligand to topol.top, as in Module 02 Part 6:
    #  - include its parameters directly after the force field line
    #  - add one ligand molecule after the protein molecule(s)
    run topology python3 - topol.top "$LIGAND_RESNAME" "$FORCEFIELD" << 'PY'
import sys
top, lig, ff = sys.argv[1], sys.argv[2], sys.argv[3]
text = open(top).read()
inc = f'#include "{ff}.ff/forcefield.itp"\n'
if inc not in text:
    sys.exit(f"force field include line not found in {top}")
add = (inc + f'\n; Ligand parameters\n#include "{lig}_GMX.itp"\n\n'
       f'; Ligand position restraints\n#ifdef POSRES\n#include "posre_{lig}.itp"\n#endif\n')
text = text.replace(inc, add, 1).rstrip("\n") + f"\n{lig:<20}1\n"
open(top, "w").write(text)
print("ligand added to topology")
PY

    # Join protein and ligand coordinates (Module 02, Part 7)
    head -n -1 protein.gro > complex.gro
    tail -n +3 "${LIGAND_RESNAME}_GMX.gro" | head -n -1 >> complex.gro
    tail -n 1 protein.gro >> complex.gro
    sed -i "2s/.*/ $(( $(wc -l < complex.gro) - 3 ))/" complex.gro
    START=complex.gro
  else
    START=protein.gro
  fi

  # --- 4. Box, water and ions (Module 01, Parts 3-5)
  say "Solvating: $BOX_TYPE box, water, ions"
  run editconf gmx editconf -f "$START" -o newbox.gro -c -d "$BOX_DISTANCE" -bt "$BOX_TYPE"
  run solvate  gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
  write_mdps
  # grompp warns about the net charge before ions are added; allow that
  # one warning only, and only if it is the net-charge warning.
  if ! gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr > "$LOGS/grompp_ions.log" 2>&1; then
    if grep -q "net charge\|non-zero total charge" "$LOGS/grompp_ions.log" && \
       [ "$(grep -c '^WARNING' "$LOGS/grompp_ions.log")" -eq 1 ]; then
      run grompp_ions gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
    else
      echo "    [FAILED] grompp_ions"; grep -A4 "^WARNING\|^ERROR\|Fatal" "$LOGS/grompp_ions.log" | sed 's/^/      | /'
      explain_error "$LOGS/grompp_ions.log"; exit 1
    fi
  fi
  run_answer "SOL" genion gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
      -pname NA -nname CL -neutral -conc "$SALT_CONC"

  # --- 5. Index groups for the thermostat: Solute (protein + ligand)
  #        and Solvent (everything else), as in Module 02 Part 10.
  local n; n=$(echo q | gmx make_ndx -f solv_ions.gro -o tmp.ndx 2>/dev/null | grep -c "^ *[0-9]* [A-Za-z_-].*: *[0-9]* atoms")
  local solute='"Protein"'
  has_ligand && solute="\"Protein\" | r ${LIGAND_RESNAME}"
  run_answer "${solute}\nname $n Solute\n! $n\nname $((n+1)) Solvent\nq" make_ndx \
      gmx make_ndx -f solv_ions.gro -o index.ndx
  rm -f tmp.ndx
  grep -q "^\[ Solute \]" index.ndx && grep -q "^\[ Solvent \]" index.ndx \
      || die "Could not create the Solute and Solvent index groups (see logs/make_ndx.log)."
  # Safety check: Solute must contain exactly the protein plus the ligand
  local n_solute n_protein n_lig=0
  n_solute=$(awk '/^\[/ {in_group = ($2 == "Solute"); next} in_group {n += NF} END {print n+0}' index.ndx)
  n_protein=$(sed -n 2p protein.gro | tr -d ' ')
  has_ligand && n_lig=$(( $(wc -l < "${LIGAND_RESNAME}_GMX.gro") - 3 ))
  if has_ligand && ! grep -q "^\[ ${LIGAND_RESNAME} \]" index.ndx; then
    die "The ligand residue $LIGAND_RESNAME is not in the system (see logs/make_ndx.log)."
  fi
  [ "$n_solute" = "$((n_protein + n_lig))" ] \
    || die "The Solute group has $n_solute atoms, expected $n_protein (protein) + $n_lig (ligand). See logs/make_ndx.log."
  ok "Solute group: $n_protein protein atoms + $n_lig ligand atoms"

  # --- Report
  local natoms; natoms=$(sed -n 2p solv_ions.gro | tr -d ' ')
  local nwat; nwat=$(awk '$1=="SOL" {print $2}' topol.top | tail -n 1)
  local nna; nna=$(awk '$1=="NA" {print $2}' topol.top | tail -n 1)
  local ncl; ncl=$(awk '$1=="CL" {print $2}' topol.top | tail -n 1)
  ok "system built: $natoms atoms, ${nwat:-0} waters, ${nna:-0} Na+, ${ncl:-0} Cl-"
  rm -f "$RES"/.summary_*
  {
    echo "System"
    echo "  input: $INPUT_PDB    ligand: ${LIGAND_RESNAME:-none}"
    echo "  molecules:"; sed -n '/\[ molecules \]/,$p' topol.top | tail -n +3 | sed 's/^/    /'
    echo "  atoms: $natoms"
  } | summary prepare
  cd "$HERE"; mark_done prepare; clear_after prepare
}

# Write the .mdp files from the settings in config.sh. These are the
# same settings as in Modules 01 and 02; only lengths and temperature vary.
write_mdps() {
  local save; save=$(nsteps "$SAVE_EVERY_PS")
  local common="cutoff-scheme  = Verlet
nstlist        = 10
coulombtype    = PME
rcoulomb       = 1.0
rvdw           = 1.0
pbc            = xyz"
  local dyn="constraint_algorithm = lincs
constraints    = h-bonds
pme_order      = 4
fourierspacing = 0.16
DispCorr       = EnerPres
tcoupl         = V-rescale
tc-grps        = Solute Solvent
tau_t          = 0.1     0.1
ref_t          = $TEMPERATURE     $TEMPERATURE"

  printf "; ions.mdp - only used to prepare gmx genion\nintegrator     = steep\nnsteps         = 0\n%s\n" "$common" > ions.mdp
  printf "; em.mdp - energy minimization\ndefine         = -DFLEXIBLE\nintegrator     = steep\nemtol          = 1000.0\nemstep         = 0.01\nnsteps         = 50000\n%s\n" "$common" > em.mdp
  cat > nvt.mdp << EOF
; nvt.mdp - NVT equilibration: ${NVT_PS} ps at ${TEMPERATURE} K, solute restrained
define         = -DPOSRES
integrator     = md
dt             = $DT
nsteps         = $(nsteps "$NVT_PS")
nstxout-compressed = $(nsteps 1)
nstenergy      = $(nsteps 1)
nstlog         = $(nsteps 1)
continuation   = no
$common
$dyn
pcoupl         = no
gen_vel        = yes
gen_temp       = $TEMPERATURE
gen_seed       = -1
EOF
  cat > npt.mdp << EOF
; npt.mdp - NPT equilibration: ${NPT_PS} ps at 1 bar, solute restrained
define         = -DPOSRES
integrator     = md
dt             = $DT
nsteps         = $(nsteps "$NPT_PS")
nstxout-compressed = $(nsteps 1)
nstenergy      = $(nsteps 1)
nstlog         = $(nsteps 1)
continuation   = yes
$common
$dyn
pcoupl         = C-rescale
pcoupltype     = isotropic
tau_p          = 2.0
ref_p          = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
gen_vel        = no
EOF
  cat > md.mdp << EOF
; md.mdp - production: ${MD_PS} ps, no restraints, a frame every ${SAVE_EVERY_PS} ps
integrator     = md
dt             = $DT
nsteps         = $(nsteps "$MD_PS")
nstxout-compressed = $save
compressed-x-grps  = System
nstenergy      = $save
nstlog         = $(( save * 10 ))
continuation   = yes
$common
$dyn
pcoupl         = C-rescale
pcoupltype     = isotropic
tau_p          = 2.0
ref_p          = 1.0
compressibility = 4.5e-5
gen_vel        = no
EOF
}

# =====================================================================
# Stage: minimize  (Module 01, Part 6)
# =====================================================================
stage_minimize() {
  need prepare
  say "MINIMIZE: energy minimization"
  cd "$RUN" || exit 1
  rm -f nvt.tpr nvt.gro nvt.cpt nvt.edr nvt.log nvt.xtc npt.tpr npt.gro npt.cpt npt.edr npt.log npt.xtc \
        md.tpr md.gro md.cpt md.edr md.log md.xtc md_noPBC.xtc md_fit.xtc
  run grompp_em gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
  run mdrun_em  gmx mdrun -deffnm em $NT_FLAG "$THREADS"
  local epot fmax
  epot=$(grep "^Potential Energy" em.log | awk '{print $4}')
  fmax=$(grep "^Maximum force" em.log | awk '{print $4}')
  python3 -c "import sys; sys.exit(0 if float('$fmax') < 1000 else 1)" \
    && ok "converged: potential energy $epot kJ/mol, maximum force $fmax" \
    || {
      python3 -c "import sys; sys.exit(0 if float('$fmax') < 1e5 else 1)" \
        || die "Maximum force $fmax is far too high: atoms overlap. Check gaps and the ligand ('bash run_md.sh check')."
      warn "maximum force $fmax is above 1000. Usually still fine; watch equilibration closely."
    }
  cis_check em.gro
  echo "Potential" | gmx energy -f em.edr -o "$RES/potential.xvg" > "$LOGS/energy_em.log" 2>&1
  plot "$RES/potential.xvg" --xlabel "Minimization step" --ylabel "Potential energy (kJ/mol)" --title "Energy minimization" -o "$RES/1_potential.png"
  echo "Minimization: potential energy $epot kJ/mol, maximum force $fmax" | summary minimize
  cd "$HERE"; mark_done minimize; clear_after minimize
}

# =====================================================================
# Stage: equilibrate  (Module 01, Parts 7-8)
# =====================================================================
stage_equilibrate() {
  need minimize
  cd "$RUN" || exit 1
  rm -f md.tpr md.gro md.cpt md.edr md.log md.xtc md_noPBC.xtc md_fit.xtc
  say "EQUILIBRATE: NVT, ${NVT_PS} ps at ${TEMPERATURE} K"
  run grompp_nvt gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
  run mdrun_nvt  gmx mdrun -deffnm nvt $NT_FLAG "$THREADS" $MDRUN_OPTIONS
  local t; t=$(avg_energy Temperature nvt.edr)
  ok "average temperature $t K (target $TEMPERATURE)"
  python3 -c "import sys; sys.exit(0 if abs(float('$t')-$TEMPERATURE) < 10 else 1)" \
    || warn "the temperature is more than 10 K from the target. Consider a longer NVT_PS."
  echo "Temperature" | gmx energy -f nvt.edr -o "$RES/temperature.xvg" > "$LOGS/energy_nvt.log" 2>&1
  plot "$RES/temperature.xvg" --avg 10 -o "$RES/2_temperature.png"

  say "EQUILIBRATE: NPT, ${NPT_PS} ps at 1 bar"
  run grompp_npt gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
  run mdrun_npt  gmx mdrun -deffnm npt $NT_FLAG "$THREADS" $MDRUN_OPTIONS
  local p d; p=$(avg_energy Pressure npt.edr); d=$(avg_energy Density npt.edr)
  ok "average pressure $p bar (target 1; large fluctuations are normal)"
  ok "average density $d kg/m3"
  python3 -c "import sys; sys.exit(0 if 950 < float('$d') < 1100 else 1)" \
    || warn "the density is unusual (expected 1000-1050). Consider a longer NPT_PS."
  echo "Pressure" | gmx energy -f npt.edr -o "$RES/pressure.xvg" > "$LOGS/energy_npt_p.log" 2>&1
  echo "Density"  | gmx energy -f npt.edr -o "$RES/density.xvg"  > "$LOGS/energy_npt_d.log" 2>&1
  plot "$RES/pressure.xvg" --avg 10 -o "$RES/3_pressure.png"
  plot "$RES/density.xvg"  --avg 10 -o "$RES/4_density.png"
  cis_check npt.gro
  echo "Equilibration: temperature $t K, pressure $p bar, density $d kg/m3" | summary equilibrate
  cd "$HERE"; mark_done equilibrate; clear_after equilibrate
}

# =====================================================================
# Stage: production  (Module 01, Part 9)
# =====================================================================
stage_production() {
  need equilibrate
  cd "$RUN" || exit 1
  say "PRODUCTION: ${MD_PS} ps, frame every ${SAVE_EVERY_PS} ps"

  # Estimate the run time from the measured speed of the NPT run
  local speed; speed=$(grep "Performance:" npt.log | awk '{print $2}')
  if [ -n "$speed" ]; then
    python3 - "$speed" "$MD_PS" << 'PY'
import sys
nsday, ps = float(sys.argv[1]), float(sys.argv[2])
hours = ps / 1000 / nsday * 24
print(f"    Speed of this system on this computer: {nsday:.1f} ns/day")
print(f"    Estimated production time: {hours*60:.0f} minutes" if hours < 2 else
      f"    Estimated production time: {hours:.1f} hours ({hours/24:.1f} days)")
if hours > 3.5:
    print("    [WARNING] Longer than a codespace session can stay idle. Keep the")
    print("    tab open, or run production on an HPC cluster (Module 04, Part B).")
    print("    If it stops, run 'bash run_md.sh all' again: it continues from the checkpoint.")
PY
  fi

  if [ -f md.cpt ] && [ -f md.tpr ]; then
    info "Checkpoint found: continuing the interrupted run."
    run mdrun_md gmx mdrun -deffnm md -cpi md.cpt $NT_FLAG "$THREADS" $MDRUN_OPTIONS
  else
    run grompp_md gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
    run mdrun_md  gmx mdrun -deffnm md $NT_FLAG "$THREADS" $MDRUN_OPTIONS
  fi
  grep -q "Finished mdrun" md.log || die "Production did not finish. Run 'bash run_md.sh all' again to continue."
  speed=$(grep "Performance:" md.log | tail -n 1 | awk '{print $2}')
  ok "production finished (${speed} ns/day)"
  echo "Production: ${MD_PS} ps at ${speed} ns/day" | summary production
  cd "$HERE"; mark_done production; clear_after production
}

# =====================================================================
# Stage: analyse  (Module 03)
# =====================================================================
stage_analyse() {
  need production
  cd "$RUN" || exit 1
  say "ANALYSE: preparing the trajectory"
  # Make the solute whole and centred (Module 03, Part 2). "-pbc cluster"
  # keeps all solute molecules (every protein chain and the ligand) in the
  # SAME periodic image. The simpler "-pbc mol" of Modules 02-03 keeps each
  # molecule whole, but can leave a surface ligand or a second chain in a
  # neighbouring image: the analysis and MM/GBSA would then see it 20+ A away.
  run_answer "Solute\nSolute\nSystem" trjconv_pbc gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_noPBC.xtc -pbc cluster -center
  # Keep only the solute (protein + ligand): smaller and faster to analyse
  run_answer "Solute" solute_traj gmx trjconv -s md.tpr -f md_noPBC.xtc -n index.ndx -o solute_noPBC.xtc
  run_answer "Solute" convert_tpr gmx convert-tpr -s md.tpr -n index.ndx -o solute.tpr
  run_answer "q" make_ndx_solute gmx make_ndx -f solute.tpr -o solute.ndx
  # Reference structure = first frame of the corrected trajectory. The
  # structure stored inside md.tpr is NOT used as reference, because in it
  # a ligand or second chain may still sit in a neighbouring periodic image.
  run_answer "System" reference gmx trjconv -s solute.tpr -f solute_noPBC.xtc -o reference.pdb -dump 0
  # Fit every frame on the backbone of the reference: removes tumbling
  run_answer "Backbone\nSystem" trjconv_fit gmx trjconv -s reference.pdb -f solute_noPBC.xtc -n solute.ndx -o md_fit.xtc -fit rot+trans
  cp reference.pdb "$RES/start.pdb"
  cp md_fit.xtc "$RES/trajectory.xtc"

  say "ANALYSE: RMSD, RMSF, radius of gyration"
  run_answer "Backbone\nBackbone" rms_bb gmx rms -s reference.pdb -f md_fit.xtc -n solute.ndx -o "$RES/rmsd_backbone.xvg" -tu ps
  run_answer "C-alpha" rmsf gmx rmsf -s reference.pdb -f md_fit.xtc -n solute.ndx -o "$RES/rmsf.xvg" -res -oq "$RES/rmsf_bfactors.pdb"
  run_answer "Protein" gyrate gmx gyrate -s solute.tpr -f md_fit.xtc -n solute.ndx -o "$RES/gyrate.xvg"
  plot "$RES/rmsf.xvg"   --unit A --xlabel "Residue number" -o "$RES/6_rmsf.png"
  plot "$RES/gyrate.xvg" --unit A --cols 1 -o "$RES/7_gyrate.png"

  mean() { grep -v '^[#@]' "$1" | awk -v f="${2:-1}" -v c="${3:-2}" '{s+=$c; n++} END {printf "%.2f", f*s/n}'; }
  local sum="Analysis
  backbone RMSD (mean): $(mean "$RES/rmsd_backbone.xvg" 10) A
  radius of gyration (mean): $(mean "$RES/gyrate.xvg" 10) A"

  if has_ligand; then
    say "ANALYSE: ligand RMSD, distance and contacts"
    run_answer "Backbone\n${LIGAND_RESNAME}" rms_lig gmx rms -s reference.pdb -f md_fit.xtc -n solute.ndx -o "$RES/rmsd_ligand.xvg" -tu ps
    plot "$RES/rmsd_backbone.xvg" "$RES/rmsd_ligand.xvg" --unit A --avg 10 -o "$RES/5_rmsd.png" --title "Backbone and ligand RMSD"
    run_answer "Protein\n${LIGAND_RESNAME}" mindist gmx mindist -s solute.tpr -f md_fit.xtc -n solute.ndx -od "$RES/mindist.xvg" -tu ps
    plot "$RES/mindist.xvg" --unit A -o "$RES/8_mindist.png"
    sum+="
  ligand RMSD (mean): $(mean "$RES/rmsd_ligand.xvg" 10) A
  protein-ligand minimum distance (mean): $(mean "$RES/mindist.xvg" 10) A"
    if [ -n "$CONTACTS" ]; then
      ( cd "$RES" && python "$CONTACTS" "$RUN/solute.tpr" "$RUN/md_fit.xtc" --ligand "$LIGAND_RESNAME" ) > "$LOGS/contacts.log" 2>&1 \
        && mv "$RES/contacts.png" "$RES/9_contacts.png" 2>/dev/null
      grep -E "[0-9] %$" "$LOGS/contacts.log" | head -n 8 | sed 's/^ */    /'
      sum+="
  most frequent contacts (% of frames):
$(grep -E "[0-9] %$" "$LOGS/contacts.log" | grep -v -- "->" | head -n 8 | sed 's/^ */    /')"
    fi
  else
    plot "$RES/rmsd_backbone.xvg" --unit A --avg 10 -o "$RES/5_rmsd.png"
  fi
  echo "$sum" | summary analyse
  ok "plots and tables saved in results/"
  cd "$HERE"; mark_done analyse; clear_after analyse
}

# =====================================================================
# Stage: mmgbsa  (Module 04, Part A)
# =====================================================================
stage_mmgbsa() {
  has_ligand || { info "No ligand: MM/GBSA skipped."; return 0; }
  need analyse
  say "MM/GBSA: binding free energy over $MMGBSA_FRAMES frames"
  local dir="$RES/mmgbsa"; mkdir -p "$dir"
  local nframes; nframes=$(gmx check -f "$RUN/md_noPBC.xtc" 2>&1 | awk '/^Coords/ {print $2}')
  local interval=$(( nframes / MMGBSA_FRAMES )); [ "$interval" -lt 1 ] && interval=1
  # Group numbers of Protein and the ligand in index.ndx (Module 04, Part A3)
  local gp gl
  gp=$(awk '/^\[/ {g=$2; n++} g=="Protein" && /^\[/ {print n-1}' "$RUN/index.ndx" | head -n 1)
  gl=$(awk -v r="$LIGAND_RESNAME" '/^\[/ {g=$2; n++} g==r && /^\[/ {print n-1}' "$RUN/index.ndx" | head -n 1)
  [ -n "$gp" ] && [ -n "$gl" ] || die "Could not find the Protein and $LIGAND_RESNAME groups in run/index.ndx."
  info "receptor = group $gp (Protein), ligand = group $gl ($LIGAND_RESNAME), one frame in every ${interval}"
  cat > "$dir/mmgbsa.in" << EOF
&general
  sys_name   = "${LIGAND_RESNAME}",
  interval   = ${interval},
  PBRadii    = 3,
/
&gb
  igb        = 5,
  saltcon    = ${SALT_CONC},
/
&decomp
  idecomp    = 2,
  dec_verbose = 0,
  print_res  = "within 4",
/
EOF
  cd "$dir" || exit 1
  for f in topol.top posre.itp "${LIGAND_RESNAME}_GMX.itp" "posre_${LIGAND_RESNAME}.itp"; do
    [ -f "$RUN/$f" ] && cp "$RUN/$f" .
  done
  cp "$RUN"/topol_*.itp "$RUN"/posre_*.itp . 2>/dev/null || true   # multi-chain proteins
  run gmx_mmpbsa gmx_MMPBSA -O -i mmgbsa.in -cs "$RUN/md.tpr" -ci "$RUN/index.ndx" \
      -cg "$gp" "$gl" -ct "$RUN/md_noPBC.xtc" -cp topol.top -nogui
  [ -n "$MMPLOT" ] && python "$MMPLOT" > "$LOGS/plot_mmgbsa.log" 2>&1
  local dg; dg=$(awk '/^ΔTOTAL/ {print $2, "+/-", $4; exit}' FINAL_RESULTS_MMPBSA.dat)
  ok "MM/GBSA binding free energy: $dg kcal/mol (see Module 04 for how to interpret it)"
  echo "MM/GBSA: dG bind = $dg kcal/mol over $MMGBSA_FRAMES frames" | summary mmgbsa
  rm -f _GMXMMPBSA_*
  cd "$HERE"; mark_done mmgbsa
}

# =====================================================================
# Main
# =====================================================================
[ "${1:-}" != "check" ] && mkdir -p "$RUN" "$RES" "$LOGS"   # "check" creates nothing
case "${1:-}" in
  check)        stage_check ;;
  prepare)      stage_prepare ;;
  minimize)     stage_minimize ;;
  equilibrate)  stage_equilibrate ;;
  production)   stage_production ;;
  analyse|analyze) stage_analyse ;;
  mmgbsa)       stage_mmgbsa ;;
  all)
    for s in "${STAGES[@]}"; do
      if [ "$s" = mmgbsa ] && [ "$RUN_MMGBSA" != "yes" ]; then continue; fi
      if is_done "$s"; then info "stage '$s' already finished: skipped"; continue; fi
      "stage_$s"
    done
    say "ALL DONE. Your results are in the results/ folder:"
    ls "$RES" | sed 's/^/      /'
    echo; cat "$RES/summary.txt" | sed 's/^/    /' ;;
  *)
    sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'
    exit 1 ;;
esac
exit 0
