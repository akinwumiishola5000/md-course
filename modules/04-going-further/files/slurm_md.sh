#!/bin/bash
# slurm_md.sh - template for running a GROMACS production simulation on an
# HPC cluster that uses the SLURM job scheduler.
#
# Submit it with:      sbatch slurm_md.sh
# Check your jobs:     squeue -u $USER
# Cancel a job:        scancel JOBID
#
# Every cluster is different. Before using this, read your cluster's
# documentation for: the partition (queue) names, the account name,
# the GPU request syntax, and how GROMACS is loaded (module names).

#SBATCH --job-name=md_complex        # name shown in squeue
#SBATCH --partition=gpu              # CHANGE: your cluster's GPU partition
#SBATCH --account=my_account         # CHANGE: your allocation/account (delete if not used)
#SBATCH --nodes=1                    # one computer
#SBATCH --ntasks=1                   # one GROMACS process
#SBATCH --cpus-per-task=8            # CPU cores for that process
#SBATCH --gres=gpu:1                 # one GPU (syntax varies between clusters)
#SBATCH --mem=16G                    # memory
#SBATCH --time=24:00:00              # maximum run time (hh:mm:ss)
#SBATCH --output=md_%j.out           # screen output; %j is replaced by the job number

# ---- Load software (CHANGE: use your cluster's module names) ----
module purge
module load gromacs/2024.5           # example only; check with: module avail gromacs

# ---- Run ----
cd "$SLURM_SUBMIT_DIR"               # the folder you submitted from

# Stop 15 minutes before the time limit so a checkpoint is written cleanly
MAXH=$(awk "BEGIN {print 24 - 0.25}")

# -cpi md.cpt continues from the checkpoint if one exists, so you can
# submit the same script again to extend a run that hit the time limit.
gmx mdrun -deffnm md -cpi md.cpt -maxh "$MAXH" \
          -ntmpi 1 -ntomp "$SLURM_CPUS_PER_TASK" \
          -nb gpu -pme gpu -bonded gpu -update gpu
