# The Insilico Lab: Molecular Dynamics Simulation Course

Hands-on training in molecular dynamics simulations of proteins and
protein-ligand complexes, from a bare PDB file to analyzed trajectories.

You do not need Linux on your computer. Everything runs in your web
browser on a free cloud machine with all the software already installed.

## Getting started (first session)

1. Create a free GitHub account at https://github.com/signup if you
   do not already have one, and verify your email address.
2. On this page, click the green **Code** button, open the
   **Codespaces** tab, and click **Create codespace on main**.
3. Wait for the environment to load. The first launch takes a few
   minutes; later launches are faster.
4. When you see the terminal at the bottom of the screen, type:

   ```
   md-check
   ```

   and press Enter. After about a minute you should see
   **All 18 checks passed**. If any check fails, copy the full output
   and post it in the course help channel.

## Rules for your free hours

Your GitHub account includes 120 core-hours of Codespaces each month,
which is 60 hours on the 2-core machine this course uses. The course
needs well under that, as long as you follow three rules:

- **Stop your codespace when you finish a session.** Go to
  https://github.com/codespaces, click the three dots next to your
  codespace, and choose **Stop codespace**. Idle codespaces also stop
  automatically after 30 minutes.
- **Keep one codespace for the whole course.** Reopen the same one each
  session from https://github.com/codespaces instead of creating a new one.
- **Download your results, then delete the codespace at the end of the
  course.** Stopped codespaces still use your free storage.

Do not change the machine type to 4 or 8 cores. Larger machines use up
your free hours two to four times faster.

## Course structure

| Module | Topic |
|---|---|
| 00 | Working in the Linux terminal |
| 01 | Protein in water: preparation, minimization, equilibration, production |
| 02 | Protein-ligand complex: ligand parameterization with ACPYPE/GAFF2 |
| 03 | Trajectory analysis: RMSD, RMSF, radius of gyration, hydrogen bonds |
| 04 | Going further: binding free energy and running on HPC clusters |

Materials for each module are added to the `modules/` folder before the
session.

## Software in this environment

GROMACS 2024.5, AmberTools 24.8, ACPYPE, Open Babel, PDBFixer, OpenMM,
MDAnalysis, NGLView, py3Dmol, and JupyterLab.

---
