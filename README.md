# ROBOT_SOUPLE — Project README ✅

## Overview
**Purpose:** This repository contains MATLAB/Octave code for numerical simulation and control experiments on flexible / beam-like systems and a set of example scripts for a soft-robot ("robot souple") case study. The codebase focuses on assembling finite-element style system matrices, time-integration with Newmark methods, Heaviside-based excitation functions, temporal noise generation, and controller/estimation experiments.

> ⚠️ Note: The descriptions below are based on the filenames and typical numerical workflows present in the folder. For precise behavior, inspect each script's header comments and input/output arguments.

---

## Aims & Goals 🎯
- Implement and test time-integration schemes for structural dynamics (Newmark family).
- Provide utilities for building system matrices (mass, stiffness, loads) for beam-like models.
- Offer smoothed Heaviside functions and their derivatives/integrals for controlled excitations.
- Generate reproducible temporal noise for Monte Carlo/robustness studies.
- Explore control and estimation strategies for a soft-robot example (Kalman, PID, optimal control).

---

## File-by-file summary 🔍
- `Bconstruct.m` — Assembles a B matrix used in the model discretization (boundary or operator matrix used in system equations).
- `Fconstruct.m` — Builds forcing / load vectors (time- or space-dependent loads) for simulations.
- `Heaviside.m` — Implementation of a smoothed Heaviside (step) function used to apply time/space-limited inputs.
- `intHeaviside.m` — Integral (or smoothed integral) of the Heaviside function useful for ramped inputs or modal forcing.
- `difHeaviside.m` — Derivative (or smoothed derivative) of the Heaviside function (e.g., an impulse-like excitation).
- `generateNoiseTemporal.m` — Generates temporal noise sequences and (optionally) uses saved random state files for reproducibility.
- `neb_beam_matrices.m` — Assembles beam finite-element-like matrices (mass, stiffness, possibly damping) used by the time integrators.
- `newmark1stepMRHS.m` — A Newmark integrator variant (one-step form) likely implementing time-stepping with a modified right-hand side.
- `Newmark2N.m` — A different Newmark implementation / variant (2N scheme) for second-order dynamics.

Example / experiment scripts for the soft-robot (tp_robot_souple):
- `tp_robot_souple_base.m` — Base simulation example that sets up the model, initial conditions and runs a base experiment.
- `tp_robot_souple_kalman_pid.m` — Combines Kalman filtering with PID control; used to test estimation + control.
- `tp_robot_souple_kalman_ua.m` — A Kalman-based estimator variant ("ua" may indicate an uncertainty-aware variant); check script for specifics.
- `tp_robot_souple_noncausal.m` — Experiments with noncausal reconstruction or reference trajectories.
- `tp_robot_souple_optimal.m` — Optimal control experiments (e.g., LQR or other optimal control strategies).

Other:
- `octave-workspace` — Saved workspace for Octave/MATLAB sessions.
- `old/` — Older versions of Newmark and other scripts (`Newmark.m`, `Newmark2.m`, `forthodir.m`) kept for reference.
- `randomstate_head/` and `randomstate_nohead/` — Saved `.mat` files containing RNG states for reproducible noise (`randomState.mat`, etc.).

---

## How to run / quick start 🔧
1. Open the repository in MATLAB or GNU Octave.
2. Add the project folder to the MATLAB/Octave path (or `cd` into the folder).
3. Start by running `tp_robot_souple_base.m` to see the typical setup and workflow.
4. Use the `randomState*.mat` files and `generateNoiseTemporal.m` if you need reproducible random inputs.
5. Inspect `neb_beam_matrices.m`, `Bconstruct.m`, and `Fconstruct.m` to understand how the model matrices and loads are assembled.

> Tip: Read top-of-file comments in each `.m` for inputs, outputs, and examples — they often include usage instructions.

---

## Dependencies
- MATLAB (recommended) or GNU Octave should be sufficient for most scripts.
- The code uses standard numerical routines; specific toolboxes are likely *not* required, but check scripts for calls to specialized toolbox functions (Control System Toolbox, Optimization Toolbox, etc.).

---

## Reproducibility & tests ✅
- Use the provided `randomState*.mat` files to restore RNG state when reproducing experiments.
- Compare results across `Newmark` variants (`newmark1stepMRHS.m`, `Newmark2N.m`, and files in `old/`) to validate integrator behavior.

---

## Next steps & suggestions 💡
- Add a short `run_all.m` or `examples/` folder with concise example scripts that produce plots and simple result comparisons.
- Add documentation headers to functions that don't have clear comments (inputs/outputs, units, assumptions).
- Optionally add unit tests (a simple script that verifies energy behavior, or convergence with decreasing dt).

---

## License & Contribution
- Add a `LICENSE` file if you plan to publish or share the repository.
- If you'd like, I can add a short `CONTRIBUTING.md` and some example runs / plots.

---

If you'd like, I can: (1) refine the file descriptions by reading top-of-file comments and adding precise inputs/outputs, or (2) add a `run_examples.m` script that demonstrates the main workflows. Which would you prefer? 
