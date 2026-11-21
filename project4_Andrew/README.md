# Project 4 - PDE-Constrained Drone Arm Optimization

This repository is a complete starter implementation for Project 4. It contains a custom sparse 3D frame finite-element solver, static/vibration/buckling analyses, adjoint/eigenvalue sensitivities, SLSQP optimization, verification tests, plots, a generated Markdown report, and a GitHub Pages Three.js viewer.

## Run everything

```bash
python -m pip install -r requirements.txt
python run_all.py
```

Generated outputs include:

- `REPORT.md` - rubric-aligned engineering report
- `index.md` - GitHub Pages version of the report
- `assets/results.json` - single source of truth for report/viewer results
- `assets/figures/*.png` - convergence, radius distribution, Taylor test
- `results/verification.json` - V1-V9 verification ladder
- `results/gradient_checks.json` - adjoint vs central-difference checks
- `results/benchmark.json` - adjoint vs finite-difference timing
- `results/optimization_history.csv` - iteration history
- `results/uniform_feasible.json` - lightest feasible constant-radius reference design

## Main source files

- `src/frame3d.py` - 3D Euler-Bernoulli frame element matrices
- `src/drone_arm.py` - drone-arm assembly, statics, stress, vibration, buckling
- `src/sensitivities.py` - reduced/adjoint gradients and eigenvalue sensitivities
- `src/verify.py` - verification ladder
- `src/optimize.py` - scaled SLSQP optimization
- `src/export_results.py` - viewer/result export
- `src/report_builder.py` - reproducible report generation

## Baseline engineering case

The model represents one arm of a heavy-lift / professional-filming multirotor with six coaxial propulsion stations. Each arm is 2 ft (0.6096 m) long, uses a 6061-T6 aluminum hollow circular tube with a fixed 1.5 mm wall, and is fully clamped at the main frame. The coaxial pair produces 30 lbf (133.4 N) combined maximum thrust, has a combined propulsion-system tip mass of 1.5 lb (0.680 kg), and operates at 3500 RPM. Counter-rotation is assumed to cancel net motor reaction torque. The vibration requirement is 65 Hz, about 11.4% above the 58.3 Hz shaft frequency. A secondary team-selected 2g horizontal-maneuver case is aligned with the arm axis, producing a conservative compressive stability load for the N1 study.

The initial uniform tube radius is 35 mm and the optimization allows 8-40 mm outer-radius variables. These wider bounds are necessary for the updated industrial case: a uniform 18 mm-radius tube has a first natural frequency of only about 27.4 Hz with the 1.5 lb tip mass, well below the 65 Hz requirement.

The current implementation assembles global stiffness and mass matrices sparsely, uses sparse direct static/adjoint solves and a sparse vibration eigensolver, and supplies all verified constraint gradients to SLSQP at every optimization point.
