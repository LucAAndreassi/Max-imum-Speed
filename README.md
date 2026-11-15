# Project 4 - PDE-Constrained Drone Arm Optimization

This repository is a complete starter implementation for Project 4. It contains a custom 3D frame finite-element solver, static/vibration/buckling analyses, adjoint/eigenvalue sensitivities, SLSQP optimization, verification tests, plots, a generated Markdown report, and a GitHub Pages Three.js viewer.

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
- `results/verification.json` - V1-V6 verification ladder
- `results/gradient_checks.json` - adjoint vs central-difference checks
- `results/benchmark.json` - adjoint vs finite-difference timing
- `results/optimization_history.csv` - iteration history

## Main source files

- `src/frame3d.py` - 3D Euler-Bernoulli frame element matrices
- `src/drone_arm.py` - drone-arm assembly, statics, stress, vibration, buckling
- `src/sensitivities.py` - reduced/adjoint gradients and eigenvalue sensitivities
- `src/verify.py` - verification ladder
- `src/optimize.py` - scaled SLSQP optimization
- `src/export_results.py` - viewer/result export
- `src/report_builder.py` - reproducible report generation

## Baseline engineering case

The model uses a 0.30 m, 20-element 6061-T6 aluminum tube arm with a fixed 1.5 mm wall. The propulsion reference is the T-MOTOR MN3110 KV700: 99 g motor mass and the 14.8 V / 12x4 propeller full-throttle test point (1290 g thrust, 7600 RPM). The 140 Hz minimum structural frequency is a project design requirement chosen to remain above the 126.7 Hz shaft rotational frequency.

Before final submission, the team should confirm that this propulsion system and the selected load case match the intended drone concept.
