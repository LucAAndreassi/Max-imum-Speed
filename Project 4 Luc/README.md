# Drone arm optimization - submission workspace

Current case: 2 ft uniform CFRP tube, 30 lbf combined thrust, **1.5 lb combined tip mass**, 3,500 RPM, wall > 2 mm, nominal material-strength SF >= 2. The optimized tube is approximately **55.285 mm OD x 2.010 mm wall**, **311.714 g**, with **1.240 mm linear tip deflection** and **58.333 Hz first bending frequency**. Generated results remain the numerical authority.

Material inputs are estimated effective properties. Shell buckling, joints, laminate failure interaction and operating dynamic stability are unverified. See `SUBMISSION_CHECKLIST.md` for course scope and submission gaps.

## Rebuild the report and website

Use Python 3.12. From this folder, create a fresh local environment:

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

One command reruns gradient checks, both optimization cases, mesh refinement, core verification, and regenerates the report, figures and viewer data:

```powershell
.venv\Scripts\python.exe optimize_arm.py --update-report
```

Verification tests:

```powershell
.venv\Scripts\python.exe -m unittest test_optimize test_vibration -v
```

Numerical timings vary by machine. Results include source/input fingerprints; report generation rejects stale optimization results.

## View or publish

```powershell
.venv\Scripts\python.exe -m http.server 8000 --bind 127.0.0.1
```

Open [the local preview](http://127.0.0.1:8000/preview.html). Use an HTTP server rather than opening the file directly. Internet access is needed for the pinned three.js, MathJax and Marked libraries.

For GitHub Pages, publish this repository with Pages building the root of the selected branch. The website uses `index.md`, `_config.yml`, `_layouts/default.html`, and `assets/`. A public repository/site and their submission are still team actions; this workspace has not been deployed.

## Files to keep

| Files | Purpose |
|---|---|
| `inputs.json`, `optimizer_inputs.json` | Loads, geometry bounds, material estimates and sources, performance limits |
| `arm.py`, `vibration.py` | Own beam FEM, finite-rotation follower-thrust model, NACA tip-mass vibration solution |
| `optimize.py`, `optimize_arm.py` | Design-dependent adjoints, multistart optimization, KKT acceptance and exports |
| `verify.py`, `test_optimize.py`, `test_vibration.py` | Verification ladder, independent checks and sensitivity tests |
| `run.py`, `report_template.md` | Report and website generation; edit template rather than generated numbers |
| `index.md`, `preview.html`, `_config.yml`, `_layouts/`, `assets/` | Published report, local preview, live 3D comparison, figures and shared results |
| `analyze.py`, `check_design.py` | Standalone FEM exports and nominal strength screening |
| `drone_arm.py`, `compare_drone_arm.py` | Supplied independent 3D frame solver and reproducible comparison runner |
| `requirements.txt`, `SUBMISSION_CHECKLIST.md` | Environment and rubric/evidence map |

## Current output only

- `output/optimization/`: optimizer results, input card for the optimum, convergence history and plot.
- `output/optimized-fe/`: current standalone FEM mesh, matrices, displacements, modes and results.
- `output/drone-arm-comparison/`: current independent 3D frame comparison and FE data.

Refresh the standalone exports and 3D comparison after changing the design:

```powershell
.venv\Scripts\python.exe analyze.py --inputs output\optimization\optimal_inputs.json --nonlinear --output output\optimized-fe
.venv\Scripts\python.exe compare_drone_arm.py
```

`drone_arm.py` is preserved as supplied and imports an absent `arm_inputs.py`. Run it through `compare_drone_arm.py`, which supplies defaults in memory and passes all inputs explicitly. Its auxiliary Poisson-ratio assumption is recorded and checked for influence on the bending results. It is a linear 3D frame comparison, not the nonlinear optimization engine.

Default vibration uses NACA TN 3459, printed p. 17, equations (13)-(14), for a stationary uniform cantilever with point tip mass. `Beam.fe_modes()` and `analyze.py --vibration-method fem` retain the independent matrix eigenproblem. The motor RPM is a sizing target, not a beam rotation speed. Analytical mode samples and discrete FE eigenvectors are exported separately.

Generated environments, caches, temporary downloads and historical runs are excluded from this submission workspace. Recreate the environment using the commands above; `.git/` is retained for version control.
