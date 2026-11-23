# Project 4 submission and evidence map

Requirements source: [Project 4 overview, report requirements and rubric](https://designinformaticslab.github.io/DesignOptimization2025/project4.html), reviewed 23 November 2026. No separate local rubric was provided. This is a file/evidence map, not a claim that all grading requirements are satisfied.

| Rubric category | Retained evidence | Status / action |
|---|---|---|
| Formulation (10) | Report sections 1-3; both input JSON cards | Variables, units, bounds, loads and estimated material sources recorded; confirm scope and justified limits |
| Analysis and verification (20) | `arm.py`, `vibration.py`, `verify.py`, report verification table; `drone_arm.py` comparison | Core ladder passes; reduced solver and separate 3D comparison must be distinguished |
| Sensitivities (20) | `optimize.py`, `arm.py`, `vibration.py`, tests, report derivations, gradient tables/plots and timing | Evidence retained; human derivation/review still needed |
| Optimization (15) | `optimize_arm.py`, `output/optimization/`, report KKT and convergence sections | Current 1.5 lb case; multistart, feasibility, active constraints, multipliers and shadow checks |
| Nonlinear component (15) | `arm.Elastica`, nonlinear adjoints, independent BVP checks, with/without comparison | N2 evidence retained; reported design change may be negligible |
| Website / 3D comparison (10) | `index.md`, `_layouts/`, `_config.yml`, `assets/`, `preview.html` | Generated report and viewer present; test live public Pages site, mobile layout and all required viewer metrics |
| Reproducibility (5) | `requirements.txt`, README, input cards, source and result fingerprints | `python optimize_arm.py --update-report` rebuilds site numbers in an installed environment |
| Presentation (5) | Report, README, this checklist | Review as a team before submission |

## Outstanding scope and submission decisions

- The physical design has **two uniform-tube sizing variables**. The course lists a 20-variable threshold for proposed structures. Confirm expectations for this reduced drone-arm scope with the instructor; do not create artificial variables just to increase the count.
- The main optimization uses a scalar bending FEM plus a planar nonlinear model. The supplied 3D frame is a separate comparison. Do not claim a fully verified general 3D-frame optimization or six-rigid-mode/rotation-invariance coverage from the scalar ladder.
- Vibration constraints currently use the exact NACA uniform-beam tip-mass equation, with the FE eigenproblem retained as a check. Confirm this reduced-model choice is acceptable for the course's own-FEM/eigenproblem expectations.
- No stiffness/deflection limit or operating frequency-separation margin was supplied. The enabled frequency threshold is f1 >= RPM/60. Confirm the intended engineering requirements.
- Applicable laminate, local shell buckling, joint and fatigue allowables remain unresolved; nominal model SF does not establish overall structural safety.
- Expected verification values and adjoint derivations must receive the team's own review. The report distinguishes agent-written checks from human work.
- Create the public team repository, enable GitHub Pages, verify the live math/viewer in a private browser and at phone width, then submit repository and site links. No publishing or submission was performed by this cleanup.

## Website data discipline

`assets/results.json` is the shared generated numerical source for the report pipeline and live viewer. `output/optimization/results.json` retains the upstream optimization audit trail. They serve different roles and are both needed for an immediately usable, reproducible handoff. `assets/report.md` is required by the local preview; `index.md` is required by GitHub Pages. Do not remove these as apparent duplicates.

## Workspace cleanup

Keep current source, input cards, verification/tests, report template, website dependencies and the three current output directories listed in the README. Historical stock-tube exports, old tip-mass studies, the earlier strength-only output, temporary paper downloads and environment/cache files are not current submission evidence. A recovery archive and manifest are saved outside this workspace; `.git/` remains intact.
