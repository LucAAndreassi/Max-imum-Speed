# Project 2 DBF Ill-Conditioning — Group Handoff

## Project Goal

This project reuses the Project 1 Design/Build/Fly aircraft optimization and studies **ill-conditioning caused by a quadratic penalty reformulation** of the Mission 2 gross-weight constraint.

The numerical mechanism is **Family G: penalty/barrier reformulation**.

The main idea is:

1. Start with the DBF aircraft optimization problem.
2. Enforce the Mission 2 gross-weight constraint with a quadratic exterior penalty.
3. Increase the penalty weight `rho`.
4. Show that the Hessian condition number grows with `rho`.
5. Show that the large condition number survives Jacobi rescaling.
6. Demonstrate that projected gradient descent slows down.
7. Replace the large fixed penalty with an augmented Lagrangian method.
8. Compare the before/after conditioning and convergence behavior.

---

## Current Main Files

The repository should contain:

```text
Project2DBFIllConditioning(1).md
Project2Models.py
Project2Diagnostics.py
project2_outputs/
```

The report currently being edited is:

```text
Project2DBFIllConditioning(1).md
```

The two code files are:

- `Project2Models.py` — smooth DBF model and optimization objectives.
- `Project2Diagnostics.py` — D1-D4 diagnostics, plots, CSV generation, gradient descent tests, and augmented Lagrangian study.

The output folder should contain:

```text
project2_outputs/D1_hessian_spectrum.png
project2_outputs/D2_kappa_vs_rho.png
project2_outputs/D3_gradient_descent_convergence.png
project2_outputs/D4_augmented_lagrangian_constraint.png
project2_outputs/D4_before_after_convergence.png
project2_outputs/conditioning_vs_rho.csv
project2_outputs/gradient_descent_summary.csv
project2_outputs/augmented_lagrangian_history.csv
project2_outputs/d4_before_after_summary.csv
```

---

## Current Formulation

The six continuous design variables are:

```math
\mathbf{x} =
\left[
S_w,\;
L_F,\;
SW_2,\;
SW_3,\;
V_2,\;
V_3
\right]^T
\in \mathbb{R}^6
```

where:

- `S_w` = wing area
- `L_F` = fuselage length
- `SW_2` = Ground Mission / Mission 2 sensor weight
- `SW_3` = Mission 3 sensor-weight reduction
- `V_2` = Mission 2 cruise velocity
- `V_3` = Mission 3 cruise velocity

The Project 2 payload upper bound was intentionally increased from 12 lb to 16 lb so the 20 lb Mission 2 gross-weight constraint becomes active before the payload box bound.

The physical inequality constraint is:

```math
g(\mathbf{x}) =
W_{\mathrm{M2}}(\mathbf{x}) - 20
\le 0
```

with:

```math
W_{\mathrm{M2}} =
W_e + 1.1\,SW_2
```

The quadratic exterior penalty formulation is:

```math
F_{\rho}(\mathbf{x}) =
f(\mathbf{x})
+
\frac{\rho}{2}
\left[
\max\left(0,\,g(\mathbf{x})\right)
\right]^2
```

---

## Why the Problem Becomes Ill-Conditioned

When the weight constraint is active, the penalty term is:

```math
P(\mathbf{x}) =
\frac{\rho}{2}\,g(\mathbf{x})^2
```

Its local Hessian contribution is:

```math
\nabla^2 P =
\rho\,\nabla g\,\nabla g^T
+
\rho\,g\,\nabla^2 g
```

Near the active boundary, `g(x) ≈ 0`, so the dominant term is approximately:

```math
\nabla^2 P
\approx
\rho\,\nabla g\,\nabla g^T
```

This introduces a very stiff curvature direction normal to the constraint while directions tangent to the constraint remain much softer.

As a result:

```math
\kappa(H) \propto \rho
```

approximately over the tested range.

---

## Key Numerical Results

### D1 — Hessian Spectrum

At `rho = 10,000`:

- smallest positive eigenvalue ≈ `0.1166`
- largest eigenvalue ≈ `1.98e6`
- Hessian condition number ≈ `1.70e7`

### D2 — Intrinsic Ill-Conditioning Test

| rho | Weight residual g (lb) | kappa(H) | kappa after Jacobi |
|---:|---:|---:|---:|
| 0.1 | 1.03646 | 2.09e2 | 1.22e1 |
| 1 | 0.10278 | 1.71e3 | 1.95e1 |
| 10 | 0.01027 | 1.70e4 | 1.60e2 |
| 100 | 0.00103 | 1.70e5 | 1.57e3 |
| 1,000 | 0.000103 | 1.70e6 | 1.57e4 |
| 10,000 | 0.0000102 | 1.70e7 | 1.57e5 |

This satisfies both required parts of the intrinsic-conditioning test:

1. `kappa` grows as the structural knob `rho` increases.
2. The large condition number survives Jacobi rescaling.

---

## D3 — Effect on Projected Gradient Descent

Convergence tolerance:

```math
\left\|\nabla_P F_{\rho}(\mathbf{x}_k)\right\|_2
\le 10^{-5}
```

Maximum iterations: `15,000`

| rho | Iterations | Final projected-gradient norm | Final objective gap |
|---:|---:|---:|---:|
| 0.1 | 849 | 9.97e-6 | 1.14e-11 |
| 1 | 9,668 | 1.00e-5 | 1.00e-10 |
| 10 | 15,000* | 3.49e-1 | 3.43e-4 |

`rho = 10` hit the iteration limit before converging.

This is the main evidence that increasing penalty-induced conditioning directly slows the baseline first-order optimizer.

---

## D4 — Augmented Lagrangian Remedy

The augmented Lagrangian formulation is:

```math
\mathcal{L}_A(\mathbf{x}, \lambda, \rho)
=
f(\mathbf{x})
+
\frac{\rho}{2}
\left[
\max\left(
0,\,
g(\mathbf{x}) + \frac{\lambda}{\rho}
\right)
\right]^2
-
\frac{\lambda^2}{2\rho}
```

Multiplier update:

```math
\lambda_{k+1}
=
\max\left(
0,\,
\lambda_k + \rho\,g(\mathbf{x}_k)
\right)
```

Comparison:

### Fixed penalty, rho = 10,000

- constraint residual ≈ `1.02e-5 lb`
- condition number ≈ `1.70e7`

### Augmented Lagrangian

- started with `rho = 5`
- feasible after 2 outer iterations
- final residual ≈ `-5.11e-6 lb`
- local condition number ≈ `8.51e3`

The augmented Lagrangian reaches essentially the same constraint accuracy while reducing the local condition number by about three orders of magnitude.

---

## D4 Direct Convergence Comparison

After 20,000 projected-gradient iterations from equal-size normalized perturbations:

| Formulation | Final projected-gradient norm | Final objective gap | Relative objective gap |
|---|---:|---:|---:|
| Fixed penalty, rho = 10,000 | 1.33e-1 | 2.10e-3 | 9.92e-6 |
| Augmented-Lagrangian local objective | 8.33e-4 | 2.05e-9 | 1.90e-8 |

The augmented-Lagrangian local formulation has:

- about 160x smaller projected-gradient norm
- more than 500x smaller relative objective gap

at the same iteration count.

---

## Final Augmented-Lagrangian Design

Approximate final design:

| Variable | Value |
|---|---:|
| Wing area | 5.96 ft^2 |
| Fuselage length | 5.06 ft |
| Mission 2 sensor weight | 12.50 lb |
| Mission 3 sensor reduction | 0.87 lb |
| Mission 2 cruise velocity | 140.0 ft/s |
| Mission 3 cruise velocity | 85.48 ft/s |
| Empty aircraft weight | 6.25 lb |
| Mission 2 gross weight | 20.00 lb |
| Mission 2 five-lap time | 123.46 s |
| Continuous Mission 3 lap estimate | 7.21 laps |

Final unpenalized objective:

```math
f(\mathbf{x}) = -3.39336
```

---

## Important Modeling Changes from Project 1

The Project 2 model is intentionally smoothed because finite-difference gradients and Hessians are required.

Changes include:

1. Mission 3 lap count is continuous instead of using `floor()`.
2. Throttle is solved continuously instead of selected from a 100-point throttle grid.
3. Turn load factor is handled continuously instead of repeated discrete decrements.
4. Smooth limiting expressions replace abrupt switching where necessary.
5. Conditioning is evaluated in normalized design coordinates.
6. The Project 2 `SW_2` upper bound is 16 lb instead of 12 lb so the 20 lb gross-weight constraint becomes active.

These should be described transparently in the report.

---

## GitHub Math Formatting

GitHub has been inconsistent with some `$$ ... $$` blocks.

Preferred display format for equations is:

````text
```math
equation here
```
````

Avoid macros GitHub rejects.

Known issue encountered:

```text
\operatorname{diag}(H)
```

was rejected.

Use:

```math
D = \mathrm{diag}(H)
```

instead.

For units, prefer:

```math
1.02 \times 10^{-5}\,\mathrm{lb}
```

or plain text outside the equation if necessary.

Before submission, open the Markdown file directly on GitHub and visually inspect every equation.

---

## Current Rubric Status

The current report is technically strong against the provided rubric.

| Category | Max | Status |
|---|---:|---|
| Motivation / real-world relevance | 10 | Strong |
| Formulation | 20 | Strong |
| Ill-conditioning mechanism + intrinsic test | 20 | Strong |
| D1 / D3 effect demonstration | 20 | Strong |
| D4 solution and before/after comparison | 20 | Strong |
| Reproducibility | 5 | Strong |
| Presentation / clarity | 5 | Main remaining cleanup area |

The biggest remaining risk is GitHub rendering, not missing technical content.

---

## Remaining Work

1. **GitHub rendering pass**
   - Open the report on GitHub.
   - Check every equation.
   - Convert broken `$$ ... $$` blocks to fenced `math` blocks.
   - Check all figures load correctly.

2. **Technical consistency pass**
   - Make sure `SW_3` is consistently described as Mission 3 sensor-weight reduction.
   - Make sure all references to the 20 lb constraint use Mission 2 gross weight.
   - Keep the distinction between the original constrained problem and the penalized problem clear.

3. **Hessian wording**
   - The report refers to the smallest positive Hessian eigenvalue.
   - If needed, add a brief explanation that conditioning is evaluated on the positive-curvature subspace near the studied minimizer.

4. **Figure check**
   - Confirm every image path is correct.
   - Confirm captions are readable.
   - Make sure axes and legends are visible at GitHub display size.

5. **Reproducibility check**
   - Clone or download the repo into a fresh folder.
   - Run:

```bash
python -m pip install numpy scipy matplotlib
python Project2Diagnostics.py
```

   - Confirm the CSVs and figures regenerate successfully.

6. **Final proofread**
   - Check grammar.
   - Remove any stale references to older filenames such as `Project2DBFIllConditioning_updated.md` if the final file has a different name.
   - Make sure the report file name used in the repository matches any references in the text.

---

## Suggested Division of Group Work

- **Person 1:** GitHub Markdown/math rendering and figure paths
- **Person 2:** Technical proofread of formulation, constraints, and notation
- **Person 3:** Run code from a clean environment and verify outputs
- **Person 4:** Review D1-D4 discussion against the rubric
- **Person 5:** Final grammar, consistency, and presentation pass

If the group is smaller, combine these roles.

---

## Submission Checklist

- [ ] Markdown renders correctly on GitHub
- [ ] All equations display properly
- [ ] All four required diagnostic figures display
- [ ] D1 spectrum shown
- [ ] D2 intrinsic test shown
- [ ] D3 baseline slowdown shown
- [ ] D4 before/after remedy shown
- [ ] All CSV outputs included or reproducible
- [ ] `Project2Models.py` runs/imports correctly
- [ ] `Project2Diagnostics.py` runs successfully
- [ ] No missing package dependencies
- [ ] Report assumptions are included
- [ ] Reproducibility instructions are included
- [ ] Final repository is public
- [ ] Repository link submitted on Canvas

---

## Important Caution

Do not change the numerical results in the report without rerunning `Project2Diagnostics.py`.

If any model equation, bound, smoothing assumption, tolerance, or optimizer setting is changed, regenerate the figures and CSV files and update the report numbers together so they remain consistent.
