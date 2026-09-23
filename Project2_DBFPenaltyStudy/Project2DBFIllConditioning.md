# Project 2 - Ill-Conditioned Optimization of a Design/Build/Fly Aircraft

## Problem Outline
Design/Build/Fly (DBF) is an annual collegiate aircraft design competition hosted by the AIAA (American Institute of Aeronautics and Astronautics). Teams are required to design, build, and fly a remotely piloted aircraft that satisfies a mission-specific ruleset. The competition requires teams to balance structural design, payload capability, aerodynamic performance, propulsion requirements, and mission completion time.

This project reuses the aircraft sizing and mission-performance model developed in Project 1. In Project 1, the objective was to maximize a normalized competition score using aircraft geometry, payload weight, and mission cruise speeds as optimization variables. The present project studies the same design problem from a numerical optimization perspective by introducing a penalty reformulation of the aircraft weight constraint.

The purpose of the reformulation is to demonstrate **intrinsic ill-conditioning**. As the penalty weight is increased, the optimizer is forced more tightly toward the weight constraint, but the local curvature of the optimization landscape becomes increasingly unequal in different directions. This produces a large Hessian condition number and causes a first-order optimization method to converge slowly.

The three scoring missions remain based on the Project 1 model:

**Ground Mission:** The aircraft carries a sensor payload. The score contribution is proportional to sensor weight.

**Mission 2:** The aircraft carries the sensor in its shipping container and completes five laps. The score rewards higher payload weight and shorter completion time.

**Mission 3:** The sensor is deployed and the aircraft flies as many laps as possible within five minutes. The score is proportional to the deployed sensor weight multiplied by the number of laps completed.

---

## Decision Variables
The optimization uses the same six main aircraft and mission variables as Project 1.

| Variable | Definition | Project 2 Bounds | Type |
|---|---|---:|---|
| $S_w$ | Wing area | $2 \le S_w \le 10$ ft² | Continuous |
| $L_F$ | Fuselage length | $5 \le L_F \le 7$ ft | Continuous |
| $SW_2$ | Sensor weight for Ground Mission and Mission 2 | $4 \le SW_2 \le 16$ lb | Continuous |
| $SW_3$ | Sensor weight reduction for Mission 3 | $0 \le SW_3 \le 6$ lb | Continuous |
| $V_2$ | Mission 2 cruise velocity | $40 \le V_2 \le 140$ ft/s | Continuous |
| $V_3$ | Mission 3 cruise velocity | $40 \le V_3 \le 140$ ft/s | Continuous |

The Project 1 upper bound on $SW_2$ was 12 lb. For this conditioning study, the search bound is intentionally increased to 16 lb. This allows the 20 lb Mission 2 gross-weight constraint to become active and control the optimum rather than having the payload box bound become active first. This is a Project 2 study modification and is not intended to represent a change to the original competition rules.

All other aircraft geometry is still derived from these variables. The wing span remains fixed at 6 ft, the sensor is modeled as a 6:1 cylindrical payload, and the empennage and fuselage dimensions are determined from the selected design variables.

---

## Objective Function
The base objective is the same normalized DBF competition score used in Project 1. Because SciPy optimizers minimize functions, the competition score is negated:

$$
f(\mathbf{x})
= -\left(
\frac{SW_2}{12}
+ \frac{SW_2/t_2}{0.15}
+ \frac{(SW_2-SW_3)L_3}{50}
\right)
$$

where:

- $t_2$ is the time required to complete five Mission 2 laps,
- $L_3$ is the predicted number of Mission 3 laps,
- $\mathbf{x}=[S_w,L_F,SW_2,SW_3,V_2,V_3]^T$.

For Project 2, a quadratic exterior penalty is added to enforce the Mission 2 aircraft gross-weight limit. Define the constraint as

$$
g(\mathbf{x}) = W_{M2}(\mathbf{x}) - 20 \le 0,
$$

where

$$
W_{M2}=W_e+1.1SW_2.
$$

The factor of 1.1 accounts for the modeled shipping-container weight used in Project 1.

The penalized objective is therefore

$$
F_\rho(\mathbf{x})
=
f(\mathbf{x})
+
\frac{\rho}{2}
\left[\max(0,g(\mathbf{x}))\right]^2.
$$

The penalty parameter $\rho$ controls how strongly violations of the 20 lb weight target are punished. Increasing $\rho$ forces the optimized design closer to the constraint boundary, but also creates the ill-conditioning studied in this project.

---

## Constraints
The main constraint used to create the ill-conditioned penalty formulation is

$$
W_{M2}\le20\text{ lb}.
$$

The optimization also retains box bounds on all six decision variables.

| Constraint | Mathematical Expression | Description |
|---|---|---|
| Mission 2 gross weight | $W_e+1.1SW_2\le20$ lb | Team design target for aircraft and payload weight |
| Wing area | $2\le S_w\le10$ ft² | Limits lifting-surface size |
| Fuselage length | $5\le L_F\le7$ ft | Limits aircraft length |
| Mission 2 payload | $4\le SW_2\le16$ lb | Project 2 study range |
| Mission 3 payload reduction | $0\le SW_3\le6$ lb | Limits payload reduction between missions |
| Mission 2 speed | $40\le V_2\le140$ ft/s | Cruise-speed range |
| Mission 3 speed | $40\le V_3\le140$ ft/s | Cruise-speed range |

---

## Aircraft Model and Project 2 Simplifications
The Project 2 model preserves the main physical relationships from Project 1, including:

- sensor sizing from payload weight and lead-shot density,
- wing geometry based on fixed 6 ft span,
- empirical fuselage, wing, electronics, and empennage weight estimates,
- induced and parasite drag,
- propulsion current as a function of required thrust and velocity,
- straight and turning flight times,
- battery-energy limitation,
- normalized DBF mission scoring.

Several numerical changes were made because Project 2 requires gradients, Hessians, and eigenvalue analysis. The original model contained discrete operations that are appropriate for predicting competition scoring but unsuitable for smooth conditioning analysis.

The following changes were made:

1. Mission 3 lap count is treated continuously instead of using `floor()`.
2. Throttle is solved continuously rather than selected from a 100-point throttle array.
3. The turn-load-factor reduction is replaced with a continuous $C_{L,max}$-limited calculation.
4. Smooth minimum approximations are used where the original model switched abruptly between limiting cases.

These changes preserve the same engineering relationships while allowing meaningful finite-difference gradients and Hessians to be computed.

---

## Problem Classification
The Project 2 formulation is a **constrained, nonlinear, nonconvex optimization problem** that is studied through a smooth penalty reformulation.

The base aircraft objective is nonlinear because aerodynamic drag, geometry, propulsion demand, mission time, energy consumption, and score all depend nonlinearly on the design variables. The optimization is nonconvex because the coupled aircraft-performance model can create multiple local optima.

For the conditioning study, the original inequality constraint is moved into the objective through a quadratic exterior penalty. This places the problem in **Family G: penalty/barrier reformulations** from the Project 2 classification.

---

## Ill-Conditioning Mechanism
The source of ill-conditioning is the quadratic penalty term

$$
P(\mathbf{x})
=
\frac{\rho}{2}g(\mathbf{x})^2
$$

when the weight constraint is active. The Hessian contribution from this penalty is approximately

$$
\nabla^2P
\approx
\rho\nabla g\nabla g^T
+
rho g\nabla^2g.
$$

Near the constraint boundary, $g\approx0$, so the dominant term is

$$
\nabla^2P
\approx
\rho\nabla g\nabla g^T.
$$

This produces a very stiff direction normal to the weight-constraint surface. Directions tangent to the constraint remain governed mostly by the original DBF objective and therefore retain much smaller curvature.

As $\rho$ increases, the largest Hessian eigenvalue grows approximately in proportion to $\rho$, while the smallest positive eigenvalue remains nearly unchanged. The Hessian condition number

$$
\kappa(H)=\frac{\lambda_{max}(H)}{\lambda_{min}(H)}
$$

therefore grows approximately linearly with the penalty weight.

This is not simply a mismatch in engineering units. The Project 2 intrinsic-conditioning test requires the large condition number to survive diagonal, or Jacobi, rescaling. The numerical results below show that it does.

---

## D1 - Hessian Eigenvalue Spectrum
The Hessian was evaluated numerically near the optimum of the penalized problem. For the largest tested penalty, $\rho=10,000$, the smallest positive Hessian eigenvalue is approximately

$$
\lambda_{min}=0.1166,
$$

while the largest is approximately

$$
\lambda_{max}=1.98\times10^6.
$$

This gives

$$
\kappa(H)\approx1.70\times10^7.
$$

The spectrum spans many orders of magnitude, which indicates a strongly elongated local optimization landscape.

![D1 Hessian Spectrum](project2_outputs/D1_hessian_spectrum.png)

---

## D2 - Intrinsic Ill-Conditioning Test
The penalty weight was varied over five orders of magnitude. At every value of $\rho$, the penalized problem was optimized and the local Hessian was evaluated at the resulting solution.

| $\rho$ | Weight Constraint Residual $g$ (lb) | $\kappa(H)$ | $\kappa$ After Jacobi Rescaling |
|---:|---:|---:|---:|
| 0.1 | 1.03646 | $2.09\times10^2$ | $1.22\times10^1$ |
| 1 | 0.10278 | $1.71\times10^3$ | $1.95\times10^1$ |
| 10 | 0.01027 | $1.70\times10^4$ | $1.60\times10^2$ |
| 100 | 0.00103 | $1.70\times10^5$ | $1.57\times10^3$ |
| 1,000 | 0.000103 | $1.70\times10^6$ | $1.57\times10^4$ |
| 10,000 | 0.0000102 | $1.70\times10^7$ | $1.57\times10^5$ |

Two effects are visible. First, increasing $\rho$ makes the weight constraint increasingly accurate. Second, the condition number grows by approximately one order of magnitude each time $\rho$ increases by one order of magnitude.

The large condition number also survives Jacobi rescaling. At $\rho=10,000$, diagonal scaling reduces the numerical value of $\kappa$, but it remains approximately

$$
\kappa_{Jacobi}\approx1.57\times10^5.
$$

Therefore, the ill-conditioning is intrinsic to the penalty formulation rather than being caused only by different variable units or coordinate scales.

![D2 Condition Number vs Penalty Weight](project2_outputs/D2_kappa_vs_rho.png)

---

## D3 - Effect on Gradient Descent
To demonstrate the practical effect of the growing condition number, projected gradient descent was applied to the penalized objective. The variables were normalized to the interval $[0,1]$ before the gradient-descent experiment so that the observed slowdown was not dominated by engineering-unit differences.

The gradient-descent step size was selected from the local Hessian eigenvalues using

$$
\alpha=\frac{2}{L+\mu},
$$

where

$$
L=\lambda_{max}(H),
\qquad
\mu=\lambda_{min}(H).
$$

The results were:

| $\rho$ | Iterations | Final Projected Gradient Norm | Final Objective Gap |
|---:|---:|---:|---:|
| 0.1 | 849 | $9.97\times10^{-6}$ | $1.14\times10^{-11}$ |
| 1 | 9,668 | $1.00\times10^{-5}$ | $1.00\times10^{-10}$ |
| 10 | 15,000* | $3.49\times10^{-1}$ | $3.43\times10^{-4}$ |

\*The $\rho=10$ case reached the 15,000-iteration limit before satisfying the convergence tolerance.

The increase from 849 iterations at $\rho=0.1$ to 9,668 iterations at $\rho=1$ demonstrates the slowdown caused by the increasingly elongated optimization landscape. At $\rho=10$, the baseline method no longer converged within the allowed iteration count.

![D3 Gradient Descent Convergence](project2_outputs/D3_gradient_descent_convergence.png)

---

## Proposed Solution - Augmented Lagrangian Method
A fixed quadratic penalty creates a tradeoff. A small value of $\rho$ gives a manageable optimization landscape but allows noticeable constraint violation. A very large $\rho$ enforces the constraint accurately but makes the problem severely ill-conditioned.

The proposed remedy is the **augmented Lagrangian method**. For the inequality constraint $g(\mathbf{x})\le0$, the implemented objective is

$$
\mathcal{L}_A(\mathbf{x},\lambda,\rho)
=
f(\mathbf{x})
+
\frac{\rho}{2}
\left[\max\left(0,g(\mathbf{x})+\frac{\lambda}{\rho}\right)\right]^2
-
\frac{\lambda^2}{2\rho}.
$$

After each inner optimization, the multiplier is updated according to

$$
\lambda_{k+1}
=
\max\left(0,\lambda_k+\rho g(\mathbf{x}_k)\right).
$$

The multiplier allows the algorithm to enforce the constraint without requiring the penalty parameter to become extremely large.

---

## D4 - Before and After Comparison
The fixed-penalty solution with $\rho=10,000$ produced

$$
g(\mathbf{x})=1.02\times10^{-5}\text{ lb}
$$

and

$$
\kappa(H)\approx1.70\times10^7.
$$

The augmented-Lagrangian method started with $\rho=5$ and reached a feasible solution after two outer iterations. Its final residual was

$$
g(\mathbf{x})=-5.11\times10^{-6}\text{ lb},
$$

with local condition number approximately

$$
\kappa(H_{AL})\approx8.51\times10^3.
$$

Thus, the augmented-Lagrangian method achieved essentially the same constraint accuracy while reducing the local condition number by roughly three orders of magnitude.

![D4 Augmented Lagrangian Constraint Convergence](project2_outputs/D4_augmented_lagrangian_constraint.png)

The final augmented-Lagrangian design was approximately:

| Variable | Final Value |
|---|---:|
| Wing area, $S_w$ | 5.96 ft² |
| Fuselage length, $L_F$ | 5.06 ft |
| Mission 2 sensor weight, $SW_2$ | 12.50 lb |
| Mission 3 sensor reduction, $SW_3$ | 0.87 lb |
| Mission 2 cruise velocity, $V_2$ | 140.0 ft/s |
| Mission 3 cruise velocity, $V_3$ | 85.48 ft/s |
| Empty aircraft weight | 6.25 lb |
| Mission 2 gross weight | 20.00 lb |
| Mission 2 five-lap time | 123.46 s |
| Continuous Mission 3 lap estimate | 7.21 laps |

The final unpenalized objective was approximately

$$
f(\mathbf{x})=-3.39336.
$$

---

## Results and Interpretation
The penalty study shows the expected Family G ill-conditioning mechanism. Increasing the quadratic penalty weight reduces the weight-constraint violation, but it also produces a rapidly growing Hessian condition number. The largest curvature direction becomes dominated by the weight penalty, while the remaining directions retain much smaller curvature from the underlying aircraft-performance objective.

The relationship is approximately

$$
\kappa(H)\propto\rho.
$$

The Jacobi-rescaling test confirms that the effect is not merely caused by the decision variables having different units. Although diagonal rescaling improves the numerical conditioning, the scaled problem still becomes increasingly ill-conditioned as $\rho$ grows.

This change in local geometry has a direct effect on a baseline first-order optimizer. Projected gradient descent requires substantially more iterations as $\rho$ increases, and for $\rho=10$ it fails to reach the target tolerance within 15,000 iterations.

The augmented-Lagrangian method resolves the main tradeoff. Instead of using an extremely large fixed penalty, it combines a moderate penalty with a Lagrange-multiplier update. This produces a nearly feasible design while keeping the local condition number far below the value produced by the $\rho=10,000$ fixed-penalty formulation.

For the DBF design problem, this means that directly forcing engineering constraints with increasingly large penalty coefficients can make an otherwise manageable optimization problem unnecessarily difficult to solve. A constraint-handling method such as an augmented Lagrangian can enforce the same physical requirement while maintaining a significantly better-conditioned numerical problem.

---

## Assumptions and Simplifications
The results depend on the aircraft model and several simplifying assumptions:

- The model is a preliminary-design approximation and is not a full CFD or flight-dynamics simulation.
- Wing span is fixed at 6 ft.
- The aircraft uses a single-wing fixed-wing configuration with an inverted T-tail.
- The sensor is represented as a 6:1 lead-shot-filled cylindrical payload with a nosecone.
- Structural weight is estimated using empirical relationships inherited from Project 1.
- Propulsion behavior is based on the Project 1 motor/propeller model.
- Mission 3 lap count is treated continuously for the conditioning analysis rather than rounded to an integer.
- Battery and time limiting behavior are smoothed for differentiability.
- Numerical gradients and Hessians are computed with finite differences, so their accuracy depends on the selected perturbation size.
- The 16 lb Project 2 payload upper bound is used only to expose the 20 lb gross-weight constraint as the active limiting mechanism.
- The augmented-Lagrangian and fixed-penalty comparisons describe the local behavior of this model and do not guarantee the global optimum of the full competition-design problem.

---

## Code and Reproducibility
Two Python scripts are used for the Project 2 study.

[Project2Models.py](Project2Models.py) contains the smooth DBF aircraft model, mission-performance calculations, base competition objective, weight constraint, quadratic penalty function, and augmented-Lagrangian objective.

[Project2Diagnostics.py](Project2Diagnostics.py) performs the numerical experiments required for the conditioning study. It optimizes the penalized objective over a range of $\rho$ values, computes finite-difference Hessians and eigenvalue spectra, performs Jacobi rescaling, runs projected gradient descent, and evaluates the augmented-Lagrangian remedy.

Running `Project2Diagnostics.py` generates the following output files:

- `project2_outputs/D1_hessian_spectrum.png`
- `project2_outputs/D2_kappa_vs_rho.png`
- `project2_outputs/D3_gradient_descent_convergence.png`
- `project2_outputs/D4_augmented_lagrangian_constraint.png`
- `project2_outputs/conditioning_vs_rho.csv`
- `project2_outputs/gradient_descent_summary.csv`
- `project2_outputs/augmented_lagrangian_history.csv`

The numerical outputs used in this report are generated directly from these scripts so that the results can be reproduced from the submitted repository.
