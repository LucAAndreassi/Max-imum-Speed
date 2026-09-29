# Project 2 - Ill-Conditioned Optimization of a Design/Build/Fly Aircraft
## 1. Problem Identification and Motivation

Design/Build/Fly (DBF) is an annual collegiate aircraft design competition hosted by the AIAA (American Institute of Aeronautics and Astronautics). Teams design, build, and fly a remotely piloted aircraft that must satisfy a mission-specific ruleset while balancing payload capability, aerodynamic performance, structural design, propulsion requirements, and mission completion time.
This project reuses the aircraft sizing and mission-performance model developed in Project 1. In Project 1, the design objective was to maximize a normalized competition score by varying aircraft geometry, payload weight, and mission cruise speeds. Project 2 studies the same design problem from a numerical optimization perspective by reformulating an aircraft-weight constraint with a quadratic penalty.
Aircraft design teams routinely optimize performance while enforcing structural and operational requirements. If those requirements are imposed using a poorly conditioned numerical formulation, an optimizer may require many more iterations or may fail to converge to a useful tolerance. This matters in DBF and in broader UAV preliminary design because design optimization is most useful when constraints can be enforced accurately without making the numerical problem unnecessarily difficult to solve.
The purpose of this project is therefore to demonstrate **intrinsic ill-conditioning** in a physically motivated aircraft optimization problem, quantify its effect on a baseline first-order optimizer, and demonstrate a more effective constraint-handling method.

The three scoring missions remain based on the Project 1 model:
- **Ground Mission:** The aircraft carries a sensor payload. The score contribution is proportional to sensor weight.
- **Mission 2:** The aircraft carries the sensor in its shipping container and completes five laps. The score rewards higher payload weight and shorter completion time.
- **Mission 3:** The sensor is deployed and the aircraft flies as many laps as possible within five minutes. The score is proportional to deployed sensor weight multiplied by the number of completed laps.

---
## 2. Formulation
### 2.1 Decision Variables

The optimization uses six continuous design variables:

$$
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
$$
| Variable | Definition | Project 2 bounds | Type |
|---|---|---:|---|
| $S_w$ | Wing area | $2 \le S_w \le 10\ \mathrm{ft}^2$ | Continuous |
| $L_F$ | Fuselage length | $5 \le L_F \le 7\ \mathrm{ft}$ | Continuous |
| $SW_2$ | Sensor weight for Ground Mission and Mission 2 | $4 \le SW_2 \le 16\ \mathrm{lb}$ | Continuous |
| $SW_3$ | Sensor-weight reduction for Mission 3 | $0 \le SW_3 \le 6\ \mathrm{lb}$ | Continuous |
| $V_2$ | Mission 2 cruise velocity | $40 \le V_2 \le 140\ \mathrm{ft/s}$ | Continuous |
| $V_3$ | Mission 3 cruise velocity | $40 \le V_3 \le 140\ \mathrm{ft/s}$ | Continuous |

The Project 1 upper bound on $SW_2$ was 12 lb. For this conditioning study, the upper bound is intentionally increased to 16 lb. This allows the 20 lb Mission 2 gross-weight constraint to become active before the payload box bound becomes active. This is a Project 2 study modification and is not intended to represent a change to the original competition rules.
All other aircraft geometry is derived from these variables. Wing span remains fixed at 6 ft, the sensor is modeled as a 6:1 cylindrical payload, and the fuselage and empennage dimensions are determined from the selected design variables.
### 2.2 Original Constrained Problem

The base objective is the negative normalized DBF competition score because the numerical optimizer minimizes its objective:

$$
f(\mathbf{x}) =
-\left[
\frac{SW_2}{12}
+
\frac{SW_2}{0.15\,t_2}
+
\frac{(SW_2-SW_3)\,\widetilde{L}_3}{50}
\right]
$$

where:
- $t_2$ is the time required to complete five Mission 2 laps,
- $\widetilde{L}_3$ is the continuous Mission 3 lap estimate used for the conditioning study,
- $\mathbf{x}$ is the six-dimensional design vector defined above.

The physical weight constraint used in this project is

$$
g(\mathbf{x}) =
W_{M2}(\mathbf{x}) - 20
\le 0
$$

with

$W_{\mathrm{M2}} = W_e + 1.1SW_2$
Here, $W_e$ is the modeled empty-aircraft weight. The factor 1.1 accounts for the modeled Mission 2 shipping-container weight used in Project 1.

The original constrained problem can therefore be written as

$$
\begin{aligned}
\min  \quad & f(\mathbf{x}) \\
\text{subject to} \quad & g(\mathbf{x}) \le 0, \\
& \mathbf{x}_{\min} \le \mathbf{x} \le \mathbf{x}_{\max}.
\end{aligned}
$$

### 2.3 Quadratic Penalty Reformulation

For the Project 2 conditioning study, the Mission 2 weight constraint is moved into the objective with a quadratic exterior penalty:

$F_{\rho}(\mathbf{x}) =
f(\mathbf{x})
+
\frac{\rho}{2}
\left[
\max\left(0,g(\mathbf{x})\right)
\right]^2$
The scalar $\rho>0$ is the penalty weight. A larger value of $\rho$ penalizes constraint violation more strongly, forcing the solution closer to the boundary $g(\mathbf{x})=0$. The same parameter also becomes the structural knob used to demonstrate ill-conditioning.
### 2.4 Constraints and Bounds
| Constraint or bound | Mathematical expression | Description |
|---|---|---|
| Mission 2 gross weight | $W_e+1.1SW_2 \le 20\ \mathrm{lb}$ | Team design target for aircraft and payload weight |
| Wing area | $2 \le S_w \le 10\ \mathrm{ft}^2$ | Limits lifting-surface size |
| Fuselage length | $5 \le L_F \le 7\ \mathrm{ft}$ | Limits aircraft length |
| Mission 2 payload | $4 \le SW_2 \le 16\ \mathrm{lb}$ | Project 2 study range |
| Mission 3 payload reduction | $0 \le SW_3 \le 6\ \mathrm{lb}$ | Limits payload reduction between missions |
| Mission 2 speed | $40 \le V_2 \le 140\ \mathrm{ft/s}$ | Cruise-speed range |
| Mission 3 speed | $40 \le V_3 \le 140\ \mathrm{ft/s}$ | Cruise-speed range |


### 2.5 Aircraft Model and Numerical Smoothing

The Project 2 model preserves the main engineering relationships from Project 1, including:
- sensor sizing from payload weight and lead-shot density,
- wing geometry based on a fixed 6 ft span,
- empirical fuselage, wing, electronics, and empennage weight estimates,
- induced and parasite drag,
- propulsion current as a function of required thrust and velocity,
- straight and turning flight times,
- battery-energy limitation,
- normalized DBF mission scoring.
Several numerical changes are required because Project 2 uses gradients, Hessians, and eigenvalue analysis. The original Project 1 simulator contained discrete switching operations that are appropriate for competition scoring but do not produce a smooth objective for local conditioning analysis.

1. Mission 3 lap count is treated continuously rather than applying `floor()`.
2. Throttle is solved continuously and is capped at 100%.
3. Requested speeds that cannot be sustained at full throttle are treated as physically infeasible through continuous propulsion-margin constraints.
4. Straight and turning propulsion feasibility are checked separately for Missions 2 and 3.
5. Mission 2 battery feasibility is imposed with a continuous energy-margin constraint for the required five laps.
6. Mission 3 remains continuously battery limited through a smooth minimum of the time-limited and energy-limited lap estimates.
7. Turn load factor is handled with a continuous $C_{L,\max}$-limited relation rather than discrete load-factor decrements.
8. Smooth limiting expressions are used where the original simulator contained abrupt switching operations.

When the design optima are located, the physical propulsion and battery constraints are enforced directly with SLSQP. The Mission 2 weight constraint remains in the fixed-penalty or augmented-Lagrangian formulation because that constraint is the deliberate source of the conditioning study.

### 2.6 Problem Classification

The Project 2 aircraft formulation is a **constrained, nonlinear, nonconvex optimization problem** studied through a smooth penalty reformulation.

The base objective is nonlinear because geometry, aerodynamic drag, propulsion demand, mission time, battery use, and competition score depend nonlinearly on the decision variables. The coupled aircraft-performance relationships can also produce multiple local optima, so the problem is nonconvex.
For the conditioning study, the Mission 2 inequality constraint is moved into the objective using a quadratic exterior penalty. This places the numerical mechanism in **Family G: penalty/barrier reformulations**.

---
## 3. Ill-Conditioning Mechanism
### 3.1 Family G Penalty Mechanism

When the weight constraint is active, the penalty contribution can be written locally as

$$
P(\mathbf{x}) =
\frac{\rho}{2}\,g(\mathbf{x})^2
$$

Its Hessian is

$$
\nabla^2 P =
\rho\,\nabla g\,\nabla g^{T}
+
\rho\,g\,\nabla^2 g
$$

Near the active constraint boundary, $g(\mathbf{x})\approx0$, so the second term becomes small and the dominant penalty curvature is approximately

$$\nabla^2P\approx\rho\,\nabla g\,\nabla g^T.$$
This term creates a stiff direction normal to the weight-constraint surface. Directions tangent to the constraint remain governed primarily by the original aircraft objective and therefore retain much smaller curvature.

The local Hessian condition number is

$$\kappa(H) =\frac{\lambda_{\max}(H)}{\lambda_{\min}(H)}$$
where $\lambda_{\max}$ and $\lambda_{\min}$ are the largest and smallest positive eigenvalues used in the local conditioning analysis. As $\rho$ increases, the penalty-dominated eigenvalue increases approximately in proportion to $\rho$, while the smaller-curvature directions change much less. The result is an increasingly elongated optimization landscape and a rapidly increasing condition number.
### 3.2 D1 - Hessian Eigenvalue Spectrum

The Hessian was evaluated numerically in normalized design coordinates near the optimum of the penalized problem. For the largest tested penalty weight, $\rho=10{,}000$, the smallest positive eigenvalue is approximately

$$
\lambda_{\min}=2.22\times10^{-2},
$$

and the largest is approximately

$$
\lambda_{\max}=1.97\times10^6.
$$

Therefore,

$$\kappa(H)=\frac{\lambda_{\max}}{\lambda_{\min}}\approx8.90\times10^7.$$
The spectrum spans many orders of magnitude, demonstrating a strongly elongated local optimization landscape.

![D1 Hessian eigenvalue spectrum](project2_outputs/D1_hessian_spectrum.png)
### 3.3 D2 - Intrinsic Ill-Conditioning Test

The required intrinsic-conditioning test has two parts:

1. The condition number must grow with a structural knob.
2. The large condition number must survive best per-coordinate Jacobi rescaling.

For this problem, the structural knob is the penalty weight $\rho$. The penalty weight was varied from $0.1$ to $10{,}000$. At each value, the penalized problem was optimized and the local Hessian was evaluated.

| Penalty weight $\rho$ | Weight residual $g$ (lb) | $\kappa(H)$ | $\kappa$ after Jacobi rescaling |
|---:|---:|---:|---:|
| 0.1 | $7.80\times10^{-1}$ | $1.48\times10^3$ | $5.22\times10^1$ |
| 1 | $9.72\times10^{-2}$ | $8.82\times10^3$ | $4.86\times10^1$ |
| 10 | $9.80\times10^{-3}$ | $8.90\times10^4$ | $5.02\times10^2$ |
| 100 | $9.84\times10^{-4}$ | $8.91\times10^5$ | $5.07\times10^3$ |
| 500 | $1.99\times10^{-4}$ | $4.46\times10^6$ | $2.54\times10^4$ |
| 1,000 | $9.83\times10^{-5}$ | $8.91\times10^6$ | $5.07\times10^4$ |
| 10,000 | $9.76\times10^{-6}$ | $8.90\times10^7$ | $5.07\times10^5$ |

The first requirement is satisfied because $\kappa(H)$ increases by approximately one order of magnitude each time $\rho$ increases by one order of magnitude. Over the tested range,

$$\kappa(H)\propto\rho$$

is a good approximation.

For the second requirement, let

$$D = \mathrm{diag}(H)$$

and apply symmetric Jacobi rescaling:

$$H_J =D^{-1/2} H D^{-1/2}$$

At $\rho=10{,}000$, the rescaled condition number is still approximately

$$\kappa(H_J)\approx5.07\times10^5.$$
Diagonal scaling reduces the numerical value of the condition number but does not remove its growth with $\rho$. Therefore, the problem passes both parts of the required intrinsic $\kappa$ test: the ill-conditioning grows with a structural parameter and survives per-coordinate rescaling.

![D2 condition number versus penalty weight](project2_outputs/D2_kappa_vs_rho.png)

---
## 4. Effect of Ill-Conditioning
### 4.1 D3 - Baseline Projected Gradient Descent

To demonstrate the practical effect of the increasing condition number, projected gradient descent was applied to the penalized objective. Before this experiment, each decision variable was normalized to the interval $[0,1]$ so that the measured slowdown would not be dominated by differences in engineering units.

The local fixed step size was chosen from the positive Hessian eigenvalues using

```math
\alpha = \frac{2}{L + \mu}
```
where
$$L=\lambda_{\max}(H)\qquad\text{and}\qquad\mu=\lambda_{\min}(H).$$

Convergence was declared when the projected-gradient norm satisfied

$$
\left\|\nabla_PF_{\rho}(\mathbf{x}_k)\right\|_2
\le
10^{-5}.
$$

Each case used a maximum of 15,000 iterations.
| Penalty weight $\rho$ | Iterations | Final projected-gradient norm | Final objective gap |
|---:|---:|---:|---:|
| 500 | 15,000* | $1.27\times10^{-1}$ | $5.58\times10^{-4}$ |
| 1,000 | 15,000* | $1.75\times10^{-1}$ | $4.31\times10^{-3}$ |
| 10,000 | 15,000* | $1.44\times10^{-1}$ | $1.04\times10^{-2}$ |

\*Iteration limit reached before satisfying the projected-gradient tolerance.

All three high-penalty cases reach the iteration budget without converging to the required tolerance. More importantly, the remaining objective gap after the same computational budget increases substantially as $\rho$ is increased. The $\rho=10{,}000$ case retains an objective gap almost twenty times larger than the $\rho=500$ case after the same 15,000 projected-gradient iterations.

The result demonstrates the practical implication of the D1/D2 conditioning analysis: as the penalty-induced Hessian becomes increasingly elongated, a fixed-step first-order method makes progressively less effective progress.

![D3 projected-gradient convergence](project2_outputs/D3_gradient_descent_convergence.png)

---
## 5. Proposed Solution and Demonstration
### 5.1 Augmented Lagrangian Method

A fixed quadratic penalty creates a numerical tradeoff. A small penalty weight produces a manageable optimization landscape but allows noticeable constraint violation. A very large penalty weight enforces the constraint accurately but makes the local problem severely ill-conditioned.

The proposed remedy is an **augmented Lagrangian method**. For the inequality constraint $g(\mathbf{x})\le0$, the implemented augmented objective is
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

After each inner optimization, the multiplier is updated as

```math
\lambda_{k+1}
=
\max\left(
0,\,
\lambda_k + \rho\,g(\mathbf{x}_k)
\right)
```
The multiplier carries information about the active constraint, allowing the method to obtain high constraint accuracy without increasing $\rho$ to the extremely large values required by a pure fixed-penalty method.

#### 5.1.1 Why the Augmented Lagrangian Was Chosen

Several numerical approaches can reduce the effect of ill-conditioning, but they do not all address the same underlying cause. In this problem, the dominant source of ill-conditioning is the large quadratic penalty coefficient required to enforce the active Mission 2 weight constraint. Therefore, the preferred remedy is one that reduces the need for a very large fixed penalty rather than only making the resulting ill-conditioned problem easier to solve.

| Method | Main idea | Advantage | Limitation for this problem |
|---|---|---|---|
| Variable scaling / Jacobi preconditioning | Rescale the design variables or Hessian so different directions have more comparable numerical magnitudes | Simple and inexpensive | D2 showed that the condition number remains large after Jacobi rescaling, so scaling does not remove the penalty-induced ill-conditioning |
| Preconditioned gradient or Newton-type methods | Modify the search direction using curvature information so the optimizer can move more efficiently through an elongated landscape | Can greatly improve convergence on an ill-conditioned objective | Improves the optimizer's response to the bad conditioning, but does not remove the large penalty curvature that caused the problem |
| Interior-point / barrier methods | Enforce feasibility by adding a barrier that prevents iterates from crossing the constraint boundary | Effective general-purpose constrained optimization method | Barrier parameters can also create increasingly large curvature near active constraints, so they can introduce a conditioning issue similar to the penalty mechanism being studied |
| Sequential Quadratic Programming (SQP) | Solve a sequence of local quadratic constrained subproblems | Very effective for smooth nonlinear constrained problems | Requires a more complex constrained subproblem framework and would change both the constraint-handling strategy and the optimization algorithm, making the D4 comparison less direct |
| Augmented Lagrangian | Combine a moderate quadratic penalty with a Lagrange multiplier that is updated using the observed constraint violation | Enforces the active constraint without requiring an extremely large penalty coefficient | Requires outer multiplier updates and inner optimization solves |

The augmented-Lagrangian method was selected because it addresses the specific mechanism identified in Section 3. The fixed quadratic penalty produces a Hessian contribution approximately proportional to

$$
\rho\,\nabla g\,\nabla g^T,
$$

so increasing $\rho$ directly creates the large-curvature direction responsible for the observed growth in $\kappa(H)$. The augmented Lagrangian introduces the multiplier $\lambda$, which carries information about the active constraint between outer iterations. This allows the constraint to be enforced accurately while keeping $\rho$ relatively small.

This distinction is important for the D4 comparison. A method such as preconditioning could make projected gradient descent perform better on the original poorly conditioned objective, but the underlying penalty formulation would remain poorly conditioned. The augmented Lagrangian instead changes the constraint-enforcement mechanism itself. Therefore, it provides a direct test of whether removing the need for a very large fixed penalty reduces the conditioning problem identified in D1-D3.

For this study, the augmented Lagrangian also provides a useful controlled comparison because projected gradient descent can still be applied to both the original fixed-penalty formulation and the augmented-Lagrangian local subproblem. This allows the D4 convergence comparison to isolate the effect of the formulation change rather than attributing the improvement to a completely different optimization algorithm.


### 5.2 D4 - Before/After Constraint Enforcement and Conditioning

The large fixed-penalty case used $\rho=10{,}000$ and produced

$$g(\mathbf{x})=9.76\times10^{-6}\ \mathrm{lb},$$

with a local Hessian condition number of

$$\kappa(H_{\mathrm{pen}})\approx8.90\times10^7.$$


The augmented Lagrangian started with $\rho=5$. By outer iteration 2,
corresponding to three augmented-Lagrangian subproblem solves, the weight
residual was reduced to approximately
$$g(\mathbf{x})=2.78\times10^{-7}\ \mathrm{lb},$$

while the local condition number is only

$$\kappa(H_{\mathrm{AL}})\approx8.28\times10^3.$$

Thus, the augmented-Lagrangian formulation obtains better weight-constraint accuracy while reducing the local Hessian condition number by approximately

$$\frac{8.90\times10^7}{8.28\times10^3}\approx1.07\times10^4.$$

The augmented-Lagrangian outer-iteration history is shown below.

![D4 augmented-Lagrangian constraint convergence](project2_outputs/D4_augmented_lagrangian_constraint.png)
### 5.3 D4 - Before/After First-Order Convergence

To satisfy the D4 requirement directly, the same projected-gradient method used for the baseline study is applied to both the large fixed-penalty formulation and the final augmented-Lagrangian local subproblem. Each run starts from the same normalized perturbation direction and magnitude around its corresponding local minimizer. For each formulation, the fixed step size is chosen from the local Hessian spectrum as

$$\alpha = \frac{2}{\lambda_{\max}+\lambda_{\min}}.$$

Both runs use the same projected-gradient tolerance, $10^{-6}$, and the same maximum budget of 20,000 iterations. This makes the comparison a direct test of how the change in conditioning affects the same first-order algorithm.

![D4 before/after convergence](project2_outputs/D4_before_after_convergence.png)

The corresponding numerical results are:

| Formulation | $\rho$ | $\kappa(H)$ | Iterations | Final projected-gradient norm | Best final absolute objective gap |
|---|---:|---:|---:|---:|---:|
| Fixed quadratic penalty | 10,000 | $8.90\times10^7$ | 20,000* | $3.07\times10^{-2}$ | $3.83\times10^{-4}$ |
| Augmented Lagrangian | 5 | $8.28\times10^3$ | 5,571 | $1.00\times10^{-5}$ | $\le10^{-14}$ plotting floor |

\*Maximum iteration budget reached.

The augmented-Lagrangian formulation reaches the prescribed
$10^{-5}$ projected-gradient tolerance after 5,571 iterations.
The fixed-penalty formulation does not reach the tolerance within the
20,000-iteration budget and finishes with a projected-gradient norm of
approximately $3.07\times10^{-2}$.

At the same time, the local Hessian condition number is reduced from
approximately $8.90\times10^7$ for the fixed penalty to
$8.28\times10^3$ for the augmented Lagrangian, a reduction of roughly
$1.07\times10^4$.

Thus, the same first-order method converges successfully on the
augmented-Lagrangian formulation while the highly ill-conditioned
fixed-penalty formulation remains unconverged after the larger
iteration budget.

![D4_weight_violation_fixed_vs_AL](project2_outputs/D4_weight_violation_fixed_vs_AL.png)

The results here show that for nearly every penalty parameter, the augmented-Langrangian formulation is able to enforce the weight constraint to a tolerable level for all values of $\rho$. The fixed penalty has the expected downward slope, which makes sense given the increasing penalty weight's implications on the optimization, however it further proves that the fixed penalty solution needs high penalty weight for constraint enforcement, putting the optimization under ill-conditioning. Overall, the plot supports that augmented-Langrangian formulation allows for sufficient constraints without the use of large penalty weights, which can impose ill-conditioning onto the optimization process. 

Together with the constraint-residual comparison in Section 5.2, this completes D4: the augmented Lagrangian achieves essentially the same engineering constraint accuracy, substantially lowers the local condition number, and produces a much faster convergence curve for the same first-order method.

The final augmented-Lagrangian design is approximately:
| Variable | Final value |
|---|---:|
| Wing area, $S_w$ | $5.62\ \mathrm{ft}^2$ |
| Fuselage length, $L_F$ | $5.21\ \mathrm{ft}$ |
| Mission 2 sensor weight, $SW_2$ | $12.55\ \mathrm{lb}$ |
| Mission 3 sensor reduction, $SW_3$ | $1.00\ \mathrm{lb}$ |
| Mission 2 cruise velocity, $V_2$ | $101.27\ \mathrm{ft/s}$ |
| Mission 3 cruise velocity, $V_3$ | $86.94\ \mathrm{ft/s}$ |
| Empty-aircraft weight | $6.19\ \mathrm{lb}$ |
| Mission 2 gross weight | $20.00\ \mathrm{lb}$ |
| Mission 2 five-lap time | $139.73\ \mathrm{s}$ |
| Continuous Mission 3 lap estimate | $7.21$ laps |

The final unpenalized objective is approximately

```math
f(\mathbf{x}) = -3.39336
```
### 5.4 Results and Interpretation

The penalty study demonstrates the expected Family G mechanism. Increasing the quadratic penalty weight reduces the weight-constraint violation, but the penalty also introduces a very stiff curvature direction normal to the constraint surface. The remaining directions retain substantially smaller curvature from the underlying aircraft-performance objective. This produces approximately linear growth of the condition number with $\rho$.
The Jacobi-rescaling test shows that the effect is not only a consequence of the variables having different units. Although diagonal rescaling improves the numerical condition number, the scaled problem still becomes increasingly ill-conditioned as $\rho$ grows.
The conditioning change has a direct effect on a baseline first-order optimizer. None of the $\rho=500$, $1{,}000$, or $10{,}000$ cases reach
the projected-gradient tolerance within the 15,000-iteration budget,
and the remaining objective gap grows substantially as the penalty
weight is increased.
The augmented-Lagrangian method addresses the mechanism directly. Instead of relying on an extremely large fixed penalty, it uses a moderate penalty together with a multiplier update. It reaches essentially the same constraint accuracy as the $\rho=10{,}000$ penalty case while keeping the local condition number far smaller and avoiding the large penalty curvature responsible for the observed first-order slowdown. The D4 same-method comparison confirms that this conditioning improvement translates into faster first-order progress: after 20,000 projected-gradient iterations, the augmented-Lagrangian objective gap is about $5.0\times10^4$ times smaller and its projected-gradient norm is about 17 times smaller than the fixed-penalty case.
For this DBF design problem, the result demonstrates that directly forcing engineering constraints with increasingly large penalty coefficients can make an otherwise manageable optimization problem unnecessarily difficult to solve. A constraint-handling approach such as an augmented Lagrangian can enforce the same engineering requirement while maintaining a better-conditioned numerical problem.

---
## 6. Assumptions and Simplifications

The results depend on the aircraft model and several simplifying assumptions:

- The model is a preliminary-design approximation and is not a full CFD, propulsion-map, structural, or flight-dynamics simulation.
- Wing span is fixed at 6 ft.
- The aircraft uses a single-wing fixed-wing configuration with an inverted T-tail.
- The sensor is represented as a 6:1 lead-shot-filled cylindrical payload with a nosecone.
- Structural weight is estimated using empirical relationships inherited from Project 1.
- Propulsion behavior is based on the Project 1 motor/propeller model.
- Propulsion throttle is capped at 100%; requested straight and turning speeds that cannot be sustained at full throttle are treated as infeasible.
- Mission 2 must complete five laps within 75% of the modeled 3300 mAh battery capacity, leaving a 25% reserve.
- Mission 3 lap count is treated continuously for the conditioning analysis rather than rounded to an integer.
- Mission 3 battery/time limiting behavior is represented with a smooth minimum to preserve differentiability.
- Numerical gradients and Hessians are computed using finite differences, so their accuracy depends on the selected perturbation size.
- Hessian conditioning is evaluated in normalized design coordinates to reduce contamination from engineering-unit scale differences.
- The 16 lb Project 2 payload upper bound is used only to expose the 20 lb Mission 2 gross-weight constraint as an active limiting mechanism.
- SLSQP is used to locate physically feasible reference optima while the weight constraint remains handled by the penalty or augmented-Lagrangian formulation.
- For the long D3 and D4 projected-gradient diagnostics, the physical propulsion and battery constraints are linearized locally about each reference optimum and projected using a low-cost half-space/box projection. These experiments therefore describe local first-order convergence behavior rather than a global feasible optimization trajectory.
- The augmented-Lagrangian and fixed-penalty comparisons describe local numerical behavior of this model and do not prove that the reported aircraft is the global optimum of the complete competition-design problem.
---
## 7. Code and Reproducibility

Two Python scripts are used for the Project 2 study:

- [`Project2Models.py`](Project2Models.py) contains the smooth DBF aircraft model, mission-performance calculations, propulsion model, continuous battery model, physical speed and battery margins, base competition objective, gross-weight constraint, quadratic penalty objective, and augmented-Lagrangian objective.
- [`Project2Diagnostics_D4.py`](Project2Diagnostics_D4.py) locates physically feasible design optima with SLSQP, computes finite-difference gradients and Hessians, evaluates Hessian spectra and condition numbers, performs Jacobi rescaling, carries out the D3/D4 local projected-gradient studies, and generates the augmented-Lagrangian and local-boundary diagnostics.

### 7.1 Software Requirements

The analysis uses Python 3 and the following third-party packages:

- NumPy
- SciPy
- Matplotlib

A typical installation command is:

```bash
python -m pip install numpy scipy matplotlib
```
### 7.2 Running the Analysis

Place the following files in the repository root:

```text
Project2DBFIllConditioning.md
Project2Models.py
Project2Diagnostics_D4.py
project2_outputs/
```

From the repository root, run:

```bash
python Project2Diagnostics_D4.py
```

The script recreates the numerical tables and figures in `project2_outputs/`.

The main generated files are:
```text
project2_outputs/D1_hessian_spectrum.png
project2_outputs/D2_kappa_vs_rho.png
project2_outputs/D3_gradient_descent_convergence.png
project2_outputs/D4_augmented_lagrangian_constraint.png
project2_outputs/D4_before_after_convergence.png
project2_outputs/D4_local_boundary_geometry.png
project2_outputs/conditioning_vs_rho.csv
project2_outputs/gradient_descent_summary.csv
project2_outputs/augmented_lagrangian_history.csv
project2_outputs/d4_before_after_summary.csv
```
### 7.3 Reproducibility Notes

No random sampling or randomized optimizer is used in the current analysis, so a random seed is not required. The starting designs, penalty weights, finite-difference step sizes, optimizer tolerances, and iteration limits are explicitly defined in the scripts. Running the diagnostic script with the same code and package versions should therefore reproduce the reported numerical results to normal floating-point tolerance.
As a small verification step, the reported condition-number trend, constraint residuals, and convergence counts should be checked against the generated CSV files after each run before the figures are committed to the repository.
