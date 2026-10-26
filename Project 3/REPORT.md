# Inverse oscillating boundary layer: a verified baseline and a diagnosed PINN limitation

**MAE 598 Project 3 | Technical experiment completed; researcher engagement and submission actions pending.**

Team members: Luc Andreassi, Micah Kawamura, Nicholas Cheng, Sivan Garteiz, Andrew Coverston. 
Researcher: Dr. Kiran Ramesh, PhD.

## Executive finding

For this one-parameter, single-frequency problem we recommend the conventional finite-difference inverse fit. The specialized PINN recovers viscosity and wall shear to accuracy comparable to the baseline, but its median training time is **88.2 times** the baseline's median inverse-fit time on the same machine (0.013 s versus 1.15 s). We do not claim that every PINN fails. We do show that a small neural training loss does not guarantee an accurate wall shear.

The thin-layer controls expose a concrete failure mode: **rescaling the network coordinate without sampling the near-wall physics adequately can give a very small training objective and an inaccurate wall derivative.** At delta/L = 0.01 and the original matched budget, adding near-wall collocation to the scaled network reduces the median held-out near-wall PDE residual (RMS) by **19.2 times** and the median wall-shear amplitude error by **3.75 times**. Training ten times longer then exposes residual overfitting between the fixed collocation points, including in the improved model. A residual-adaptive refinement, which adds collocation points where the current PDE residual is largest, reduces the median held-out near-wall PDE residual by **12 times** relative to the longer-trained vanilla PINN. Sections 7–8 present the controlled Path B diagnosis, the initial remedy, the counterexample and the final Path C refinement, with unfavorable results kept.

This is an exploratory three-seed study. No researcher has supplied an accuracy or cost target. All numerical findings below are measured; the interview, researcher review and communication status are reported as they are and are not inferred.

## 1. Outreach log and research problem brief
### 1.1 Outreach log

### 1.1 Outreach log

Researchers contacted, in chronological order (all dates 2026). Copies of the outreach emails: [team to attach].

| Date | Researcher | Action | Outcome |
|---|---|---|---|
| 8 Oct | Prof. Yulia Peet | Email requesting a research interview | No response |
| 12 Oct | Prof. Jay Oswald | Email requesting a research interview | Replied 22 Oct, after the team had committed to a project; not pursued |
| 17 Oct | Dr. Kiran Ramesh, Associate Professor of Aerospace Engineering, SEMTE, Arizona State University | Email requesting a research interview | Accepted; interview scheduled for 18 Oct |
| 18 Oct | Dr. Ramesh | First interview (Zoom, 6 PM), prepared with the team's interview guide (`Interview_Question.md`) and a review of his publications (`ramesh-research-review.md`). Micah Kawamura could not attend | Research context and candidate problems discussed |
| 19 Oct | Dr. Ramesh, with Nicholas Cheng, Sivan Garteiz and Luc Andreassi | Second interview (Zoom) to choose the problem | Project secured: oscillating boundary-layer inverse problem with a PINN |

**Consent.** Dr. Ramesh consented during the Zoom meetings to being recorded and to being named and reported on in this project.

**Problem brief review.** Pending: the brief in Section 1.2 was sent for his correction on [date]; response: [pending].

**Recommendation memo.** Draft in Section 9; not yet sent.

**Interview findings.** The points below come mainly from the 19 Oct interview; notes from the 18 Oct interview: [team to add].

- *Research objective.* Dr. Ramesh is developing an unsteady version of XFOIL. Like XFOIL, it uses viscous–inviscid interaction: the problem is split into a potential-flow solution and an integral boundary-layer solution.
- *Where the difficulty lies.* Integral boundary-layer methods collapse the wall-normal dimension, which leaves more unknowns than equations, so closure models must supply the missing information. Extending these methods to unsteady flow requires new unsteady integral boundary-layer formulations and new unsteady closure models. This is his current research focus.
- *Data and solvers.* He is using OpenFOAM to generate CFD data for building unsteady closure models. That work is at an early stage and outside this project's scope.
- *Role of this problem.* He described the oscillating boundary-layer inverse problem (estimating viscosity and wall shear stress from noisy sensor data with a PINN) as a demonstration: it has an exact solution, so it validates the method rather than serving as a production application.
- *Not established in the interviews.* A required wall-shear accuracy, the solver's run cost and hardware, and the group's prior ML or PINN experience. [team to add if these were discussed]




### 1.2 Provisional problem brief for researcher review

Oscillating external flow creates an unsteady viscous layer near an airfoil or rotor-blade surface. Wall shear can lead the external velocity, so a closure based only on instantaneous outer velocity misses part of the dynamics. The supplied task card connects this minimal example to unsteady integral boundary-layer modeling in Kiran's lab. This experiment asks whether sparse noisy velocities can identify viscosity and recover the unsteady wall-shear history.

The minimal case is laminar, one-dimensional in wall-normal position, periodic and driven at one frequency. It retains a thin near-wall region and the wall-shear phase relation, while removing streamwise variation, airfoil geometry, separation and nonlinear unsteady closure physics. We compare a neural inverse approach with a conventional solver-based inverse fit. The numerical baseline's cost is measured here; the researcher's production solver, actual computational bottleneck, number of cases and previous PINN attempts remain unknown until the interview evidence is added.

**Questions to resolve in the researcher check:** Is viscosity identification representative of the information gap in the lab's closure work? What wall-shear amplitude/phase errors are useful? Are frequency content, probe placement, separation, noisy forcing or model discrepancy the next important difficulty?

## 2. Complete mathematical formulation

### 2.1 Dimensional problem

Let y be distance from a stationary wall. The imposed outer velocity is Ue(t)=U0 cos(omega t). Incompressible, unidirectional laminar flow satisfies

$$
u_t=\frac{dU_e}{dt}+\nu u_{yy},\qquad u(0,t)=0,\qquad u(y,t)\to U_e(t)\text{ as }y\to\infty.
$$

The pressure gradient is represented by dUe/dt. Viscosity nu is constant but unknown to the inverse methods. Inertia associated with streamwise variation is absent in this minimal problem.

### 2.2 Nondimensional variables, domain and unknown

Use yhat=y/L, that=omega t, uhat=u/U0 and kappa=nu/(omega L^2). Hats are omitted below for readability. The computational domain is y in [0,1] and t in [0,2 pi), representing one complete periodic cycle:

$$
u_t=-\sin t+\kappa u_{yy},\qquad u(0,t)=0,\quad u(1,t)=\cos t,\quad u(y,t+2\pi)=u(y,t).
$$

The inverse unknown is the scalar kappa>0. Recover dimensional viscosity by nu=kappa omega L^2 if dimensional frequency and length are known. We study the periodic state after startup, so no independent initial condition is imposed on the PINN. Both inverse methods restrict kappa to [1e-8,0.2]; these broad bounds are shared and contain every tested truth.

### 2.3 Exact reference, wall shear and finite-domain consistency

The half-space reference supplied in the task card is

$$
u_\infty(y,t)=\cos t-e^{-y/\delta}\cos(t-y/\delta),\qquad \delta=\sqrt{2\kappa}.
$$

Here delta means delta/L. Its dimensionless wall shear, a velocity gradient rather than dimensional stress, is

$$
g(t)=u_y(0,t)=\frac{\sqrt{2}}{\delta}\cos(t+\pi/4).
$$

The target shear amplitude is sqrt(2)/delta and the phase lead is 45 degrees. Dimensional stress would be mu U0/L times this gradient. If viscosity is inferred, dimensional stress requires consistently converting the inferred nu into mu=rho nu.

The half-space field is not exactly equal to cos(t) at y=1. Its discrepancy is bounded by exp(-1/delta): 0.00248 U0 at delta=1/6, 4.54e-5 U0 at delta=0.1 and negligible in the thinnest cases. We retain the task-card half-space reference for data generation and scoring, while enforcing its approximate far-field condition in both inverse methods. We separately verify discretization against the **boundary-consistent finite-domain analytic amplitude**

$$
a_F(y)=1-\frac{\sinh[q(1-y)]}{\sinh q},\qquad q=(1+i)/\delta,\qquad u_F=\Re[a_F e^{it}].
$$

The implementation uses a stable exponential form to avoid overflow. This finite-domain expression is used for solver verification only, not for fitting measurements. The finite-domain approximation is an explicitly bounded model discrepancy, not a coding error or PINN pathology.

### 2.4 Difficulty knobs and experimental design

Primary knob: delta in {1/6,1/10,1/20,1/50,1/100}; smaller delta means a thinner layer. Equivalently use L/delta in {6,10,20,50,100}, which increases with difficulty. At fixed dimensional nu and L, increasing omega creates the thinner layers. The corresponding kappa truths are delta^2/2.

Secondary knobs: noise sigma in {0,0.01,0.02,0.05,0.10} U0 and Nt in {64,16,8} samples per probe per cycle. We vary one knob at a time about delta=0.1, sigma=0.02, Nt=16. Removing duplicate reference settings gives 11 settings; seeds 17,29,41 give 33 measurement datasets. There are 117 primary method/case runs, including paired diagnosis controls, plus nine longer-budget controls and 15 adaptive continuation fits. These supplementary tests reuse the existing thin-case data. This is not a full factorial sweep and does not characterize every interaction among knobs.

The three probes are at y={0.5 delta,delta,2 delta}. Their positions move with the true layer between synthetic experiments, following the task card. Within an inverse fit the supplied physical probe positions are fixed and do not move with candidate kappa. This deliberately maintains measurement coverage inside the layer; it does not test the harder scenario of fixed instruments falling outside a thinning layer.

## 3. Reproducibility and implementation map

### 3.1 Run instructions

Requires Python 3.12 and the pinned packages in requirements.txt. No PyTorch, GPU, XFOIL or experimental data are required.

```bash
python -m venv .venv
# Activate .venv for your operating system, then:
python -m pip install -r requirements.txt
python -m unittest -v test_numerics.py
python oscillating_bl.py --suite full --seeds 17 29 41 --width 16 --collocation 160 --max-nfev 600 --threads 1 --output final_results
python budget_controls.py --results final_results
python adaptive_fix.py --results final_results --rounds 5
python audit_submission.py --results final_results
python make_report.py --results final_results
```

For a two-setting smoke run, use --suite pilot --seeds 17 --output pilot_results. For an interrupted full run, repeat the identical command with --resume. Resume refuses changes to code hash or experiment settings. A resumed run retains original per-model training times, while environment.json records the resumed invocation's wall duration.

The provided figures and numerical tables come from the completed full command, not the smoke run. Near-identical floating-point results are expected with the pinned versions; wall times and exact optimizer trajectories can differ across machines/BLAS implementations.

### 3.2 Main functions and data flow

| Stage | Function(s) | Inputs -> outputs |
|---|---|---|
| Reference and noise | exact_amplitude, exact_velocity, make_measurements | Truth, probe positions, seed -> clean/noisy samples |
| Numerical forward solver | grid_operator, fd_periodic | Candidate kappa -> numerical periodic velocity amplitude |
| Time-integration verification | mol_check | Numerical periodic state -> BDF integrated cycle |
| Conventional inverse estimate | fit_baseline | Fixed noisy probes -> fitted kappa |
| Neural representation | Jet2, neural_fields | Trainable weights and y -> A,B and spatial derivatives |
| Physics training | residual_jacobian, train_pinn | Data and collocation -> neural weights and kappa |
| Scoring | evaluate, field_metrics | Trained method plus withheld exact reference -> errors and timing |
| Diagnostics | p1_verification, plot_summary, plot_case, plot_diagnosis | Saved results -> P1--P5 figures |
| Report | make_report.py | Completed numerical CSV/JSON files -> this report |

### 3.3 Files and audit trail

- oscillating_bl.py: complete numerical study and figure generation.
- test_numerics.py: six scientific verification tests.
- requirements.txt: exact package versions.
- make_report.py: report generator; tables use measured results.
- budget_controls.py: nine tenfold-budget controls, retained separately from the primary runs.
- adaptive_fix.py: five fixed residual-adaptive refinement rounds for each thin-case seed.
- audit_submission.py: independent field-error audit and supplementary figures.
- final_results/data/: all 33 clean/noisy measurement datasets, standard Gaussian draws and probe coordinates.
- final_results/models/: every fitted model and per-run metrics/optimizer status.
- final_results/histories/: objective evaluations, elapsed time and separate losses for every neural run.
- final_results/metrics.csv: all 117 method/case records, including all seeds and negative outcomes.
- final_results/budget_metrics.csv and adaptive_metrics.csv: longer-budget controls and every refinement round.
- final_results/baseline_grid_convergence.csv and mol_verification.json: independent baseline evidence.
- final_results/fields/: representative full fields used for the displayed moderate/thin-case plots; all fields can be regenerated from the code.
- final_results/figures/: diagnostic PNGs linked below.
- final_results/environment.json: command, versions, hardware, thread settings and source hash.

The archive includes every model, measurement dataset, loss history and summary metric. To keep it compact, it includes full field arrays for the six plotted baseline/vanilla/combined runs (first seed, delta=0.1 and 0.01) and the first-seed final adaptive field; the full commands regenerate every other field array. No input PDF is redistributed.

Source SHA-256 at experiment execution: `567777beb993db101a9ade366ad1180b6f76d36a02b9561fd4324d60a24d7f20`.

## 4. P1: conventional forward and inverse baseline

### 4.1 Spatial discretization

Use a fixed grid yj=sinh(6 j/N)/sinh(6), with N=400 intervals for inverse fitting. The grid is chosen without knowing the true kappa. For hm=yj-yj-1 and hp=yj+1-yj,

$$
(D_2u)_j=\frac{2u_{j-1}}{h_m(h_m+h_p)}-\frac{2u_j}{h_m h_p}+\frac{2u_{j+1}}{h_p(h_m+h_p)}.
$$

Boundary values are incorporated in the right-hand side. The smooth stretched grid clusters points at the wall while retaining the full domain.

### 4.2 Periodic solver and time-integration check

Because the PDE and boundary conditions are linear with one known forcing frequency, insert u=Re[a(y) exp(it)] into the discretized PDE:

$$
(\kappa D_2-iI)a=-i\mathbf1,
$$

with a(0)=0 and a(1)=1. A tridiagonal complex banded solve produces the numerical periodic state. It does not insert the exact spatial formula. This is an efficient conventional solver for exactly the PDE being compared, not a weakened baseline.

**Departure from the task-card time-marching recipe:** the main inverse baseline solves the periodic finite-difference system directly instead of marching from u=1 for each trial viscosity. A separate BDF solve_ivp check integrates that same method-of-lines system for a complete cycle, initialized from the numerical periodic state. The maximum measured relative difference from that state is **1.37e-07**. This confirms the time and harmonic solvers agree; it is not a startup-from-uniform-flow validation. The task card's suggested 10--20 startup cycles are not a convergence guarantee: the slowest error mode decays approximately as exp(-kappa pi^2 t). At delta=0.01, a factor 1e-6 decay can require about 4,456 cycles. Solving the periodic state avoids charging an unnecessarily expensive transient to the baseline.

Grid refinement N={100,200,400,800} is checked against the finite-domain exact amplitude. At N=400 the maximum tested relative L2 amplitude error is **1.3e-05**. The independent unit test confirms approximately second-order convergence.

![P1 grid convergence](final_results/figures/P1_grid_convergence.png)

### 4.3 Conventional inverse fit

For candidate kappa, solve the finite-difference system, interpolate its complex amplitude to the fixed probes, and minimize the mean squared velocity mismatch to the same noisy measurements used by the PINN. Optimize log(kappa) over the shared broad bounds. A 61-point coarse scan precedes three local bounded scalar refinements. The global scan is necessary because trial layers much thinner than all probes produce almost flat misfit regions; a single broad local scalar search can stop there incorrectly.

The baseline's measured inverse time includes the scan, every numerical forward solve, interpolation and scalar optimization. Exact reference values never enter its objective. The independent noiseless inverse test recovers kappa within 0.1%; the actual default-grid error is much smaller.

## 5. Synthetic measurements and the PINN

### 5.1 Measurement generation

At each of the three heights, sample one cycle at t=2 pi j/Nt, j=0,...,Nt-1. The duplicate cycle endpoint is omitted. Generate

$$
u^{obs}_{ij}=u_\infty(y_i,t_j;\kappa_{true})+\sigma\epsilon_{ij},\qquad\epsilon_{ij}\sim\mathcal N(0,1).
$$

Sigma is a fraction of U0, not a percentage of the instantaneous velocity. Zero crossings therefore do not create a relative-noise singularity. Noise is added only to observations. The reference, PDE and boundary conditions remain clean. Save the clean samples and standard Gaussian draws. The same dataset is reused for every paired method. Seeds also change neural initialization; these three replications therefore combine measurement and initialization variability rather than independently estimating each contribution.

### 5.2 Trainable neural representation and its scope

The network is a fully trainable, single-hidden-layer, 16-unit tanh MLP with two outputs:

$$
z_j(y)=\tanh(w_j y/s+b_j),\quad A(y)=\sum_jv^A_jz_j+c_A,\quad B(y)=\sum_jv^B_jz_j+c_B,
$$

$$u_\theta(y,t)=A(y)\cos t+B(y)\sin t.$$

All hidden weights, hidden biases, output weights and output biases are optimized together with log(kappa): 67 trainable parameters. Time periodicity is exact. Neither A(y) nor B(y) is supplied as the Stokes exponential; their spatial shapes must be learned from sparse probes and physics.

**Departure from a generic MLP taking (y,cos t,sin t):** the harmonic output form uses the known single-frequency linear response to restrict time dependence. This is a specialized PINN and gives it useful physics information beyond a generic space-time network. Both neural variants and the strong numerical baseline exploit the same known forcing frequency. Conclusions do not establish how a generic PyTorch space-time PINN would train, or how nonlinear/multifrequency airfoil flows would behave.

Spatial first/second derivatives are propagated with a minimal forward-mode automatic differentiator (Jet2) using affine and tanh primitive rules. Optimizer parameter Jacobians are analytic. Complex-step checks independently verify every Jacobian column and the second spatial derivatives in both raw/scaled coordinates. No finite differences of noisy observations are used by the PINN.

### 5.3 Loss, optimizer and initialization

For the harmonic trial function the residual is

$$
r(y,t)=[B-\kappa A'']\cos t+[1-A-\kappa B'']\sin t.
$$

The PDE loss integrates its squared value **exactly in time**, giving

$$
L_{PDE}=\frac1{2N_c}\sum_{j=1}^{N_c}\left[(B_j-\kappa A_j'')^2+(1-A_j-\kappa B_j'')^2\right].
$$

The data loss is ordinary mean squared velocity error. The boundary loss is the mean of A(0)^2, B(0)^2, [A(1)-1]^2 and B(1)^2. Train with L=Ldata+LPDE+LBC, equal component weights and no separate initial-condition term. Boundary conditions are soft rather than embedded exactly in the spatial network.

SciPy's trust-region reflective nonlinear least-squares method uses the supplied full Jacobian, float64 arithmetic, ftol=xtol=gtol=1e-9 and at most 600 function evaluations. This is an optimization choice suited to a small residual-based inverse problem; it is not an Adam experiment. Wall time is also reported because equal evaluation counts do not imply identical compute cost.

Initialization is away from truth: seeds 17 and 41 start kappa at 5 times true, and seed 29 at 0.2 times true, following the task card's synthetic initialization check. This is a controlled benchmark initialization, not a procedure for selecting a starting viscosity from unknown real data. Each matched method uses the same initial weights and estimate. An additional experiment with data-selected initialization would be needed for a deployment claim.

**Stopping status:** 84 of 84 neural runs reached the 600-evaluation budget rather than a convergence tolerance. Their measured errors are valid finite-budget results; no claim of fully converged global optima is made. Status, final objective and optimality are retained per run.

## 6. P2--P5 diagnostics and measured results

### 6.1 P2: separate loss histories

Every neural run records data, PDE and boundary losses separately on each objective evaluation. Histories include rejected trust-region trial steps, so curves need not be monotone. The supervised oracle control's optimized history is its supervised field loss, not a physics-training history; its final physics residual is scored separately.

![P2 moderate case](final_results/figures/P2_losses_d0.1.png)

![P2 thin case](final_results/figures/P2_losses_d0.01.png)

### 6.2 P3: errors against each knob

Relative field L2 is the square root of integrated mean-in-time squared error divided by integrated mean-in-time squared reference, evaluated on an independent uniform grid with 1,201 heights and 128 times. Also score the region y<=min(6 delta,1) separately: a whole-domain norm can conceal errors in a very thin layer. Parameter error is abs(kappa_hat/kappa_true-1). Shear amplitude and phase are derived from the learned wall derivative, not imposed as targets. Three-seed medians and full observed ranges are shown; the ranges are not confidence intervals.

![P3 knob sweeps](final_results/figures/P3_knob_sweeps.png)

P3 shows the main physical difficulty knob, \(L/\delta\), in the first column. Increasing this ratio produces a thinner boundary layer with sharper velocity changes near the wall, which can make the flow more difficult to resolve. The plots compare the conventional numerical baseline (black), vanilla PINN (red), and combined PINN (green). The combined PINN uses coordinate scaling and concentrates more physics collocation points near the wall.

Two additional knobs control the measurement data: the standard deviation of the simulated noise and the number of samples per cycle per probe. Higher noise makes the observations less reliable, while fewer samples provide less information for estimating viscosity and reconstructing the flow. Each knob is varied separately while the other settings remain fixed.

The results show that errors increase with measurement noise and generally decrease with more samples. However, the thickness sweep does not show a uniform increase in error as the boundary layer becomes thinner. Whole-field velocity error decreases partly because the near-wall region occupies a smaller fraction of the domain, potentially concealing localized errors.

The PINNs produce viscosity and wall-shear estimates comparable to the baseline, with occasional improvements, but the baseline generally achieves lower whole-field velocity error. The combined PINN does not show a consistent accuracy advantage. All three methods remain sensitive to measurement noise, indicating that improved near-wall physics enforcement alone does not eliminate uncertainty caused by noisy observations.

| delta/L | kappa true | Baseline kappa error (%) | Vanilla kappa error (%) | Combined kappa error (%) | Vanilla near-wall L2 (%) | Combined near-wall L2 (%) |
|---|---:|---:|---:|---:|---:|---:|
| 0.167 | 0.0139 | 2.13 | 2.2 | 2.2 | 0.355 | 0.355 |
| 0.1 | 0.005 | 2.12 | 2.16 | 2.19 | 0.422 | 0.421 |
| 0.05 | 0.00125 | 2.12 | 2.06 | 2.16 | 0.502 | 0.476 |
| 0.02 | 0.0002 | 2.12 | 1.56 | 2.1 | 0.578 | 0.505 |
| 0.01 | 5e-05 | 2.12 | 1.88 | 2.01 | 0.611 | 0.605 |

### 6.3 Thin-case numerical comparison

At delta=0.01, sigma=0.02 and Nt=16, the following values are medians over seeds 17,29,41. Controls marked supervised and known_kappa receive additional oracle information and **do not estimate viscosity**.

| Method | Viscosity error (%) | Whole-field L2 (%) | Near-wall L2 (%) | Shear amplitude error (%) | Shear phase error (deg) | Fit/train (s) |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 2.12 | 0.0759 | 0.325 | 1.08 | 0.00237 | 0.0131 |
| vanilla | 1.88 | 0.144 | 0.611 | 0.571 | 0.833 | 1.09 |
| sampling | 2.1 | 0.121 | 0.517 | 1.18 | 0.332 | 1.14 |
| scaling | 10.2 | 0.192 | 0.82 | 4.32 | 1.68 | 1.07 |
| combined | 2.01 | 0.142 | 0.605 | 1.15 | 0.322 | 1.17 |
| supervised | -- | 0.0163 | 0.0288 | 0.282 | 0.894 | 1.31 |
| known_kappa | -- | 0.145 | 0.619 | 0.0252 | 0.343 | 1.02 |

### 6.4 P4: where the errors occur and their spectra

The following plots use seed 17, chosen as the first configured seed before inspecting errors. Field panels have separately labeled color scales to reveal each method's structure, not to claim equal magnitudes. The tables and spectra provide quantitative comparisons. Baseline continuous residuals are not inferred by differentiating interpolated FD values; that panel states it is unavailable. The FD system itself is solved to numerical precision.

![P4 moderate near-wall fields](final_results/figures/P4_fields_d0.1.png)

![P4 thin near-wall fields](final_results/figures/P4_fields_d0.01.png)

![P4 thin shear, hysteresis and spatial spectrum](final_results/figures/P4_shear_spectrum_d0.01.png)

The spatial error spectrum is a type-II discrete cosine transform, averaged over time. It is a Fourier cosine spectrum appropriate to a nonperiodic wall-normal interval; a periodic spatial FFT would introduce artificial endpoint jumps. The temporal FFT is also saved and plotted:

![P4 temporal error spectrum](final_results/figures/P4_temporal_spectrum_d0.01.png)

The first temporal harmonic dominates by construction. This experiment cannot diagnose temporal causality violations or broadband temporal spectral bias. Spatial high-mode error is descriptive evidence, not proof that spectral bias alone caused a failure. A changing wall derivative can matter even when the integrated velocity error is small.

### 6.5 P5: total inverse cost, query cost and hardware

Across the full primary sweeps, median inverse fit/train times are 0.013 s for the conventional baseline, 1.15 s for the vanilla PINN, and 1.15 s for the combined PINN. Each plotted point represents the median timing across the three seeds at that boundary-layer thickness, so the curves emphasize the typical computational cost at each difficulty setting rather than the timing of an individual run. The vanilla PINN therefore requires 88.2 times longer than the baseline to obtain an inverse solution. In the left panel below, the x-axis is the inverse boundary-layer thickness, $L/\delta$, so moving to the right represents an increasingly thin boundary layer. Fit/train time remains nearly constant as the layer becomes thinner for all three methods. The dominant cost difference is therefore the choice of solution method rather than boundary-layer thickness: both PINN approaches consistently require roughly two orders of magnitude more computation than the conventional inverse fit. The logarithmic y-axis in the left panel makes the nearly two-order-of-magnitude difference between the baseline and PINN fitting times easier to compare; the vertical separation therefore represents a substantial computational-cost difference rather than a small absolute offset.

The right panel measures inference time after each method has already been fitted or trained. Each timing is the median of nine warmed evaluations of the same 1,201 by 128 velocity grid, corresponding to 153,728 predicted velocity values. Evaluation is fast for all three methods, but the conventional baseline remains faster at approximately 0.7–1.0 ms, compared with roughly 2–3 ms for the PINNs. As with fitting cost, there is no clear increase in inference time as the boundary layer becomes thinner. Because these measurements are only a few milliseconds, the small point-to-point fluctuations should not be interpreted as a meaningful trend. Derivatives are not requested during this query, although the current neural implementation still computes spatial derivative intermediates, so the PINN timings should not be interpreted as an optimized deployment benchmark.

![P5 costs](final_results/figures/P5_cost.png)

Lower values indicate lower computational cost in both panels. The left panel compares the cost of obtaining each inverse solution, while the right panel compares the cost of evaluating the resulting fitted models.

Hardware/software: **AMD64 Family 25 Model 33 Stepping 0, AuthenticAMD / Windows-11-10.0.26200-SP0**; one BLAS thread; Python **3.12.1**; NumPy **1.26.4**; SciPy **1.14.1**; Matplotlib **3.9.3**. No GPU is used. See environment.json for full provenance and linked BLAS information. Timing values should be rerun on the intended deployment hardware.

The surrogate break-even query number $N^*$ is not applicable to this experiment because these models perform inverse estimation rather than learning an operator that can be reused across a family of PDE inputs. A new measurement dataset requires a new conventional inverse fit or PINN training run. Although the fitted fields from either method can subsequently be evaluated cheaply, the right panel shows no inference-time advantage that could offset the PINNs' much larger fitting cost. For this minimal case, the appropriate comparison is therefore the total cost of obtaining the inverse solution together with its resulting accuracy, rather than a train-once, query-many break-even calculation.

## 7. Path B: controlled diagnosis

### 7.1 Observed limitation, without overstating it

The vanilla harmonic PINN does not catastrophically collapse over the prescribed thickness sweep. It often obtains useful viscosity estimates. It remains more expensive than the verified conventional fit, and its thin-layer field/derivative errors are not controlled by training loss alone. A plausible coordinate-scaling-only treatment is particularly revealing: it can have small sampled PDE and boundary losses while its held-out near-wall PDE residual and wall-shear error are substantially larger. This negative outcome is retained rather than hidden.

### 7.2 Candidate explanations and interventions

| Candidate explanation | Controlled test | What it establishes and what it does not |
|---|---|---|
| Neural architecture cannot represent the spatial solution | Fit the same raw-coordinate architecture directly to dense exact A,B targets | Low errors demonstrate practical representability under this budget; the oracle fit receives extra information and is not an inverse competitor |
| Unknown viscosity alone creates the field problem | Keep kappa fixed to truth, retaining sparse noisy data and raw-coordinate physics training | Separates parameter-estimation difficulty from field-learning difficulty; oracle control is not a valid viscosity estimator |
| Scaling alone resolves thin spatial features | Compare vanilla and scaling with matched uniform collocation | Scaling alone can worsen derivative/generalization behavior; rescaling is not automatically a remedy |
| Near-wall physics is insufficiently constrained by the chosen sampling | Compare scaling vs combined with identical weights, scale, data, optimizer, parameter count and evaluation budget | Adding near-wall samples reduces held-out residual/shear errors; this supports a sampling/physics-enforcement mechanism rather than insufficient probe information alone |
| Merely increasing model capacity or training longer explains improvement | Keep width=16, count=160 and max_nfev=600 fixed in the sampling comparison | Those quantities cannot explain this paired improvement; optimizer trajectories still change |

The sampling intervention changes both where residual constraints are evaluated and the effective spatial weighting of the empirical physics loss. It does **not** separately identify coverage and loss reweighting as exclusive causes. Consequently the supported diagnosis is inadequate **near-wall physics enforcement**, not a proof of pure spectral bias or of one unique root cause.

| Seed | Vanilla kappa error (%) | Combined kappa error (%) | Scaling-only shear error (%) | Combined shear error (%) | Scaling-only layer PDE RMS | Combined layer PDE RMS |
|---|---:|---:|---:|---:|---:|---:|
| 17 | 6.07 | 3.89 | 3.88 | 1.64 | 0.0378 | 0.0035 |
| 29 | 1.82 | 2.01 | 12.8 | 1.15 | 0.165 | 0.00834 |
| 41 | 1.88 | 1.05 | 4.32 | 0.206 | 0.101 | 0.00526 |

![B/C controlled experiments](final_results/figures/BC_controlled_diagnosis.png)

The three paired seeds show that near-wall sampling improves the scaled network's held-out layer residual and wall-shear amplitude error in every tested replication. The combined method does not improve every metric of the raw vanilla network for every seed. Observation noise, initialization and limited optimization remain possible contributors. Three seeds are insufficient for a broad statistical superiority claim.

## 8. Path C: initial remedy, falsification check and residual-adaptive refinement

### 8.1 Initial matched-budget remedy

Define an initial layer estimate delta0=sqrt(2 kappa0). The combined remedy uses y/delta0 as the network coordinate and moves half the stratified collocation points into [0,min(6 delta0,1)], leaving the other half distributed across the whole domain. The number of points remains 160. Delta0 is computed from the **initial estimate**, not from true thickness or post-hoc error maps. The network coordinate scale is fixed during that run, avoiding an extra moving-coordinate optimization effect.

The mechanism is direct: coordinate scaling supplies useful spatial scales, while near-wall residual points constrain the derivatives at those scales. The far field remains represented by the uniform half. No analytic exponential spatial basis or fitted reference solution is inserted. The remedy is inexpensive to implement and the all-seed before/after evidence appears in the paired table and diagnosis plot.

Relative to scaling alone, median layer PDE RMS improves by 19.2 times and shear amplitude error by 3.75 times at the original budget. This is a measured improvement, but a further check is needed before calling it robust.

### 8.2 Longer training: a useful counterexample

We retrain vanilla, scaling and combined variants from their same initial states with up to 6,000 evaluations, keeping the same 160 fixed points. Training longer does not remove the held-out residual problem; some runs terminate on a tight objective tolerance yet have worse residuals between collocation points. This separates the hypothesis of merely stopping too early from finite-point enforcement: a better sampled objective is not necessarily a better solution. The longer runs are all retained in budget_metrics.csv; no checkpoint is selected using exact errors.

### 8.3 Final refinement: use PDE residuals to add constraints

Starting from each original combined model, evaluate the current PDE residual on a fixed monitor consisting of 4,097 uniform wall-normal locations and 4,097 locations spanning [0,min(6 delta0,1)], with duplicates removed. True solution values and true kappa are not used in monitoring or point selection. At each of five fixed rounds, add 64 separated high-residual points while retaining all existing constraints, then train for at most 600 more evaluations. The final allocation is 480 points and at most 3,600 total evaluations including initial training. Every round is saved; intermediate reversals are shown, and the final fifth round is reported rather than the best truth-scored checkpoint.

| Method / allocation | Max evaluations | Collocation points | Kappa error (%) | Shear amplitude error (%) | Held-out layer PDE RMS | Fit/refinement time (s) |
|---|---:|---:|---:|---:|---:|---:|
| vanilla | 600 | 160 | 1.88 | 0.571 | 0.00857 | 1.09 |
| combined | 600 | 160 | 2.01 | 1.15 | 0.00526 | 1.17 |
| vanilla, longer | 6000 | 160 | 3.79 | 4.34 | 0.0384 | 10.7 |
| scaling, longer | 6000 | 160 | 10.8 | 4.67 | 0.109 | 3.07 |
| combined, longer | 6000 | 160 | 9.01 | 1.15 | 0.0918 | 7.62 |
| residual-adaptive, final round | 3600 | 480 | 2.2 | 1.21 | 0.0032 | 10.2 |

![Longer training and adaptive refinement](final_results/figures/BC_budget_and_adaptive.png)

Compared with the longer-trained vanilla control, the final adaptive model reduces median held-out layer PDE RMS by 12 times. All three seeds improve that residual relative to their longer-trained vanilla counterparts. The final refinement is not uniformly better than the original 600-evaluation combined model on every seed/metric, and residuals remain nonzero. The increased point count and monitoring also cost work; this is a targeted remedy, not an equal-compute superiority claim.

An independent residual-resolution audit uses 24,001 heights inside y in [0,0.06], without changing training. It confirms an adaptive improvement factor of **11.4 times** over the longer-trained vanilla model. This finer check guards against the coarser plotted grid missing narrow neural residual spikes; per-seed values are in dense_residual_audit.csv.

![Adaptive histories and every refinement round](final_results/figures/BC_adaptive_losses.png)

Losses across refinement boundaries use a changing point set and should not be read as a single unchanged quadrature objective. Adaptive timing includes the initial fit, residual monitoring, point selection and refinement optimization, but excludes subsequent scoring and plotting. In the following first-seed spatial profiles, residuals are evaluated at independent plotted locations, not just at training points.

![Near-wall error and residual before/after](final_results/figures/BC_near_wall_profiles.png)

### 8.4 Practical conclusion

The supported mechanism is **finite-point near-wall physics underconstraint and residual overfitting**, with optimization and noise also contributing. A fixed allocation can improve the original training budget yet become insufficient under longer training. Residual-adaptive constraint enrichment addresses that mechanism. It does not establish that a PINN is preferable to the conventional inverse fit, nor a unique root cause across every architecture. Further tests should independently vary quadrature weighting, coverage, model width and fixed probe placement.

## 9. Recommendation memo (draft; send after team review)

**To:** [Confirmed researcher; supplied task card refers to Kiran's lab]  
**From:** [Team names], MAE 598  
**Subject:** Sparse-probe inverse oscillating boundary layer: retain the conventional fit; enforce near-wall physics in neural trials

We tested the single-frequency oscillating boundary layer as a minimal example retaining unsteady wall-shear phase and a thin viscous layer. Synthetic probes at three heights were generated from the known reference, with controlled noise and sampling. A verified finite-difference inverse fit and a 67-parameter harmonic PINN used the same observations.

For this one-parameter problem, we recommend the conventional inverse fit. The PINN generally recovered viscosity, but its median training cost was 88.2 times the baseline's measured inverse cost, without a consistent accuracy advantage. This result does not rule out neural methods for unknown closures, nonlinear unsteady airfoils or other inverse problems.

The neural controls identified an actionable limitation. Small sampled physics losses can conceal poor derivatives between collocation points. Near-wall allocation initially improved the scaled network, but longer fixed-point training exposed residual overfitting. A further five-round residual-adaptive refinement, using only the current PDE residual to add points, reduced median held-out layer PDE RMS by 12 times relative to the longer-trained vanilla control at delta/L=0.01. We recommend monitoring held-out near-wall physics and enriching constraints rather than trusting training loss or simply training longer.

Confidence is strongest in the numerical verification and the paired intervention within this small study. Generalization remains uncertain: only three seeds, one forcing harmonic, laminar periodic flow and probes deliberately placed inside each layer were tested. Many neural optimizations exhausted their evaluation budget. No laboratory-specific accuracy requirement or production computational bottleneck was available to us.

The next useful experiment is a fixed-probe placement sweep across frequency, followed by a two-harmonic forcing case. Please confirm which wall-shear amplitude/phase errors and missing-physics scenarios would change your modeling decision. If the actual goal is learning an unsteady closure rather than estimating viscosity, we suggest defining that gap before scaling up the neural model.

**Delivery status:** draft prepared, not sent. Add actual delivery date and researcher response after sending it.

## 10. Assumptions, limitations and interpretation

1. Synthetic observations satisfy the half-space model apart from independent additive Gaussian noise. Real sensor biases, correlated noise and model discrepancy are not represented.
2. The finite top boundary is approximate; its mismatch is quantified above and matters most at delta=1/6.
3. Harmonic time structure is supplied to both methods. No claim is made about generic space-time PINNs, nonlinear flow, multiple frequencies, turbulence, transition, separation or airfoil geometry.
4. True thickness determines synthetic probe locations and the controlled initial-estimate factors. The inverse algorithms receive only fixed probe coordinates, observations and starting estimates, not the true kappa in their objectives. Oracle controls are explicitly labeled.
5. A one-at-a-time sweep does not characterize interactions such as very thin layers plus 10% noise plus only eight samples.
6. Noise and initialization vary together across three seeds. Full ranges are descriptive, not statistical confidence bounds.
7. Equal network size and evaluation budgets allow controlled comparisons; many finite-budget fits are not optimizer-converged.
8. Spatial spectra cannot by themselves identify spectral bias. Sampling interventions change coverage and effective residual weighting together. The adaptive remedy is demonstrated only at the thinnest case; it has not been swept over all knobs.
9. The optional unknown-outer-flow/two-harmonic extension is not implemented. The one-page recommendation prioritizes the next informative experiment over adding an unvalidated extension.
10. The interview, researcher review, public-sharing consent, repository publication, memo delivery and Canvas submission remain team actions. This report does not invent them.

## 11. Verification and a hand-checkable example

All six tests in test_numerics.py pass: exact boundary/PDE consistency; full network parameter Jacobian by complex step; second spatial derivatives by complex step; second-order FD convergence; noiseless inverse recovery; and reproducible noise/probe coordinates. BDF/periodic-state checks and grid-refinement results are provided separately. Verify again after modifying the code.

At delta=0.1, kappa=0.005, y=delta/2=0.05 and t=0,

$$u=1-e^{-0.5}\cos(0.5)=0.4677192697843292.$$

The dimensionless shear amplitude is sqrt(2)/0.1=14.142135623730951 with phase lead 45 degrees. At delta=0.01 the amplitude is 141.4213562373095. These are independent arithmetic checks of the reference and the expected scaling.

The field files provide the exact and predicted values used for error scoring. The packaged audit below independently recomputes one full relative L2 error from a saved field without calling the neural evaluator or report generator.

## 12. Rubric mapping and final team actions

| Rubric category | Points | Evidence/location | Status |
|---|---:|---|---|
| Outreach/interview | 10 | Section 1.1 and actual correspondence | Team evidence required |
| Problem brief | 15 | Section 1.2 | Provisional; researcher check required |
| Minimal formulation | 15 | Section 2 | Implemented; simplifications stated |
| Baseline/diagnostics | 20 | Sections 4--6 and saved data | Executed P1--P5 |
| Path findings | 25 | Sections 7--8 and paired controls | Measured diagnosis and targeted remedy; limits stated |
| Recommendation memo | 5 | Section 9 | Draft complete; actual sending pending |
| Reproducibility | 5 | Section 3, pinned dependencies, code/data/models | Executed and checked |
| Presentation/clarity | 5 | Single Markdown report, equations and figures | Check rendering in the final GitHub repository |
| Path C bonus | up to 10 | Section 8 | Targeted before/after evidence; award is grader's decision |

Before submitting: replace team/researcher placeholders, add the actual outreach/interview findings and consent, obtain a researcher correction of the brief, send the reviewed memo, and record any reply. Put this Markdown file and its relative figure paths in the team's public GitHub repository with runnable code, keeping restricted interview material out of it. Submit the repository URL on Canvas. No repository was created or message sent by this task.

## References

1. Supplied task card, *Oscillating boundary layer inverse problem with a PINN*, class-pinn-oscillating-bl.pdf, three pages. Source of the physical example, parameter/noise/sample grids and requested comparisons. Not redistributed in the archive.
2. [MAE 598 Project 3 assignment and rubric](https://designinformaticslab.github.io/DesignOptimization2025/project3.html), accessed 2026-10-19. Source of report structure, P1--P5 and Path B/C expectations.
3. [SciPy 1.14.1 least_squares documentation](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.optimize.least_squares.html). Trust-region reflective nonlinear least squares and analytic Jacobian interface.
4. [SciPy solve_ivp documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). BDF integration and sparse Jacobian support; executed version is pinned in the environment.
5. [SciPy 1.14.1 solve_banded documentation](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.linalg.solve_banded.html). Numerical tridiagonal harmonic solver.
6. [Krishnapriyan et al. (2021), *Characterizing possible failure modes in physics-informed neural networks*](https://arxiv.org/abs/2109.01050). Motivation for separating representability and physics-training issues; this study's numerical conclusions come from its own saved experiments.
