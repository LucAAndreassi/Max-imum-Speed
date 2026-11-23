---
layout: default
title: Carbon-fiber coaxial drone arm - verified preliminary design study
---

# A uniform carbon-fiber arm for coaxial rotors

**Working engineering report, 23 November 2026.** Core numerical verification: **@@COUNTS@@ checks passed**. The current uniform-tube optimizer uses an explicitly estimated effective material profile. Its safety-factor claim applies to nominal material-strength checks; shell failure, joints, laminate interaction and dynamic stability remain unverified.

@@OPTIMIZATION_HEADLINE@@

All numerical tables and the viewer are generated from `assets/results.json`. The legacy 30/40 mm stock comparison is retained below as source-backed preliminary analysis. Both stock examples have 1.5 mm walls and **do not meet the user's updated wall > 2 mm requirement**.

<div id="arm-viewer"></div>
<script type="module" src="assets/viewer.js"></script>

## 1. Design brief and evidence

The team specified a uniform hollow circular carbon-fiber tube, 2 ft long, supporting 30 lbf **combined** coaxial thrust and a @@TIP_MASS_LB@@ lb combined motor assembly. Both rotors operate at 3,500 RPM with two blades each. Applied torque is zero. Thrust rotates with the motor mount. The objective is minimum total weight while maintaining safety factor 2, with **wall thickness strictly greater than 2 mm**. The arm is represented as horizontal with an ideal fixed root. The listed @@TIP_MASS_LB@@ lb is used once as the tip mass; any uncounted propeller, mount, wiring, or fitting mass must be added.

Using SI units, $L=0.6096\,\mathrm{m}$, $T=133.446648\,\mathrm{N}$, $m_t=@@TIP_MASS_KG@@\,\mathrm{kg}$ and $g=9.80665\,\mathrm{m/s^2}$. Thrust alone produces a root bending moment of **@@ROOT_MOMENT@@ N m**, despite zero applied torque. The motor weight and tube self-weight oppose the upward thrust. [NIST conversion reference](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8).

The shaft frequency is **@@SHAFT@@ Hz** and each rotor's two-blade passing frequency is **@@BPF@@ Hz**. Two coaxial rotors do not automatically double either frequency; their relative phase and aerodynamic interaction determine additional forcing. The supplied `optimize_arm.py` explicitly requires the first bending frequency to reach the 1x shaft rate; the integrated study retains this as an editable rotor multiple. No additional frequency separation margin or deflection limit is assumed. Meeting an unloaded frequency threshold does not establish operating dynamic stability.

### Proposed commercial material

The Easy Composites CFT-WF family uses roll-wrapped 0/90 reinforcement with a woven outer finish. Its datasheet reports **typical** longitudinal tensile modulus 64 GPa and tensile strength 650 MPa. We use that modulus as a preliminary effective bending modulus. This is an assumption about $EI$, not a ply-level laminate calculation. The sheet explicitly disallows treating its typical values as specification limits. [CFT-WF datasheet](https://media.easycomposites.co.uk/datasheets/EC-TDS-Carbon-Fibre-Woven-Finish-Tube.pdf).

The reference is 40 mm OD / 37 mm ID, published line mass 275.6 g/m; the lighter comparison is 30 mm OD / 27 mm ID, 193.8 g/m. Both have 1.5 mm nominal walls. Published line masses are used directly, without inventing one common density. [40 mm product](https://www.easycomposites.co.uk/40mm-woven-finish-carbon-fibre-tube), [30 mm product](https://www.easycomposites.co.uk/30mm-woven-finish-carbon-fibre-tube).

Comparable roll-wrapped products provide **comparison data**. ACP publishes compressive strength 132 ksi (about 910 MPa), but its modulus and laminate differ. Graphtek's GR-CFRWT sheet lists compression 910.1 MPa and torsional strength 103.4 MPa. Torsional strength is not automatically a transverse or interlaminar shear allowable. ACP's web page prints a strength in modulus units for its torsional entry; that ambiguous entry was excluded. The user authorized estimated properties for optimization: `optimizer_inputs.json` records their transfer assumptions explicitly, while `inputs.json` preserves the original material-data gaps. [ACP technical data](https://acpcomposites.com/shop/tubes-rods/roll-wrapped-tubing/carbon-fiber-uni-roll-wrapped-tubes), [Graphtek sheet](https://www.graphitestore.com/core/media/media.nl?c=4343521&h=22DqRe5XcXtOFJUqMtPHOuVeC4o8PE9P-gnJuKBYWPrJ3COD&id=9929).

## 2. Formulation and unresolved requirements

For a uniform tube the physical sizing variables are $d=(D,t)$, with inner diameter $D_i=D-2t$ and $D>2t>0$. The section formulas are

$$A=\frac{\pi}{4}(D^2-D_i^2),\qquad I=\frac{\pi}{64}(D^4-D_i^4).$$

For stock tubes the objective is $\min_j g(m_t+\mu_j L)$ over an explicitly documented catalogue. For a custom uniform laminate with validated density and manufacturing bounds it becomes $\min_{D,t} g[m_t+\rho A(D,t)L]$. These are different optimization problems: a continuously optimized dimension is not automatically an available stock tube.

The intended nominal strength constraints are

$$g_t=2\sigma_{t,\max}/X_t-1\le0,\qquad g_c=2|\sigma_{c,\max}|/X_c-1\le0,\qquad g_s=2\tau_{\max}/S-1\le0.$$

Here $X_t,X_c,S$ must be applicable laminate/tube allowables. This maximum-stress screen alone does not cover interaction, delamination, crushing at clamps, ovalization, or local buckling. A full safety factor is the smallest valid reserve factor across the relevant failure modes. Von Mises yield stress is not a carbon-fiber failure criterion. The reported $650\,\mathrm{MPa}/\sigma_b$ is only a typical tensile-property ratio.

Configurable constraints are $|w(L)|/w_{allow}-1\le0$, $f_{min}/f_1-1\le0$ and $|z_{NL}(L)|/w_{allow}-1\le0$. The integrated default study enables $f_{min}=\mathrm{RPM}/60$ from the supplied script; displacement limits remain unset. No limits are inferred from the bridge example. Explicit preliminary geometry bounds are recorded in the optimization section. Using twice the loads and also halving the allowable would double-count the requested factor.

Diameter and wall bounds are mathematically necessary for a strength-only continuous study. In the thin-wall limit, $A\approx\pi Dt$ and the bending section modulus is $Z\approx\pi D^2t/4$. Holding bending strength fixed permits $t\propto D^{-2}$, so tube mass decreases as $D^{-1}$. Without a diameter bound, minimum wall, restricted stock catalogue, or a validated local-failure constraint, that model has no finite minimum-mass tube. A numerical optimizer reaching an arbitrary search boundary would not resolve this modeling issue.

The [course assignment](https://designinformaticslab.github.io/DesignOptimization2025/project4.html) expects a substantial differentiable design problem and a nonlinear-design comparison. The chosen uniform tube has **two** physical sizing variables. The team should confirm this scope with the instructor; adding artificial independent dimensions to the mesh would contradict the uniform-tube brief. The current work includes a two-variable constrained optimizer, KKT checks and a with/without-N2 comparison. The general 3D-frame and laminate-analysis expectations still need scope agreement.

## 3. Linear analysis model

In the loaded bending plane, the Euler-Bernoulli equation is $(EIw'')''=q$, where $q=-\mu g$. The root has $w(0)=w'(0)=0$. The initial net transverse tip force is $P=T-m_tg$ and the applied tip bending moment is zero. With $+z$ upward,

$$w(L)=\frac{PL^3}{3EI}+\frac{qL^4}{8EI},\qquad M(x)=P(L-x)+\frac{q}{2}(L-x)^2,\qquad \sigma_b(x)=\frac{|M(x)|D}{2I}.$$

`arm.py` assembles an original sparse cubic-Hermite bending FE solver, with consistent distributed loading and consistent translational mass, using $[w,w']$ at each node. The production linear mesh contains 20 elements. The fixed root is removed from the solve; support reactions are recovered from the unmodified system. A point tip mass enters the translational mass matrix exactly once. The independent FE modal comparison solves $K\phi=\omega^2M\phi$. Default vibration frequencies and shapes now use the uniform cantilever with a tip mass in Yntema, NACA TN 3459, printed page 17, equations (13)-(14), implemented in `vibration.py`.

Define $r=m_t/(\mu L)$ and take each positive root $\beta_n$ of

$$H(\beta,r)=1+\cos\beta\cosh\beta-r\beta(\sin\beta\cosh\beta-\cos\beta\sinh\beta)=0.$$

Then $f_n=\beta_n^2\sqrt{EI/(\mu L^4)}/(2\pi)$. Equation (14) supplies mode shapes and slopes sampled at the FE nodes; normalization uses $\int_0^L\mu y_n^2\,dx+m_t y_n(L)^2=1$. These samples are continuum shapes, not exact eigenvectors of the discrete matrices. The implementation divides $H$ by $\cosh\beta$ for stable root finding. This is the paper's **nonrotating** beam case: the arm is fixed, while motor RPM sets only the sizing target. No rotating-blade centrifugal stiffening is inferred. [NACA TN 3459, equations (13)-(14)](https://ntrs.nasa.gov/api/citations/19930084187/downloads/19930084187.pdf#page=21).

This reduced model is suitable for preliminary bending of a straight circular tube under planar loads. It does not include axial stretch, transverse shear deformation, torsional modes, rotor inertia, gyroscopic effects, clamp compliance, laminate coupling, or shell deformation. The displayed frequency is a **bending frequency**, not a verified lowest frequency of every possible aircraft-arm mode. Out-of-plane bending is degenerate only under the assumed circumferentially uniform effective stiffness. Zero net torque does not remove possible torsional dynamics.

The homogeneous-section shear screen uses $\tau_{NA}=V(R^3-r_i^3)/(3It)$, derived from $VQ/(Ib)$ with $Q=2(R^3-r_i^3)/3$ and $b=2t$. It is reported to expose a missing failure check, not as a validated composite shear prediction.

## 4. Legacy stock-tube comparison (1.5 mm walls, excluded from new design)

@@RESULTS@@

Dimensions and material values are nominal. No tolerance, damage, fatigue, environmental or joint knockdown has been applied. The candidates were chosen to compare two sourced uniform products; the catalogue has not been exhaustively searched. Neither is designated a final selection.

The comparison compression ratios use the separately sourced 910.1 MPa value only as an explicit estimate. The required-strength rows invert the SF 2 screen: they state what applicable material strengths would have to exceed, before other failure modes are considered. They do not establish those strengths for the selected tube. This makes new test data easy to assess without changing the structural calculation.

## 5. Numerical verification

`verify.py` uses deliberately synthetic beam properties: $L=1$ m, $E=70$ GPa, $D=40$ mm, $t=2$ mm, $\mu=0.4$ kg/m and a 100 N tip-force fixture. These are numerical fixtures, not claimed carbon-fiber material data. Closed-form targets are separate from the FE assembly. The nonlinear reference is an independently implemented continuum boundary-value solve; it is a numerical reference, not an exact closed-form value.

For a cantilever, integrating the bending equation gives $w(L)=PL^3/(3EI)$ and $w'(L)=PL^2/(2EI)$; uniform load gives $qL^4/(8EI)$. The simply supported uniform-load midpoint is $5qL^4/(384EI)$. [University mechanics derivation](https://mechanics.ju.se/SolidMechanics/EnergyMethods.html).

The zero-tip-mass limit recovers $1+\cos\beta\cosh\beta=0$. The independent FE eigenproblem is compared with the analytical roots both with and without a tip mass, including the first three modes of the tip-mass fixture. The FE modal refinement plot retains the matrix eigensolver so that it measures discretization error; the default analytical frequencies are independent of mesh size.

@@VERIFICATION@@

The free-free check correctly has two rigid modes in the scalar bending subspace; it is not claimed to demonstrate six rigid modes of a 3D frame. This model does not use coordinate transformation matrices, so a 3D rotation-invariance test remains necessary if the solver is extended to a 3D frame. A team member still needs to derive or independently review the benchmark targets, as the assignment asks; these agent-written checks do not substitute for that human review.

![Mesh convergence](assets/mesh-convergence.png)

The uniform-load displacement is quartic, so the cubic FE field has nonzero interpolation error even when nodal deflections are exact. The observed field convergence is approximately fourth order. The modal refinement also approaches fourth order; neither check relies on a tip value that is already exact.

## 6. Adjoint and eigenvalue sensitivities

With $h(d,u)=K(d)u-F(d)=0$, any scalar output $Q(d,u)$ has

$$K^Ta=Q_u^T,\qquad \frac{dQ}{dd_j}=Q_{d_j}-a^T(K_{d_j}u-F_{d_j}).$$

The tip-displacement output selects the tip translation. The root-stress output also differentiates its section and element-force dependence. Complex-step calculations supply partial derivatives with the state held fixed; the state dependence is handled by the transposed linear solve, not by differencing converged states. The original solver fixtures hold line mass fixed. The new optimizer instead uses $\mu=\rho(t)A(D,t)$, differentiating mass, self-weight, inertia, section stiffness and any tabulated properties. Its separate gradient checks and unit tests cover those dependencies, including nonlinear adjoints and optional frequency/deflection constraints.

The default vibration sensitivities differentiate the characteristic root implicitly: $\partial\beta/\partial r=-H_r/H_\beta$. Hence

$$d\ln f=\tfrac12d\ln(EI)-\tfrac12d\ln\mu-2d\ln L+\frac{2}{\beta}\frac{\partial\beta}{\partial r}\,dr,\qquad r=\frac{m_t}{\mu L}.$$

The optimizer includes both section/material stiffness changes and the density-dependent line-mass changes in this expression. Central-difference checks include tabulated thickness-dependent modulus and density. For comparison, a mass-normalized FE bending eigenvector obeys

$$\frac{d\omega^2}{dd_j}=\phi^T(K_{d_j}-\omega^2 M_{d_j})\phi,\qquad \frac{df}{dd_j}=\frac{1}{8\pi^2 f}\frac{d\omega^2}{dd_j}.$$

The check sweeps central-difference step sizes in two coordinate directions and three seeded random directions for displacement, stress, and frequency. Coarse-step errors fall approximately quadratically; cancellation appears at small steps, especially for eigenvalues. The table records the best step for each direction, while the plot exposes the whole sweep.

![Directional gradient checks](assets/gradient-check.png)

For the synthetic eight-element, two-variable fixture, measured tip-gradient time is **@@ADJOINT_MS@@ ms** with the adjoint and **@@FD_MS@@ ms** with central differences. This small problem does not demonstrate the large-variable advantage required in a substantial adjoint study. Timings are machine-dependent and are regenerated on each run.

## 7. N2: large rotation with motor-following thrust

Let $s$ be undeformed arc length, $\theta(s)$ the centerline angle, $x'=\cos\theta$ and $z'=\sin\theta$. The thrust vector is $T[-\sin\theta_L,\cos\theta_L]$, normal to the tip tangent. Gravity remains vertical. The inextensible planar equilibrium is

$$ (EI\theta')'+T\cos(\theta_L-\theta)+[-m_tg-\mu g(L-s)]\cos\theta=0,\quad \theta(0)=0,\quad EI\theta'(L)=0.$$

Linear angle elements and four-point Gauss integration produce the residual. Newton iteration uses load stepping, residual line search and the consistent **nonsymmetric** tangent. A follower force has no general potential energy; a symmetric energy Hessian would omit its dependence on the unknown tip angle. The nonlinear adjoint therefore solves $R_\theta^Ta=Q_\theta^T$, then $dQ/dd=Q_d-a^TR_d$. The ladder checks this transposed solve and the full follower tangent.

A pure end-moment fixture has the exact circular solution $x(L)=L\sin\alpha/\alpha$, $z(L)=L(1-\cos\alpha)/\alpha$, $\alpha=ML/EI$. Dead transverse-load and follower-load solutions are also compared with separate continuum BVP implementations. The [Abaqus nonlinear cantilever benchmarks](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-nlgeocantilever.htm) motivate the dead-load and pure-moment checks; they do not validate this code's follower-force implementation by themselves.

@@REFINEMENT@@

![Nonlinear follower-force gradient check](assets/nonlinear-gradient.png)

The nonlinear EI-sensitivity step sweep checks the transposed tangent adjoint against separate converged solves. The optimization section adds complete OD/wall geometry-variable checks with density-dependent mass and loading.

The results table compares the linear solution with the continuum nonlinear reference to separate physical nonlinearity from the lower-order angle-element error. The final nonlinear FE mesh has 160 elements. The exploratory path below multiplies thrust **and both gravity contributions** together; it is a numerical load-factor experiment, not an RPM/throttle prediction or a proof of strength at increased loads.

The 3D linear state is sampled with cubic Hermite interpolation including nodal slopes. The nonlinear state is drawn from the 160-element FE angle field, integrated within each element; the independent BVP curve is retained separately in the results file. The table distinguishes the continuum reference from the FE values.

![Nonlinear load path](assets/nonlinear-path.png)

Nonconservative follower-force systems can require a dynamic flutter analysis. A converged static solution does not prove dynamic stability. Local shell buckling and Brazier ovalization are also absent. The next section reports optimized designs with and without the follower-force constraints under the stated estimated material and preliminary bounds.

## 8. Constrained mass optimization and KKT checks

@@OPTIMIZATION@@

## 9. Limitations and work remaining

- Obtain applicable compression, shear, local shell-failure and joint data; retain comparison estimates as estimates.
- Set stiffness/deflection and resonance-separation requirements; confirm provisional OD and upper-wall bounds and actual manufacturing tolerances.
- Confirm the two-variable uniform-tube scope and reduced bending model with the instructor. The present model is not a general 3D frame.
- Check rotational/translational inertia and actual support stiffness. Verify the included @@TIP_MASS_LB@@ lb covers every carried item.
- Extend the preliminary nominal-strength optimum to a validated laminate/shell/joint design before treating it as a safe physical arm.
- Publish the team repository and enable GitHub Pages. No Git remote or deployment was provided, so no public submission has been made.

## 10. Reproduction and AI use

Run `python optimize_arm.py --update-report` in an environment installed from `requirements.txt` to rerun sizing and regenerate the report. The optimizer verifies its gradients before using them. The report step runs the core verification ladder and refuses regeneration on failure. Optimizer results are imported only when source/input fingerprints match. Edit `report_template.md` for prose; `index.md` is generated. Edit `inputs.json` for the original loads and `optimizer_inputs.json` for optimizer bounds, estimated properties, and optional limits. `README.md` gives Windows commands. `python -m unittest test_optimize test_vibration -v` runs eleven optimizer tests and seven dedicated vibration tests, including limiting cases, four-mode FE comparisons, continuum mass orthogonality and implicit-root derivatives.

The coding agent researched the linked sources, wrote the reduced FE solver, nonlinear solver, sensitivity routines, checks, viewer and report generator. The team supplied the engineering brief and clarified follower thrust. No human derivation or physical validation is claimed. A real defect caught during verification was a floating-point-only nonlinear tangent accumulator: complex-step perturbations caused a casting failure. The accumulator now preserves the perturbation dtype and the nonlinear tangent and adjoint checks pass. The dead-load formulation was also revised when the team clarified that thrust follows the motor mount; the old assumption is not used for the drone results.
