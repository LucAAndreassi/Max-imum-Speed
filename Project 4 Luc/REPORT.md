# PDE-Constrained Optimization of a Uniform Carbon-Fiber Coaxial Drone Arm

## 1. Summary

This project minimizes the mass of a **0.6096 m (2 ft) uniform carbon-fiber drone arm** that supports a coaxial two-motor propulsion assembly. The arm carries **30 lbf (133.4466 N) combined thrust**, a **1.5 lb (0.68039 kg) tip mass**, and operates at **3500 RPM**. The structural model uses an in-house Euler-Bernoulli finite-element solver for linear bending, an analytical/FE eigenvalue model for vibration with tip mass, and an **N2 finite-rotation follower-thrust model** solved with Newton iteration and adjoint sensitivities.

The optimized uniform tube is **55.2848 mm outer diameter with a 2.010 mm wall**, giving a structural tube mass of **311.714 g**. Its first bending frequency is **58.3333 Hz**, exactly at the configured requirement of one shaft-order frequency, so **vibration is the governing performance constraint**. The wall lower bound is also active. The modeled nominal material-strength safety factor is **36.81**, and the linear and nonlinear tip deflections are approximately **1.240 mm**. Adding the N2 follower-load constraints changes the optimized mass by **0.0%**, showing that geometric nonlinearity is correctly modeled but inactive for this load case.

The original 40 mm OD / 2.5 mm-wall starting design is lighter (**272.904 g**) but is **infeasible** because its first bending frequency is only **38.714 Hz**. A more meaningful comparison is the lightest feasible fixed-wall reference, **52.3 mm OD / 2.5 mm wall at 362.416 g**. Optimizing both diameter and wall thickness reduces tube mass by **13.990%** relative to that feasible reference.

<div id="arm-viewer"></div>
<script type="module" src="assets/viewer.js"></script>

The live viewer compares the starting and optimized tubes using one linked camera, true hollow-section geometry, deformation control, follower-load state display, first-bending-mode coloring and animation, and element inspection with nominal stress information. The metrics panel reports the mass objective and every enabled constraint with PASS/FAIL status. The live viewer draws its numerical content from the repository's generated assets/results.json; the polished report text and reported numerical results are maintained separately from the automated numerical-generation workflow.

The principal engineering result is that satisfying static strength alone is not sufficient for this arm. The structure can remain well below its material-strength limit while still being dynamically unsuitable because of insufficient vibrational performance.

All code and results can be found in the Github directory: https://github.com/LucAAndreassi/Max-imum-Speed/tree/main/Project%204%20Luc

## 2. Structure and formulation

### 2.1 Engineering motivation and physical model

The structure is one arm of a large coaxial multirotor. Reducing structural mass is beneficial because it increases useful payload capacity and reduces vehicle inertial mass; however, removing material also reduces bending stiffness and can move the arm’s natural frequencies closer to rotor-induced excitation. The arm is modeled as horizontal and perfectly clamped at the vehicle frame. The tip carries the complete **1.5 lb combined motor assembly**. Two coaxial rotors together generate **30 lbf total thrust**. Because the coaxial rotors counter-rotate, the mean applied shaft-reaction torque is assumed to cancel and is therefore set to zero in the present structural model. In the nonlinear model, thrust follows the orientation of the motor mount and therefore rotates with the deformed tip tangent.

The shaft frequency is

$$
f_{shaft}=\frac{3500}{60}=58.3333\ \text{Hz},
$$

and the two-blade passing frequency of each rotor is **116.667 Hz**. The vibration limit is derived from the specified 3500 RPM operating speed using the current design criterion that the first bending frequency must be at least one shaft-order frequency.

$$
f_1\ge f_{shaft}=58.3333\ \text{Hz}.
$$

This is a preliminary sizing threshold, not a complete dynamic-stability or resonance-separation analysis.

### 2.2 Effective carbon-fiber material model

The optimization uses a preliminary effective roll-wrapped CFRP property set:

| Property | Value used | Basis recorded in repository |
|---|---:|---|
| Effective longitudinal/bending modulus $E$ | 64 GPa | Easy Composites CFT-WF typical longitudinal tensile modulus |
| Density $\rho$ | 1520 $kg/m^3$ | Graphtek comparison laminate, used as estimate |
| Tensile strength $X_t$ | 650 MPa | Easy Composites typical tensile strength |
| Compressive strength $X_c$ | 910.1 MPa | Graphtek comparison laminate, used as estimate |
| Shear strength proxy $S$ | 103.4 MPa | Graphtek torsional strength, used as nominal longitudinal-shear proxy |

These are **effective preliminary properties**, not certified laminate design allowables. The model does not establish ply-level interaction, shell buckling, delamination, clamp crushing, impact, fatigue, environmental knockdown, or joint strength. Therefore, the reported material safety factor applies only to the modeled nominal tension/compression/shear checks; **overall structural safety remains unresolved**.

### 2.3 Design variables and bounds

The physical design vector is

$$
\mathbf d=[D,t]^T,
$$

where $D$ is uniform tube outer diameter and $t$ is uniform wall thickness. The inner diameter is

$$
D_i=D-2t.
$$

The provisional study bounds are

$$
10\ \text{mm}\le D\le60\ \text{mm},
$$

$$
2.01\ \text{mm}\le t\le3.00\ \text{mm},
$$

with an additional minimum-bore requirement

$$
D_i\ge1\ \text{mm}.
$$

The 2.01 mm lower wall bound implements the requirement that wall thickness be strictly greater than 2 mm. The 10–60 mm outer-diameter interval, 3 mm upper wall-thickness bound, and minimum-bore requirement are team-selected study bounds used to define a practical optimization domain rather than validated manufacturing or installation limits. A minimum nominal material safety factor of 2 is also adopted as a team-selected preliminary design requirement.

The hollow circular section properties are

$$
A=\frac{\pi}{4}(D^2-D_i^2),\qquad
I=\frac{\pi}{64}(D^4-D_i^4).
$$

### 2.4 Objective and explicit constraints

The optimization seeks the lightest admissible structural tube while enforcing material-strength, vibration, and geometric requirements. The nonlinear optimization evaluates the same strength requirements using the finite-rotation follower-load state. The objective is minimum tube mass:

$$
\boxed{\min_{D,t}\ m(D,t)=\rho A(D,t)L}.
$$

The constant tip propulsion mass is included in the structural/dynamic model but not in the tube-mass objective. Because tip mass is constant, minimizing total modeled mass gives the same tube geometry.

At 21 stations along the arm, the model evaluates outer-fiber axial/bending normal stress and nominal section shear. The safety-factor requirement is implemented using smooth KS aggregation with $\rho_{KS}=80$:

$$
g_t=KS_{80}\!\left(\frac{2\sigma_t}{X_t}\right)-1\le0,
$$

$$
g_c=KS_{80}\!\left(\frac{2|\sigma_c|}{X_c}\right)-1\le0,
$$

$$
g_s=KS_{80}\!\left(\frac{2\tau}{S}\right)-1\le0.
$$

The N2 optimization adds the same three strength constraints evaluated from the follower-load nonlinear state. The vibration requirement is

$$
g_f=\frac{f_{min}}{f_1}-1\le0,\qquad f_{min}=58.3333\ \text{Hz}.
$$

No deflection constraint is currently configured. Tip displacement is reported as an engineering response but is **not** an active requirement in this data set.

## 3. Analysis model

### 3.1 Linear Euler-Bernoulli finite-element model

The loaded bending plane is modeled with

$$
(EIw'')''=q,
$$

where $q=-\mu g$ is distributed tube self-weight. The fixed root satisfies

$$
w(0)=0,\qquad w'(0)=0.
$$

With net transverse tip force $P=T-m_tg$, the corresponding closed-form cantilever responses are

$$
w(L)=\frac{PL^3}{3EI}+\frac{qL^4}{8EI},
$$

$$
M(x)=P(L-x)+\frac{q}{2}(L-x)^2,
$$

$$
\sigma_b(x)=\frac{|M(x)|D}{2I}.
$$

`arm.py` implements an in-house sparse cubic-Hermite Euler-Bernoulli beam element with nodal degrees of freedom $[w,w']$. The production linear mesh has **20 elements**. Distributed weight and consistent translational mass are assembled element-wise. The root DOFs are removed from the solve, and support reactions are recovered from the unreduced system. The tip mass is added once to the translational mass matrix.

### 3.2 Loads and equilibrium

The case includes:

- combined upward coaxial thrust: **133.446648 N**;
- tip assembly weight: $m_tg$;
- distributed tube self-weight: $\mu g$;
- zero applied shaft torque.

Thrust alone produces a root moment of **81.3491 N m**. In the nonlinear model, gravity remains vertical while thrust rotates normal to the deformed tip tangent.

### 3.3 Stress recovery and composite screening

For the linear model, internal resultants are recovered from FE reactions and the exact uniform-load moment field. The outer-fiber normal stresses include axial and bending contributions,

$$
\sigma_{\pm}=\frac{N}{A}\pm\frac{MD}{2I}.
$$

A homogeneous-section shear screen is evaluated using

$$
\tau_{NA}=\frac{V(R^3-r_i^3)}{3It}.
$$

### 3.4 Vibration with concentrated tip mass

The independent FE eigenproblem is

$$
K\phi=\omega^2M\phi.
$$

For the production frequency used by the optimizer, `vibration.py` implements the exact uniform Euler-Bernoulli cantilever-with-tip-mass characteristic equation from Yntema, NACA TN 3459:

$$
H(\beta,r)=1+\cos\beta\cosh\beta-r\beta(\sin\beta\cosh\beta-\cos\beta\sinh\beta)=0,
$$

where

$$
r=\frac{m_t}{\mu L}.
$$

Each positive root gives

$$
f_n=\frac{\beta_n^2}{2\pi}\sqrt{\frac{EI}{\mu L^4}}.
$$

The FE matrix eigenproblem is retained as an independent discretization check. The model represents an unloaded stationary beam with a point tip mass and no tip rotary inertia. RPM defines the design threshold; the beam itself is not treated as rotating.

### 3.5 N2 finite-rotation follower-thrust model

The nonlinear component uses an inextensible planar elastica. Let $s$ be undeformed arc length and $\theta(s)$ the local centerline angle:

$$
x'=\cos\theta,\qquad z'=\sin\theta.
$$

The rotor thrust remains normal to the motor mount, so the load direction depends on the unknown tip rotation $\theta_L$. The governing equilibrium is

$$
(EI\theta')'+T\cos(\theta_L-\theta)+[-m_tg-\mu g(L-s)]\cos\theta=0,
$$

with

$$
\theta(0)=0,\qquad EI\theta'(L)=0.
$$

Linear angle elements with four-point Gauss integration form the nonlinear residual. Newton iteration uses load stepping, line search, and a consistent **nonsymmetric tangent**. Because follower force is nonconservative, using a symmetric energy Hessian would be incorrect. The nonlinear adjoint therefore solves the transposed tangent system

$$
R_\theta^Ta=Q_\theta^T,
$$

followed by

$$
\frac{dQ}{dd}=Q_d-a^TR_d.
$$

### 3.6 Modeling assumptions

The principal modeling assumptions are:

- the arm is straight, uniform, and perfectly clamped at the vehicle frame;
- Euler-Bernoulli beam theory is applicable, so transverse shear deformation is neglected in the primary bending model;
- effective homogeneous CFRP properties represent the tube in bending and strength screening;
- the motor assembly is modeled as a concentrated translational tip mass without rotary inertia;
- tube self-weight is uniformly distributed along the arm;
- the two coaxial rotors are represented by their combined thrust, with zero net applied shaft-reaction torque;
- the linear vibration model is evaluated about the unloaded stationary configuration;
- the N2 model is planar and inextensible, with thrust remaining normal to the deformed motor mount; and
- local shell effects, laminate interaction, joints, fatigue, damage, aerodynamic coupling, and control-system coupling are outside the model.

### 3.7 Independent 3D frame comparison

A separate linear 3D Euler-Bernoulli frame implementation (`drone_arm.py` through `compare_drone_arm.py`) provides an independent check that the reduced bending model reproduces the corresponding 3D frame response for this straight circular tube.

| Quantity | Reduced-model reference | Independent 3D frame | Difference |
|---|---:|---:|---:|
| Linear tip deflection (mm) | 1.240167843059 | 1.240167843070 | 1.163e-11 |
| First bending frequency (Hz) | 58.333333333333 | 58.333333363984 | 3.065e-08 |

The 3D frame contains **126 total DOFs (120 free)**. Its first two modes form the expected nearly degenerate two-plane bending pair at 58.3333 Hz. The very small deflection/frequency differences support the reduced bending model for the modeled planar load case, while not converting the main optimizer into a general 3D-frame optimization.

## 4. Analysis-code verification

The verification checks in this section use the prescribed base geometry, $D=40$ mm outer diameter and $t=2.5$ mm wall thickness. Here $D$ denotes outer diameter; the inner diameter is 35 mm. The full arm is a horizontal cantilever with combined upward thrust, a concentrated motor/rotor assembly mass and distributed tube self-weight. All calculations below are analysis or derivative checks to validate the finite element model. 

| Base-case input | Value |
|---|---:|
| Outer diameter $D$ | 40 mm |
| Wall thickness $t$ | 2.5 mm |
| Inner diameter $D-2t$ | 35 mm |
| Clamp-to-tip length $L$ | 0.6096 m |
| Effective bending modulus $E$ | 64 GPa |
| Effective density $\rho$ | 1520 $kg/m^3$ |
| Combined thrust $T$ | 133.446648457815 N |
| Concentrated tip mass $m_t$ | 0.680388555 kg |
| Gravity $g$ | 9.80665 $m/s^2$ |
| Net upward tip force $P=T-m_tg$ | 126.7743160349 N |
| Downward tube self-weight $q=-\rho Ag$ | -4.390211192 N/m |
| Linear FE mesh | 20 elements, 21 nodes, 40 free DOFs |
| Nonlinear verification mesh | 80 angle elements; 40/80/160 refinement also checked |

The isolated V1 and V2 tests deliberately change the active loads or supports to obtain classical closed-form benchmarks while retaining the base-case section, material and span. V3-V7 use the base-case properties and the configurations stated below. The nonlinear model is static; modal frequencies are those of an unloaded stationary beam with point tip mass.

### Verification ladder

| Test | Setup using the 40 mm OD x 2.5 mm wall base case | Exact or independent answer |
|---|---|---|
| **V1 Cantilever tip load** | $L=0.6096$ m, $E=64$ GPa, $P=126.7743160349$ N; suppress distributed self-weight to isolate the tip-load term | $\delta_L=PL^3/(3EI)$, $\theta_L=PL^2/(2EI)$, $\lvert M_0\rvert=PL$ |
| **V2 Simply supported beam** | Same section and length; apply only $q=-\rho Ag=-4.390211192$ N/m; pin both ends and leave rotations free | $\delta_{mid}=5qL^4/(384EI)$ |
| **V3 Mesh convergence** | Same cantilever with net tip load and tube self-weight; use 2, 4, 8, and 16 cubic-Hermite elements | Fourth-order convergence of the displacement-field $L_2$ error |
| **V4 Beam frequency** | Same cantilever, $m_t=0.680388555$ kg, $\rho=1520$ kg/m$^3$; compare the 20-element mass/stiffness eigenproblem with the characteristic equation | Cantilever-with-point-tip-mass characteristic solution |
| **V5 Global checks** | Full base-case cantilever loading for equilibrium; release the supports for a separate planar rigid-mode count | Zero force/moment imbalance and two free-free planar rigid modes |
| **V6 Design gradients** | Evaluate all 18 mass/constraint derivatives at $(D,t)=(40,2.5)$ mm; use 20 linear and 40 nonlinear elements | Central finite differences; first-order Taylor remainder proportional to $h^2$ |
| **V7 N2 follower-load benchmark** | Same physical arm, gravity and follower thrust; compare 80 nonlinear angle elements with an independent continuum BVP | Continuum follower-thrust tip displacement under the same loads |

The principal numerical comparisons for this base case are:

| Test | Numerical result | Exact / independent result | Error | Status |
|---|---:|---:|---:|---|
| **V1** Isolated tip-load cantilever displacement | 0.0028763736649 m | 0.0028763736649 m | $1.1795\times10^{-11}$ | PASS |
| **V2** Simply supported self-weight midspan displacement | -2.3719440104e-06 m | -2.3719440105e-06 m | $8.0025\times10^{-12}$ | PASS |
| **V3** Displacement convergence order | 4.000323851 | 4.000000000 | $8.0963\times10^{-5}$ | PASS |
| **V4** First frequency with tip mass | 38.7142512070 Hz | 38.7142507467 Hz | $1.1889\times10^{-8}$ | PASS |
| **V5** Global force balance | $1.9408\times10^{-9}$ N residual | 0 N | $1.9408\times10^{-9}$ N (absolute) | PASS |
| **V6** Maximum normalized derivative discrepancy | $6.3342\times10^{-8}$ | Central differences, 18 entries | $6.3342\times10^{-8}$ | PASS |
| **V7** N2 follower-load tip displacement, 80 elements | 0.0028534721428 m | 0.0028535833141 m | $3.8959\times10^{-5}$ | PASS |

The error column is a dimensionless relative error except for V5, where the exact target is zero and an absolute force residual is reported. For V3, the error is $|p-4|/4$. V6 reports the maximum discrepancy across all 18 checked entries, using the derivative checker's normalization. PASS means that the numerical verification tolerance is met.

The displacement checks use relative tolerance $10^{-9}$; V3 permits a 1% relative difference from order 4; V4 uses relative tolerance $2\times10^{-6}$; V5 requires force and moment imbalance below $10^{-7}$ in N and N m, respectively, and exactly two planar rigid modes. V6 uses $3\times10^{-4}$ for constraint derivative discrepancies and $10^{-6}$ for mass derivatives. V7 uses relative tolerance $10^{-4}$. All seven rungs pass these checks.

### 4.1 Classical beam verification and MATLAB comparison

For the annular section, $A=\pi[D^2-(D-2t)^2]/4=2.9452431127\times10^{-4}$ $m^2$ and $I=\pi[D^4-(D-2t)^4]/64=5.2001948709\times10^{-8}$ $m^4$. The MATLAB hand calculations and Python model use the same $E$, density, length, thrust, tip mass, gravity and exact pound conversion. The tube mass is $m_b=\rho AL=0.2729038706$ kg and its total weight is $W_b=m_bg=2.6762727429$ N.

V1 isolates the net tip force by suppressing the tube's distributed gravity load. V2 isolates the tube self-weight and changes the supports to a simply supported beam. Their formulas are listed in the verification ladder. Thus V1's tip-load-only displacement is larger than the displacement of the actual arm, where downward tube self-weight reduces the net upward deflection.

For the actual cantilever under both loads, the MATLAB hand calculations use superposition:

$$
\delta_L=\frac{PL^3}{3EI}-\frac{W_bL^3}{8EI},
\qquad
\theta_L=\frac{PL^2}{2EI}-\frac{W_bL^2}{6EI},
$$

$$
|M_0|=PL-\frac{W_bL}{2},
\qquad
\sigma_{\max}=\frac{|M_0|D}{2I}.
$$

The actual-arm Python solution gives $\delta_L=2.8536030024$ mm, $\theta_L=0.0070278865$ rad, $|M_0|=76.4658951219$ N m and $\sigma_{\max}=29.4088575601$ MPa. The signed root reaction moment is negative under the chosen convention; the expression above reports its magnitude. These actual-arm results are separate from the isolated V1 and V2 checks.

Errors are computed from the unrounded FE result and the supplied MATLAB number:

```math
\varepsilon_{\%}
=100\frac{\lvert r_{\mathrm{FE}}-r_{\mathrm{MATLAB}}\rvert}
{\lvert r_{\mathrm{MATLAB}}\rvert}
```

$$
\begin{array}{lrrr}
\hline
\text{Quantity} & \text{Python FEA} & \text{MATLAB output} & \text{Percentage Difference}\\
\hline
\text{Tip deflection (mm)} & 2.8536030024 & 2.8536030024 & 1.6743\times10^{-9} \\
\text{Peak bending stress (MPa)} & 29.4088575601 & 29.4088575604 & 1.1443\times10^{-9} \\
f_{1}\;\text{(Hz)} & 38.7142512070 & 38.7142507467 & 1.1889\times10^{-6} \\
f_{2}\;\text{(Hz)} & 582.6353810913 & 582.6347793544 & 1.0328\times10^{-4} \\
f_{3}\;\text{(Hz)} & 1859.4944976727 & 1859.4744067567 & 1.0805\times10^{-3} \\
f_{4}\;\text{(Hz)} & 3864.2776029757 & 3864.0965664807 & 4.6851\times10^{-3} \\
f_{5}\;\text{(Hz)} & 6598.5039727188 & 6597.6046499335 & 1.3631\times10^{-2} \\
\hline
\end{array}
$$

The static results agree to about $10^{-9}\%$. The largest frequency discrepancy is $0.013631\%$ in mode 5. This validates that our finite element model concurs with classical theory and can be used in the optimization for modeling the structure. 

### 4.2 Mesh convergence

![Base-case mesh convergence](assets/section4-basecase-mesh-convergence.png)

The same base-case cantilever is solved with 2, 4, 8 and 16 cubic-Hermite elements. Displacement-field error is integrated against the exact superposed tip-load and uniform-load solution using five-point Gauss quadrature on each element. The error is the unnormalized $L_2$ norm, $\sqrt{\int_0^L(w_h-w)^2\,dx}$, with units $m^{3/2}$. The point-load contribution is cubic and exactly representable; the quartic self-weight term drives the interpolation error.

| Linear elements | Displacement-field $L_2$ error ($m^{3/2}$) | FE first frequency (Hz) | Relative frequency error |
|---:|---:|---:|---:|
| 2 | $1.4757\times10^{-8}$ | 38.7143930749 | $3.6764\times10^{-6}$ |
| 4 | $9.2229\times10^{-10}$ | 38.7142604852 | $2.5155\times10^{-7}$ |
| 8 | $5.7643\times10^{-11}$ | 38.7142513598 | $1.5837\times10^{-8}$ |
| 16 | $3.6000\times10^{-12}$ | 38.7142508007 | $1.3944\times10^{-9}$ |

The fitted displacement convergence order is **4.000323851**, close to the theoretical value of 4. The first-frequency convergence order, measured from the 2- and 8-element errors to avoid the fine-mesh eigenvalue roundoff floor, is **3.929437419**. These are base-case measurements.

### 4.3 Natural-frequency verification

The uniform cantilever carries a concentrated tip mass but no tip rotary inertia. The independent characteristic equation used by the MATLAB hand calculation is

$$
1+\cos\beta\cosh\beta+
\eta\beta(\cos\beta\sinh\beta-\sin\beta\cosh\beta)=0,
\qquad \eta=\frac{m_t}{\rho AL}=2.4931436605.
$$

The positive roots give

$$
f_n=\frac{\beta_n^2}{2\pi L^2}\sqrt{\frac{EI}{\rho A}}.
$$

For the base case, the 20-element FE model gives $f_1=38.7142512070$ Hz, compared with the supplied MATLAB result $38.7142507467$ Hz. The relative difference is approximately $1.19\times10^{-8}$, or $1.19\times10^{-6}\%$. The five-mode comparison is shown in Section 4.1; the maximum frequency discrepancy is $0.013631\%$ in mode 5. Both calculations omit load-dependent tangent stiffness, gyroscopic terms and operating follower-force dynamic stability.

### 4.4 Global equilibrium and rigid-body behavior

Releasing the supports of the same reduced planar beam produces exactly **two rigid-body modes**, corresponding to transverse translation and in-plane rotation. The clamped full-load arm is checked independently for force and moment balance using assembled nodal forces and recovered root reactions.

$$
\left|R_z+\sum_i F_{z,i}\right|=1.9408\times10^{-9}\;\mathrm{N},
\qquad
\left|M_R+\sum_i x_iF_{z,i}+\sum_i M_i\right|=9.3303\times10^{-10}\;\mathrm{N\,m}.
$$

Both residuals are below the $10^{-7}$ absolute verification tolerance.

### 4.5 Gradient verification

All 18 objective/constraint derivative entries are checked at the base point $(D,t)=(40,2.5)$ mm using the supplied design-dependent material/mass model, 20 linear beam elements and 40 nonlinear angle elements. This evaluates the derivative code without executing an optimizer. Constraint derivatives have units $mm^{-1}$ and mass derivatives kg/mm. Central differences use $\Delta d_j=10^{-4}\max(|d_j|,1)$ in mm. The normalized discrepancy divides by $\max(|g_{FD}|,|g_{adj}|,10^{-7})$ for constraint derivatives; the mass check uses its mass-derivative normalization.

| Response | Variable | Analytic / adjoint derivative | Central difference | Normalized relative error |
|---|---|---:|---:|---:|
| linear_tension | OD_mm | -0.00440790410738 | -0.00440790424314 | $3.0798\times10^{-8}$ |
| linear_compression | OD_mm | -0.00295088435809 | -0.00295088445813 | $3.3902\times10^{-8}$ |
| linear_shear | OD_mm | -0.00044014558882 | -0.000440145573599 | $3.4581\times10^{-8}$ |
| nonlinear_tension | OD_mm | -0.00440689171423 | -0.00440689183294 | $2.6938\times10^{-8}$ |
| nonlinear_compression | OD_mm | -0.0029514929869 | -0.00295149307125 | $2.8579\times10^{-8}$ |
| nonlinear_shear | OD_mm | -0.000440130793445 | -0.000440130798349 | $1.1142\times10^{-8}$ |
| bending_frequency | OD_mm | -0.0583495554261 | -0.0583495563762 | $1.6283\times10^{-8}$ |
| inner_diameter | OD_mm | -0.1 | -0.0999999999999 | $6.6530\times10^{-13}$ |
| mass | OD_mm | 0.00727743655019 | 0.00727743655018 | $1.3405\times10^{-12}$ |
| linear_tension | wall_mm | -0.0262272457859 | -0.0262272460227 | $9.0293\times10^{-9}$ |
| linear_compression | wall_mm | -0.0175539197569 | -0.0175539198606 | $5.9086\times10^{-9}$ |
| linear_shear | wall_mm | -0.006239459337 | -0.00623945973222 | $6.3342\times10^{-8}$ |
| nonlinear_tension | wall_mm | -0.0262207744575 | -0.0262207747501 | $1.1161\times10^{-8}$ |
| nonlinear_compression | wall_mm | -0.0175579170121 | -0.017557917217 | $1.1670\times10^{-8}$ |
| nonlinear_shear | wall_mm | -0.00623937743334 | -0.00623937749555 | $9.9717\times10^{-9}$ |
| bending_frequency | wall_mm | -0.219520891001 | -0.219520892593 | $7.2535\times10^{-9}$ |
| inner_diameter | wall_mm | 0.2 | 0.200000000001 | $6.5511\times10^{-12}$ |
| mass | wall_mm | 0.101884111703 | 0.101884111703 | $1.5801\times10^{-13}$ |

The largest normalized discrepancy is **6.334179e-08**. The following directional Taylor test uses $d(h)=d_0+h\,\mathrm{diag}(40,2.5)(0.6,0.8)^T$ and the combined mass/constraint output vector. Each remainder component is normalized by $\max(|f_i(d_0)|,1)$ before taking the vector norm.

| Normalized perturbation $h$ | Normalized first-order Taylor remainder |
|---:|---:|
| $1.0\times10^{-2}$ | $1.8130\times10^{-4}$ |
| $3.0\times10^{-3}$ | $1.6417\times10^{-5}$ |
| $1.0\times10^{-3}$ | $1.8273\times10^{-6}$ |
| $3.0\times10^{-4}$ | $1.6456\times10^{-7}$ |

The measured remainder order is **1.997661810**, consistent with the expected quadratic behavior. These tables are the checks at the 40 mm base point; any separately reported study-wide derivative checks in Section 5 retain their own evaluation points.

### 4.6 N2 nonlinear benchmark

The same physical base case is evaluated with finite rotations and follower thrust normal to the deformed tip tangent. Tip mass gravity and tube self-weight retain their global vertical directions. An independent continuum boundary-value problem integrates the angle/curvature equations and deformed coordinates, with fixed root angle, zero free-end bending moment, and the follower direction tied to the unknown tip angle. It does not reuse the nonlinear FE residual.

| Nonlinear angle elements | FE tip displacement (mm) | Continuum tip displacement (mm) | Relative error |
|---:|---:|---:|---:|
| 40 | 2.8531386287 | 2.8535833141 | $1.5583\times10^{-4}$ |
| 80 | 2.8534721428 | 2.8535833141 | $3.8959\times10^{-5}$ |
| 160 | 2.8535555213 | 2.8535833141 | $9.7396\times10^{-6}$ |

The 80-element nonlinear benchmark gives **0.0028534721428 m** tip displacement versus **0.0028535833141 m** from the continuum BVP. Its relative error is **3.895851e-05**, satisfying the $10^{-4}$ criterion. Refinement from 40 to 80 to 160 elements reduces this error by approximately a factor of four each time, consistent with the linear angle interpolation.

A complex-step check of the 80-element follower-load tangent gives relative discrepancy **1.991757e-16**. Newton residuals for that solve fall from 9.657654e-01 to 2.388773e-05 to 1.080025e-12. This small-rotation base case reaches numerical precision rapidly; a separate large-rotation fixture remains available in `verify.py` for broader nonlinear verification.

The continuum nonlinear displacement differs from the linear result by **-0.000690%**. All seven base-case verification rungs pass, while the base-case first frequency remains below the study's $3500/60=58.333333$ Hz sizing requirement.

## 5. Sensitivity analysis

### 5.1 Linear reduced-space adjoint

For the linear residual

$$
h(d,u)=K(d)u-F(d)=0,
$$

and scalar output $Q(d,u)$, the adjoint is defined by

$$
K^Ta=Q_u^T.
$$

The total derivative is

$$
\boxed{\frac{dQ}{dd_j}=Q_{d_j}-a^T\left(K_{d_j}u-F_{d_j}\right)}.
$$

The tip-displacement response uses the tip translation as $Q$. Root stress includes both state dependence and explicit section dependence. Complex-step derivatives evaluate fixed-state partials; the adjoint handles the state sensitivity.

### 5.2 Eigenvalue/frequency sensitivities

The exact tip-mass characteristic root is differentiated implicitly:

$$
\frac{\partial\beta}{\partial r}=-\frac{H_r}{H_\beta}.
$$

Thus

$$
d\ln f=\tfrac12d\ln(EI)-\tfrac12d\ln\mu-2d\ln L+
\frac{2}{\beta}\frac{\partial\beta}{\partial r}dr.
$$

For the independent FE eigenproblem, a mass-normalized mode gives

$$
\frac{d\omega^2}{dd_j}=\phi^T(K_{d_j}-\omega^2M_{d_j})\phi,
$$

$$
\frac{df}{dd_j}=\frac{1}{8\pi^2f}\frac{d\omega^2}{dd_j}.
$$

### 5.3 N2 follower-load adjoint

For the converged nonlinear residual $R(\theta,d)=0$ and scalar response $Q(\theta,d)$, the implicit-differentiation adjoint satisfies

$$
R_\theta^T a = Q_\theta^T,
$$

and the total design derivative is

$$
\frac{dQ}{dd}=Q_d-a^T R_d.
$$

Because the follower load is nonconservative, $R_\theta$ is nonsymmetric and the transposed tangent is required in the adjoint solve.


![Nonlinear follower-force gradient check](assets/nonlinear-gradient.png)

The nonlinear sensitivity check compares the transposed-tangent adjoint against independently converged finite-difference solves as the EI perturbation is reduced. The V9 follower-force adjoint has relative error **8.25e-09** at its accepted check point, and the tangent itself matches complex-step differentiation to approximately **2.0e-16** relative error.

### 5.4 Gradient checks for every response and design variable
The following derivatives are evaluated at D = 35.000 mm and
t = 2.505 mm using 20 linear elements and 40 nonlinear elements.
They are separate from the D = 40 mm, t = 2.5 mm checks in Section 4.5.

![Directional gradient checks](assets/gradient-check.png)

The plot compares adjoint/eigenvalue directional derivatives against central differences for tip displacement, root stress, and bending frequency in the two coordinate directions plus three seeded random directions. At coarse steps, the error decreases approximately quadratically; at very small steps, subtraction/roundoff begins to dominate. Every recorded directional check passes the specified tolerance.



All **18** mass/constraint derivatives used by the OD/wall optimizer pass central-difference checks, including density-dependent self-weight, linear and nonlinear strength constraints, frequency, bore constraint, and mass:

| Response | Variable | Analytic/adjoint | Central difference | Relative error |
|---|---|---:|---:|---:|
| linear_tension | OD_mm | -0.007111329566 | -0.00711132946 | 1.48e-08 |
| linear_compression | OD_mm | -0.004853819732 | -0.004853819661 | 1.46e-08 |
| linear_shear | OD_mm | -0.0005832705646 | -0.0005832705453 | 3.32e-08 |
| nonlinear_tension | OD_mm | -0.007109079903 | -0.007109080093 | 2.67e-08 |
| nonlinear_compression | OD_mm | -0.004855049785 | -0.00485504992 | 2.77e-08 |
| nonlinear_shear | OD_mm | -0.0005832243874 | -0.0005832243939 | 1.13e-08 |
| bending_frequency | OD_mm | -0.08307331552 | -0.0830733169 | 1.66e-08 |
| inner_diameter | OD_mm | -0.1 | -0.1 | 8.41e-13 |
| mass | OD_mm | 0.007291991423 | 0.007291991423 | 2.34e-12 |
| linear_tension | wall_mm | -0.0352071494 | -0.03520714179 | 2.16e-07 |
| linear_compression | wall_mm | -0.02402615264 | -0.0240261472 | 2.27e-07 |
| linear_shear | wall_mm | -0.007101747556 | -0.007101747743 | 2.63e-08 |
| nonlinear_tension | wall_mm | -0.03519515103 | -0.03519515142 | 1.1e-08 |
| nonlinear_compression | wall_mm | -0.02403312197 | -0.02403312224 | 1.13e-08 |
| nonlinear_shear | wall_mm | -0.007101531948 | -0.007101532019 | 1e-08 |
| bending_frequency | wall_mm | -0.2631748937 | -0.2631748957 | 7.36e-09 |
| inner_diameter | wall_mm | 0.2 | 0.2 | 1.16e-12 |
| mass | wall_mm | 0.08730012886 | 0.08730012886 | 2.05e-13 |

### 5.5 First-order Taylor remainder

A separate first-order Taylor test uses the same verified synthetic two-variable beam fixture. For a direction $v$,

$$
R(h)=|Q(d+hv)-Q(d)-h\nabla Q(d)^Tv|.
$$

For a correct first derivative, $R(h)=O(h^2)$. It shows the normalized tip-displacement Taylor remainder for an outer-diameter perturbation direction. Both axes are plotted as base-10 logarithms; therefore a slope of 2 is the expected first-order Taylor-test behavior.

```mermaid
xychart-beta
    title "First-order Taylor remainder: tip displacement"
    x-axis "log10(h)" [-5.000, -4.523, -4.000, -3.523, -3.000, -2.523, -2.000, -1.523, -1.000]
    y-axis "log10(normalized remainder)" -10 --> 0
    line [-9.179, -8.225, -7.179, -6.225, -5.180, -4.227, -3.187, -2.247, -1.250]
```

The fitted log-log slopes from the same calculation are **1.988 for tip displacement, 1.990 for root stress, and 1.999 for bending frequency**. These values are essentially the expected slope of 2 and independently confirm the first-order derivative implementation. The underlying normalized remainder values for the plotted tip-displacement series are:

| $h$ | Normalized Taylor remainder |
|---:|---:|
| $1.0\times10^{-5}$ | $6.625\times10^{-10}$ |
| $3.0\times10^{-5}$ | $5.958\times10^{-9}$ |
| $1.0\times10^{-4}$ | $6.620\times10^{-8}$ |
| $3.0\times10^{-4}$ | $5.955\times10^{-7}$ |
| $1.0\times10^{-3}$ | $6.609\times10^{-6}$ |
| $3.0\times10^{-3}$ | $5.927\times10^{-5}$ |
| $1.0\times10^{-2}$ | $6.507\times10^{-4}$ |
| $3.0\times10^{-2}$ | $5.660\times10^{-3}$ |
| $1.0\times10^{-1}$ | $5.621\times10^{-2}$ |

### 5.6 Cost comparison with finite differences

The computational cost was measured on the production drone-arm model at the optimizer gradient-check point, $D=35.000$ mm and $t=2.505$ mm. The production discretization uses **20 linear beam elements** and **40 nonlinear follower-load elements**. The analytic/adjoint evaluation computes the complete objective and constraint Jacobian for all eight enabled constraints, while the central-difference calculation perturbs both physical design variables and reruns the production analyses.

| Method | Median measured time |
|---|---:|
| Full analytic/adjoint Jacobian | 119.980 ms |
| Central-difference Jacobian, two variables | 319.529 ms |
| Measured speedup | 2.663x |

The two methods also produce essentially identical derivatives. The maximum relative discrepancy among the constraint gradients is $7.09\times10^{-8}$, while the maximum relative discrepancy in the mass gradient is $2.34\times10^{-12}$.

Because the current uniform-tube optimization has only two physical design variables, the measured computational advantage of the adjoint method is modest but still significant: the complete analytic/adjoint Jacobian is approximately **2.66 times faster** than central differences on the production model. Central differences require separate perturbed production analyses for each design variable, whereas the adjoint formulation scales primarily with the number of scalar responses rather than the number of design variables. The computational advantage would therefore become more significant for higher-dimensional sizing or shape parameterizations.
## 6. Optimization

### 6.1 Algorithm, scaling, and starting designs

The design is optimized with **SciPy SLSQP**. Objective mass is normalized by the reference mass, and outer diameter and wall thickness are scaled by their reference values so that the optimizer works with dimensionless quantities of comparable magnitude. A candidate is accepted only if the solver reports success and the final point passes feasibility, stationarity, dual-feasibility, and complementarity checks.

The same problem was solved from six deliberately different starting geometries spanning the design space. The table below shows the starts for the N2 optimization; the linear-only optimization was run from the same six geometries and converged to the same final design.

| Start | Initial OD (mm) | Initial wall (mm) | Initial tube mass (g) | Iterations | Final tube mass (g) | Result |
|---|---:|---:|---:|---:|---:|---|
| Reference design | 40.0 | 2.50 | 272.904 | 6 | 311.714 | PASS |
| Configured start 2 | 30.0 | 2.20 | 178.035 | 7 | 311.714 | PASS |
| Configured start 3 | 18.0 | 2.80 | 123.891 | 8 | 311.714 | PASS |
| Small design-box corner | 10.0 | 2.01 | 46.750 | 10 | 311.714 | PASS |
| Large design-box corner | 60.0 | 3.00 | 497.777 | 6 | 311.714 | PASS |
| Feasible fixed-wall reference | 52.3 | 2.50 | 362.416 | 5 | 311.714 | PASS |

All six starts converge to the same accepted optimum, providing evidence that the reported solution is robust within the two-variable study region.

### 6.2 Convergence history

![Optimization convergence](assets/optimization-convergence.png)

The left panel shows tube mass versus SLSQP iteration for all six starts. Although the starting masses span a wide range, every path converges to approximately **311.714 g**. The right panel shows positive constraint violation falling to numerical tolerance. The design is therefore both converged and feasible; it is not accepted solely because SLSQP returned a success flag.

### 6.3 Optimized design and feasible reference

The optimized N2 design is

$$
\boxed{D=55.2848\ \text{mm},\qquad t=2.0100\ \text{mm}},
$$

with

$$
\boxed{m_{tube}=311.714\ \text{g}}.
$$

The baseline 40mm x 2.5mm tube is infeasible because $f_1=38.714$ Hz is below the required 58.333 Hz, so comparing mass directly to that start is not meaningful as an optimization benefit. The lightest **feasible fixed-wall** reference found on the 0.1 mm OD grid is

$$
52.3\ \text{mm OD}\times2.5\ \text{mm wall},\qquad m=362.416\ \text{g}.
$$

Allowing both diameter and wall thickness to vary reduces tube mass by 13.990% relative to that feasible fixed-wall design.

A rounded custom-tube candidate of **55.3 mm OD × 2.1 mm wall** has mass **325.214 g**, minimum modeled SF **38.31**, linear tip deflection **1.191 mm**, and first bending frequency **59.378 Hz**.

### 6.4 Active constraints, KKT multipliers, and engineering interpretation

| Active requirement | Scaled multiplier | Residual |
|---|---:|---:|
| bending_frequency | 0.78795151 | -9.99e-16 |
| wall_lower_bound | 0.39090122 | 0 |

The KKT stationarity infinity norm is **4.44e-16**, complementarity infinity norm is **7.87e-16**, and maximum violation is **0**. The active set is full-rank for the two-variable problem.

The optimum is controlled by two facts:

1. **Bending frequency is active.** Any meaningful reduction in bending stiffness/mass distribution would push the first mode below 58.333 Hz.
2. **Wall thickness is at its 2.01 mm lower bound.** Once the minimum wall is reached, the optimizer adjusts diameter until the frequency requirement is just met.

Material strength is therefore not controlling the current optimum. The design is governed by the interaction between the minimum permitted wall thickness and the required bending frequency.

The shadow-price checks validate the multiplier interpretation by tightening each active requirement and reoptimizing:

| Tightened requirement | Step | Predicted mass change (g) | Reoptimized mass change (g) |
|---|---:|---:|---:|
| wall_lower_bound | 0.001 mm | 0.106678 | 0.106671 |
| bending_frequency | 0.01 Hz | 0.036863 | 0.036862 |

The predicted and reoptimized mass changes agree to within the shown numerical precision. For example, increasing the wall lower bound by 0.1 mm corresponds locally to about **10.668 g** additional tube mass while the active set remains unchanged.

## 7. N2 nonlinear component

### 7.1 Physical nonlinear effect

The nonlinear load follows the motor-mount rotation, so the force direction depends on the unknown deformed state. This makes the residual genuinely nonlinear and its tangent nonsymmetric. Newton's method, load stepping, and a transposed-tangent adjoint are therefore required. The verification ladder independently checks the nonlinear solution, tangent, Newton convergence, and implicit gradient.

### 7.2 Nonlinear load path

![Nonlinear load path](assets/nonlinear-path.png)

The exploratory path scales thrust and both gravity contributions together for two retained reference geometries used only to illustrate nonlinear load-response behavior. Solid curves are the finite-rotation elastica response; dashed curves are the corresponding linear predictions. At the design load multiplier of 1, the nonlinear correction is very small. At much larger artificial load factors, curvature increases and the nonlinear response separates more visibly from the linear prediction. This plot is a numerical nonlinearity demonstration, **not** an RPM/throttle map or a proof of strength at those amplified loads.

### 7.3 Design with and without N2

| Quantity | Linear-only optimization | Linear + N2 follower constraints | Change |
|---|---:|---:|---:|
| OD (mm) | 55.284816 | 55.284816 | 0 |
| Wall (mm) | 2.010000 | 2.010000 | 0 |
| Tube mass (g) | 311.714092 | 311.714092 | 0 |
| First bending frequency (Hz) | 58.333333 | 58.333333 | 0 |
| Linear tip deflection (mm) | 1.240168 | 1.240168 | 0 |
| N2 tip deflection (mm) | Not evaluated | 1.239973 | - |
| Minimum modeled SF | 36.810241 | 36.810241 | 0 |

Adding the follower-load strength constraints changes optimized mass by
0%. At the optimized geometry, the linear tip displacement is
1.240167843 mm, while the independent nonlinear continuum reference gives
1.240166229 mm, a reduction of approximately 0.000130%.

The production 40-element nonlinear model gives 1.239973043 mm. Most of
its difference from the linear result is discretization error rather
than physical geometric nonlinearity. Refining to 160 nonlinear elements
gives 1.240154155 mm and approaches the continuum reference. The nonlinear
constraints remain inactive, so including them does not change the optimum.
### 7.4 Optimized-design mesh refinement

| Linear elements | Nonlinear elements | Min modeled SF | Nonlinear tip (mm) | First frequency (Hz) | Max constraint residual |
|---:|---:|---:|---:|---:|---:|
| 20 | 40 | 36.81024116 | 1.239973 | 58.333333 | -9.99e-16 |
| 40 | 80 | 36.81024116 | 1.240118 | 58.333333 | -9.99e-16 |
| 80 | 160 | 36.81024118 | 1.240154 | 58.333333 | -9.99e-16 |

The optimized responses remain essentially unchanged as the linear mesh is refined from 20 to 80 elements and the nonlinear mesh from 40 to 160 elements. This supports the adequacy of the production discretization for the reported quantities.

## 8. Initial-versus-optimized results

The following table reports the objective and **every enabled constraint family** for the starting and final designs. Deflection is not included as a constraint because no deflection limit is configured.

| Metric | Starting design | Optimized design | Interpretation |
|---|---:|---:|---|
| Tube mass | 272.904 g | 311.714 g | Start is lighter but infeasible |
| Total modeled mass incl. tip assembly | 953.292 g | 992.103 g | Constant tip mass included |
| Linear tip deflection | 2.8536 mm | 1.2402 mm | Reported response, no limit configured |
| N2 tip deflection | 2.8531 mm | 1.2400 mm | Follower-load response |
| First bending frequency | 38.7143 Hz | 58.3333 Hz | Requirement: >= 58.3333 Hz |
| Minimum modeled nominal SF | 22.1022 | 36.8102 | Requirement: >= 2 |

| Constraint | Starting design | Optimized | Requirement | Start | Optimized |
|---|---:|---:|---|---|---|
| Linear tension safety factor | 22.102185 | 36.810241 | >= 2 | PASS | PASS |
| Linear compression safety factor | 30.946459 | 51.540001 | >= 2 | PASS | PASS |
| Linear shear safety factor | 120.465702 | 137.321615 | >= 2 | PASS | PASS |
| Nonlinear tension safety factor | 22.104677 | 36.812798 | >= 2 | PASS | PASS |
| Nonlinear compression safety factor | 30.943247 | 51.536508 | >= 2 | PASS | PASS |
| Nonlinear shear safety factor | 120.465545 | 137.321582 | >= 2 | PASS | PASS |
| First bending frequency (Hz) | 38.714251 | 58.333333 | >= 58.333333 | FAIL | PASS |
| Outer diameter (mm) | 40.000000 | 55.284816 | 10 <= D <= 60 | PASS | PASS |
| Wall thickness (mm) | 2.500000 | 2.010000 | 2.01 <= t <= 3.00 | PASS | PASS (active lower bound) |
| Inner diameter (mm) | 35.000000 | 51.264816 | >= 1 | PASS | PASS |

The frequency constraint is the only enabled performance requirement failed by the starting design. The optimized design satisfies all enabled nominal constraints, but this should not be interpreted as full structural certification because shell, joint, laminate-interaction, fatigue, and dynamic-stability failure modes are outside the model.

## 9. Fully nonlinear N2 extension

The follower-load analysis uses a genuinely nonlinear state equation rather than a linearized post-processing correction. The equilibrium equations are solved with Newton iteration, load stepping, and line search; the converged nonlinear state is differentiated implicitly with the transposed tangent; and the implementation is checked against an independently formulated nonlinear continuum boundary-value problem. The resulting optimization is therefore supported by a verified nonlinear solve and verified nonlinear sensitivities. At the design load, the N2 response is close to the linear response and does not change the optimum, which is itself an engineering result rather than a failure of the nonlinear model.

## 10. Reproducibility

Run the following from `Project 4 Luc/` with the dependencies in `requirements.txt` installed:

```bash
python optimize_arm.py --update-report
python -m unittest test_optimize test_vibration -v
```

The optimizer reads the loads from `inputs.json` and the design study from `optimizer_inputs.json`. It checks the objective and constraint derivatives, solves both the linear-only and N2 problems from six starting designs, and checks the optimum through mesh refinement and shadow-price reoptimization. It saves the numerical results and convergence history in `output/optimization/`. A rerun reproduced the reported optimum of **55.284816 mm outer diameter, 2.010 mm wall thickness, 311.714 g tube mass, and 58.333333 Hz first bending frequency**. All **18 unit tests passed**. The `--update-report` option also invokes `run.py` to regenerate the verification outputs, figures, and `assets/results.json`. The published page loads the written report and tables from `REPORT.md` and the 3D viewer data from `docs/assets/results.json`. After rebuilding, copy `assets/results.json` to `docs/assets/results.json`.
If the results changed, update the affected values in `REPORT.md` and copy
any regenerated figures it references to `docs/assets/`. Commit those
updates to refresh the published site.

## 11. Limitations and AI-use statement

### 11.1 Engineering limitations

This is a preliminary structural sizing study, not a flight-certified composite design. The main limitations are:

- effective beam properties are used instead of a ply-level laminate/shell model;
- local shell buckling, Brazier ovalization, delamination, clamp crushing, joint/fastener behavior, fatigue, impact, and environmental effects are not modeled;
- the shear allowable is a transferred comparison value and not a validated interlaminar allowable;
- the vibration model considers unloaded bending modes with a point tip mass but omits rotor rotary inertia, gyroscopic coupling, blade-passing harmonics, controller interaction, and flutter;
- the frequency threshold is only $1\times$ shaft frequency and has no additional separation margin;
- no deflection/serviceability limit is currently configured;
- the main optimizer is a reduced planar bending/N2 model rather than a general 3D-frame optimization.

### 11.2 Reduced design parameterization

The drone arm is a structure listed in the Project 4 specification. This implementation uses a deliberately reduced uniform-tube parameterization with two physical sizing variables, outer diameter and wall thickness. This makes the optimization problem inexpensive and transparent, but it also means that the computational speed advantage of adjoint sensitivities is smaller than it would be for a many-variable sizing or shape problem. The report therefore presents the measured adjoint-versus-finite-difference cost directly rather than claiming a large scaling benefit for this two-variable implementation.

### 11.3 AI-use statement

A coding agent researched the linked sources and drafted the FE solver, finite-rotation follower-force solver, vibration implementation, sensitivity routines, optimizer, verification checks, viewer, and report-generation workflow. The numerical work was not accepted solely because it executed: it is checked against closed-form beam solutions, independent eigenvalue formulas, mesh convergence, global equilibrium, independent nonlinear boundary-value solutions, complex-step tangent checks, central finite differences, Taylor remainder behavior, multistart convergence, KKT conditions, shadow-price reoptimization, and a separate 3D frame comparison. Additionally, classical theory implemented by hand confirmed the original base case and all optimized results concurred with the finite element model. This script was written independent of any coding agent and can be found 'hand_calcs_prj4.m' in the repository.  

One concrete bug caught during verification was a floating-point-only nonlinear tangent accumulator: complex-step perturbations caused a casting failure. The accumulator was changed to preserve the perturbation dtype, after which the nonlinear tangent and adjoint checks passed. The load model was also revised after the engineering brief clarified that thrust follows the motor mount rather than remaining a dead load.

## References

1. Course Project 4 specification, *PDE-Constrained Design Optimization - An Introduction to Optimization for Mechanical Engineers*.
2. Easy Composites, *CFT-WF Carbon Fibre Woven Finish Tube* technical data sheet: https://media.easycomposites.co.uk/datasheets/EC-TDS-Carbon-Fibre-Woven-Finish-Tube.pdf
3. Graphtek / GraphiteStore, roll-wrapped CFRP comparison property sheet: https://www.graphitestore.com/core/media/media.nl?c=4343521&h=22DqRe5XcXtOFJUqMtPHOuVeC4o8PE9P-gnJuKBYWPrJ3COD&id=9929
4. Yntema, NACA TN 3459, uniform cantilever with concentrated tip mass, equations (13)-(14): https://ntrs.nasa.gov/api/citations/19930084187/downloads/19930084187.pdf
5. NIST, SI conversion factors: https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
6. Abaqus nonlinear cantilever benchmarks: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-nlgeocantilever.htm
7. Jönköping University mechanics beam-deflection reference: https://mechanics.ju.se/SolidMechanics/EnergyMethods.html
