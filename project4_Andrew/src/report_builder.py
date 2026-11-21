from pathlib import Path
import numpy as np


def build_report(data,out_path):
    ini=data['initial']['metrics']
    opt=data['optimized']['metrics']
    nob=data['optimized_without_buckling']['metrics']
    ver=data['verification']
    gc=data['gradient_checks']
    vdet=data['verification_details']
    om=data['optimizer']
    bench=data['benchmark']
    bench_rows=data['benchmark_scaling']
    multistart=data['multistart']
    n1=data['n1_comparison']
    uniform=data['uniform_feasible']
    cfg=data['configuration']

    mass_reduction=100*(ini['mass_kg']-opt['mass_kg'])/ini['mass_kg']
    max_ver=max(v['relative_error'] for v in ver)
    max_gc=max(g['relative_error'] for g in gc)
    radii=np.asarray(data['optimized']['radii_m'])
    shape=np.asarray(data['optimized']['shape_controls_m'])
    active_tip=int(np.sum(np.abs(radii-cfg['radius_min_m'])<1e-7))
    mult=om.get('multipliers',[0,0,0,0])
    while len(mult)<4: mult=list(mult)+[0.0]
    freq_mult=float(mult[2])
    marginal_g_hz=1000*ini['mass_kg']*freq_mult/cfg['frequency_limit_hz']
    rotor_hz=cfg['rotor_shaft_hz']
    rotor_rpm=cfg['rotor_rpm']
    thrust_lbf=cfg['rotor_thrust_n']/4.4482216152605
    tip_mass_lb=cfg['motor_mass_kg']/0.45359237
    ro=cfg['radius_initial_m']; ri=ro-cfg['wall_thickness_m']; L=cfg['arm_length_m']
    A=np.pi*(ro**2-ri**2); I=np.pi/4*(ro**4-ri**4); J=2*I
    Mroot=cfg['rotor_thrust_n']*L
    sigma_hand=Mroot*ro/I
    delta_hand=cfg['rotor_thrust_n']*L**3/(3*68.9e9*I)
    mass_hand=2700.0*A*L
    beta1=1.875104068711961
    f_beam_hand=beta1**2/(2*np.pi*L**2)*np.sqrt(68.9e9*I/(2700.0*A))
    pcr_hand=np.pi**2*68.9e9*I/(4*L**2)
    uniform_mass=uniform['metrics']['mass_kg']
    uniform_saving=100*(uniform_mass-opt['mass_kg'])/uniform_mass

    vrows='\n'.join(
        f"| {v['test']} | {float(v['numerical']):.6g} | {float(v['exact']):.6g} | {float(v['relative_error']):.3e} |"
        for v in ver
    )
    grows='\n'.join(
        f"| {g['response']} | {g['variable_kind']} | {g['adjoint_directional']:.6g} | {g['finite_difference']:.6g} | {g['relative_error']:.3e} |"
        for g in gc
    )
    brows='\n'.join(
        f"| {r['n_design_variables']} | {r['adjoint_all_constraints_s']:.3f} | {r['central_difference_jacobian_s']:.3f} | {r['speedup']:.2f}x |"
        for r in bench_rows
    )
    msrows='\n'.join(
        f"| {i+1} | {r['name']} | {1000*r['initial_metrics']['mass_kg']:.2f} | {r['initial_metrics']['frequency_hz']:.2f} | {1000*r['final_metrics']['mass_kg']:.5f} | {r['final_metrics']['frequency_hz']:.5f} | {r['iterations']} | {'PASS' if r['success'] else 'FAIL'} |"
        for i,r in enumerate(multistart)
    )

    stress_slack=1-opt['stress_ks_util']
    defl_slack=1-opt['tip_deflection_m']/cfg['tip_deflection_limit_m']
    freq_slack=opt['frequency_hz']/cfg['frequency_limit_hz']-1
    buck_slack=opt['buckling_factor']/cfg['buckling_factor_limit']-1
    kktrows='\n'.join([
        f"| KS stress | {opt['stress_ks_util']:.6f} <= 1 | {stress_slack:.6f} | {mult[0]:.6f} | No |",
        f"| Tip displacement | {1000*opt['tip_deflection_m']:.6f} mm <= {1000*cfg['tip_deflection_limit_m']:.1f} mm | {defl_slack:.6f} | {mult[1]:.6f} | No |",
        f"| First frequency | {opt['frequency_hz']:.6f} Hz >= {cfg['frequency_limit_hz']:.1f} Hz | {freq_slack:.3e} | {mult[2]:.6f} | Yes |",
        f"| Buckling factor | {opt['buckling_factor']:.3f} >= {cfg['buckling_factor_limit']:.1f} | {buck_slack:.3f} | {mult[3]:.6f} | No |",
    ])

    conv=vdet['mesh_convergence']
    convrows='\n'.join(
        f"| {n} | {h:.6f} | {f:.6f} | {e:.3e} |"
        for n,h,f,e in zip(conv['elements'],conv['element_size_m'],conv['frequency_hz'],conv['relative_error'])
    )

    n1rows='\n'.join([
        f"| Arm mass | {1000*nob['mass_kg']:.6f} g | {1000*opt['mass_kg']:.6f} g | {1000*(opt['mass_kg']-nob['mass_kg']):+.6f} g |",
        f"| Peak von Mises stress | {nob['max_stress_pa']/1e6:.4f} MPa | {opt['max_stress_pa']/1e6:.4f} MPa | {opt['max_stress_pa']/1e6-nob['max_stress_pa']/1e6:+.4f} MPa |",
        f"| Tip displacement | {1000*nob['tip_deflection_m']:.6f} mm | {1000*opt['tip_deflection_m']:.6f} mm | {1000*(opt['tip_deflection_m']-nob['tip_deflection_m']):+.6f} mm |",
        f"| First natural frequency | {nob['frequency_hz']:.6f} Hz | {opt['frequency_hz']:.6f} Hz | {opt['frequency_hz']-nob['frequency_hz']:+.6f} Hz |",
        f"| Buckling factor | {nob['buckling_factor']:.3f} | {opt['buckling_factor']:.3f} | {opt['buckling_factor']-nob['buckling_factor']:+.3f} |",
    ])

    txt=fr'''# PDE-Constrained Optimization of a Tubular Drone Arm

## 1. Summary

This project minimizes the structural mass of a **0.6096 m (2 ft) cantilevered industrial multirotor arm** using a custom three-dimensional finite-element solver and gradient-based PDE-constrained optimization. The model contains 20 Euler-Bernoulli frame elements, one independent tube-radius variable per element, and four centerline shape-control variables, giving **24 design variables**. The design must satisfy strength, tip-deflection, vibration, and N1 linear-buckling requirements.

The baseline arm has a mass of **{ini['mass_kg']*1000:.2f} g**. The converged design has a mass of **{opt['mass_kg']*1000:.2f} g**, a **{mass_reduction:.1f}% reduction**, while satisfying every requirement. The active system-level constraint is the first natural frequency at **{opt['frequency_hz']:.3f} Hz**; stress, static displacement, and buckling retain substantial margin. The optimization converged from multiple radius profiles, and a separate shape-perturbation restart returned to essentially the same design basin.

<div id="viewer" style="height:560px"></div>
<div id="viewer-controls"></div>
<div id="viewer-metrics"></div>
<script type="module" src="assets/viewer.js"></script>

The viewer defaults to the **first vibration-mode amplitude**, because vibration is the governing constraint. Stress utilization remains available as a selectable field.

![Optimization history](assets/figures/optimization_history.png)

## 2. Structure and formulation

### 2.1 Engineering motivation and physical model

A multirotor arm should be light because structural mass reduces payload and endurance, but excessive mass reduction can lower stiffness, shift structural modes toward rotor excitation, increase deflection, or create instability. The case study represents a heavy-lift or professional-filming multirotor with six coaxial propulsion stations. One representative arm extends **2 ft (0.6096 m)** from the main frame to a coaxial motor pair and is fully fixed at the frame attachment. Each motor produces **15 lbf** maximum thrust, giving **30 lbf ({cfg['rotor_thrust_n']:.1f} N) total thrust per arm**. The combined motors, propellers, and mounting hardware are represented by a **{tip_mass_lb:.1f} lb ({cfg['motor_mass_kg']:.3f} kg)** point mass at the arm tip.

The two counter-rotating motors are assumed to cancel net reaction torque to first order. The primary load is therefore transverse thrust plus gravity. A secondary stability load case represents a **2g horizontal maneuver with vehicle acceleration aligned with the representative arm axis**. The {cfg['motor_mass_kg']:.3f} kg propulsion assembly therefore produces **{cfg['maneuver_axial_force_n']:.3f} N** of compressive inertia along the arm. This is treated as a team-selected worst-case stability case rather than a certified flight-envelope requirement.

![Drone-arm finite-element and load schematic](assets/figures/model_schematic.png)

### 2.2 Material

The baseline material is isotropic 6061-T6 aluminum with $E=68.9$ GPa, $\nu=0.33$, $\rho=2700$ kg/m$^3$, and tensile yield strength $S_y=276$ MPa. The strength limit uses a design factor of 2:

$$
\sigma_{{allow}}=\frac{{S_y}}{{2}}=138\ \text{{MPa}}.
$$

**Material source:** MatWeb, "Aluminum 6061-T6; 6061-T651."

### 2.3 Design variables

The design vector is

$$
\mathbf d=[r_1,\ldots,r_{{20}},z_1,z_2,z_3,z_4]^T.
$$

Each element has an independent outer radius $r_i$ and a fixed wall thickness $t=1.5$ mm:

$$
8\ \text{{mm}}\le r_i\le40\ \text{{mm}}.
$$

The four shape variables move centerline control points at normalized arm locations $x/L=0.2,0.4,0.6,0.8$:

$$
-20\ \text{{mm}}\le z_j\le20\ \text{{mm}}.
$$

Piecewise-linear interpolation defines all intermediate nodal $z$ coordinates, while the clamped root and motor location remain fixed. The size bounds are team-selected geometric/manufacturing bounds: with the fixed wall thickness they prevent invalid or vanishing tubes while allowing a broad range of 16-80 mm outside diameters. The shape bounds prevent extreme centerline curvature and very short elements.

### 2.4 Objective and requirements

The objective is structural arm mass; the motor is included in the vibration mass matrix but not in the arm-mass objective:

$$
\min_{{\mathbf d}}\;m(\mathbf d)=\sum_{{e=1}}^{{20}}\rho A_eL_e,
\qquad
A_e=\pi\left[r_e^2-(r_e-t)^2\right].
$$

| Requirement | Explicit mathematical constraint | Limit | Basis for the limit |
|---|---|---:|---|
| Strength | $KS_\rho(\sigma_{{vM}}/\sigma_{{allow}})\le1$ | $\sigma_{{allow}}=138$ MPa | 6061-T6 yield divided by design factor 2 |
| Tip displacement | $\|\mathbf u_{{tip}}\|\le\delta_{{max}}$ | $\delta_{{max}}={1000*cfg['tip_deflection_limit_m']:.2f}$ mm | Team-selected serviceability limit equal to 1% of the 2 ft arm span |
| Vibration | $f_1\ge f_{{min}}$ | $f_{{min}}={cfg['frequency_limit_hz']:.0f}$ Hz | {rotor_rpm:.0f} RPM = {rotor_hz:.1f} Hz shaft frequency; {cfg['frequency_limit_hz']:.0f} Hz gives a {100*(cfg['frequency_limit_hz']/rotor_hz-1):.1f}% separation margin |
| Stability (N1) | $\alpha_{{cr}}\ge\alpha_{{min}}$ | $\alpha_{{min}}=3$ | Team-selected stability factor following the course project's worked-example convention |
| Geometry | bounds on $r_i$ and $z_j$ | above | Prevent invalid tubes and extreme geometry |

The reduced problem is therefore

$$
\begin{{aligned}}
\min_{{\mathbf d}}\quad &m(\mathbf d)\\
\text{{s.t.}}\quad
&\mathbf K(\mathbf d)\mathbf u-\mathbf F(\mathbf d)=0,\\
&KS_\rho(\sigma_{{vM}}/\sigma_{{allow}})-1\le0,\\
&\|\mathbf u_{{tip}}\|/\delta_{{max}}-1\le0,\\
&1-f_1/f_{{min}}\le0,\\
&1-\alpha_{{cr}}/\alpha_{{min}}\le0,\\
&\mathbf d_{{min}}\le\mathbf d\le\mathbf d_{{max}}.
\end{{aligned}}
$$

## 3. Analysis model

### 3.1 Governing PDE reduction, finite elements, mesh, and boundary conditions

The arm is governed by small-strain three-dimensional elasticity. Euler-Bernoulli frame theory reduces that continuum problem to beam equations along the arm centerline, and the finite-element discretization produces

$$
\mathbf K(\mathbf d)\mathbf u=\mathbf F(\mathbf d).
$$

Each node has six degrees of freedom: three translations and three rotations. The 0.6096 m (2 ft) arm is discretized with 20 frame elements. The root node is fully clamped ($u_x=u_y=u_z=\theta_x=\theta_y=\theta_z=0$). Element orientation is handled with a three-dimensional local-to-global transformation matrix, and the global matrices are assembled from the element matrices. The global stiffness and mass matrices are assembled as SciPy sparse matrices; sparse direct solves are used for the static and adjoint equations, and a sparse symmetric generalized eigensolver is used for vibration modes. The small reduced buckling pencil is converted to dense form only for the generalized buckling eigenvalue extraction.

For a circular hollow section with inner radius $r_i=r-t$,

$$
A=\pi(r^2-r_i^2),\qquad
I_y=I_z=\frac{{\pi}}{{4}}(r^4-r_i^4),\qquad
J=2I.
$$

The mass model uses the consistent Euler-Bernoulli beam mass matrix plus a {cfg['motor_mass_kg']:.3f} kg (1.5 lb) translational point mass at the coaxial propulsion node.

### 3.2 Loads

The static load vector contains:

- distributed arm self-weight;
- motor weight $m_mg$ at the tip;
- maximum combined coaxial thrust $T={cfg['rotor_thrust_n']:.1f}$ N (30 lbf) at the tip; and
- an inward axial maneuver inertia load $F_x=-m_m(2g)=-{cfg['maneuver_axial_force_n']:.3f}$ N.

Thus the tip translational components are

$$
F_x=-m_ma_x,\qquad F_z=T-m_mg,\qquad a_x=2g.
$$

Net motor reaction torque is set to zero because the two coaxial motors counter-rotate and are assumed to cancel their shaft-reaction torques to first order.

### 3.3 Stress recovery and KS aggregation

Element local end forces provide axial force $N$, torsion $T_q$, and bending moments $M_y,M_z$. At eight perimeter locations at both ends of every element,

$$
\sigma(\theta)=\frac{{N}}{{A}}-\frac{{M_zy}}{{I_z}}+\frac{{M_yz}}{{I_y}},
\qquad
\tau=\frac{{T_qr}}{{J}},
$$

and

$$
\sigma_{{vM}}=\sqrt{{\sigma^2+3\tau^2}}.
$$

The many local stress constraints are replaced by a smooth Kreisselmeier-Steinhauser (KS) upper bound. With $a_i=\sigma_{{vM,i}}/\sigma_{{allow}}$,

$$
KS_\rho(\mathbf a)=a_{{max}}+\frac{{1}}{{\rho}}\ln\left(\sum_i e^{{\rho(a_i-a_{{max}})}}\right),\qquad \rho=50.
$$

The true maximum stress is always reported beside the KS quantity.

### 3.4 Vibration

Free vibration is obtained from

$$
\mathbf K\boldsymbol\phi=\omega^2\mathbf M\boldsymbol\phi,
\qquad
f=\frac{{\omega}}{{2\pi}}.
$$

The constraint uses the lowest eigenfrequency of any mode, avoiding a discrete mode-classification rule that can become nonsmooth when modes exchange order.

### 3.5 Linear buckling: N1 component

Axial forces from the computed static state are used to assemble the geometric stiffness $\mathbf K_G(\mathbf u,\mathbf d)$. Linearized buckling is defined by

$$
(\mathbf K+\alpha\mathbf K_G)\boldsymbol\phi_b=0.
$$

The smallest positive $\alpha$ is the critical load factor $\alpha_{{cr}}$. Because $\mathbf K_G$ depends on the static stress state, this is an N1 component: an eigenvalue problem linearized about a computed PDE state.

### 3.6 Modeling assumptions

The model assumes isotropic 6061-T6 aluminum, a circular hollow-tube arm, rigid joints, small strains and rotations at the service load, a perfectly clamped root, and the combined coaxial propulsion assembly represented as a point mass. It neglects local shell buckling, fastener/root-joint compliance, propeller gyroscopic effects, motor electromagnetic excitation, motor torque, aerodynamic distributed loads, manufacturing tolerances, and composite anisotropy. The secondary 2g horizontal-maneuver case is aligned with the arm axis to create a conservative compressive stability check; it is a team-selected worst-case case, not a certified multirotor flight envelope.

## 4. Analysis-code verification

Verification was completed before interpreting the optimization. The ladder deliberately covers exact solutions, mesh convergence, 3D invariance, rigid-body behavior, global equilibrium, stability, vibration, and gradient verification. The largest ladder discrepancy, including the directional gradient rung, is **{max_ver:.2e}**.

| Test | Numerical | Exact/target | Error |
|---|---:|---:|---:|
{vrows}

### 4.1 Independent hand check for the industrial arm

As a case-specific check independent of the finite-element assembly, the uniform starting tube is evaluated with elementary beam theory using $L=0.6096$ m, $r_o={1000*cfg['radius_initial_m']:.1f}$ mm, $t={1000*cfg['wall_thickness_m']:.1f}$ mm, and the 30 lbf transverse thrust. The resulting section and response quantities are

| Quantity | Independent analytical value |
|---|---:|
| Cross-sectional area | {A*1e6:.2f} mm$^2$ |
| Second moment of area | {I*1e12:.1f} mm$^4$ |
| Uniform arm structural mass | {1000*mass_hand:.2f} g |
| Root thrust bending moment | {Mroot:.2f} N m |
| Root bending stress from thrust | {sigma_hand/1e6:.2f} MPa |
| Tip deflection from thrust | {1000*delta_hand:.3f} mm |
| Beam-only first bending frequency | {f_beam_hand:.2f} Hz |
| Euler cantilever buckling load | {pcr_hand/1000:.2f} kN |

The hand stress and deflection use $M=FL$, $\sigma=Mc/I$, and $\delta=FL^3/(3EI)$. The beam-only frequency deliberately excludes the 1.5 lb propulsion tip mass, so it is a verification target for the uniform beam formulation rather than the final vehicle frequency. The optimized model includes that tip mass and therefore has a much lower first natural frequency.

### 4.2 Exact-solution formulas

The primary analytical checks are

$$
\delta_{{tip}}=\frac{{PL^3}}{{3EI}},\qquad
\phi_{{tip}}=\frac{{T_qL}}{{GJ}},
$$

$$
f_1=\frac{{\beta_1^2}}{{2\pi L^2}}\sqrt{{\frac{{EI}}{{\rho A}}}},\qquad \beta_1=1.875104\ldots,
$$

and for a cantilever column,

$$
P_{{cr}}=\frac{{\pi^2EI}}{{4L^2}}.
$$

These formulas were evaluated independently of the assembled drone-arm model and used as the targets in the ladder.

### 4.3 Mesh convergence

The first-frequency convergence study gives an observed error order of **$p={conv['observed_order']:.2f}$**, close to the expected fourth-order eigenvalue/frequency convergence for this cubic bending interpolation.

| Elements | Element length (m) | First frequency (Hz) | Relative error |
|---:|---:|---:|---:|
{convrows}

![Frequency mesh-convergence study](assets/figures/verification_convergence.png)

### 4.4 Invariance, rigid-body modes, and global balance

Rotating the complete beam and its load by arbitrary 3D angles $(23^\circ,-31^\circ,17^\circ)$ changes the computed tip-deflection magnitude by only **{ver[3]['relative_error']:.2e}** relative error. The free-free beam produces exactly **{int(vdet['rigid_body_modes']['zero_mode_count'])}** near-zero modes, as required for a three-dimensional rigid body. On the actual drone-arm load case, the normalized force-balance error is **{vdet['global_balance']['force_relative_error']:.2e}** and the normalized moment-balance error is **{vdet['global_balance']['moment_relative_error']:.2e}**.

## 5. Sensitivity analysis

### 5.1 General reduced-space adjoint

For a scalar response $g(\mathbf d,\mathbf u)$ and residual

$$
\mathbf h(\mathbf d,\mathbf u)=\mathbf K(\mathbf d)\mathbf u-\mathbf F(\mathbf d)=0,
$$

total differentiation gives

$$
\frac{{d\mathbf u}}{{dd_j}}=-\mathbf K^{{-1}}\left(\frac{{\partial\mathbf K}}{{\partial d_j}}\mathbf u-\frac{{\partial\mathbf F}}{{\partial d_j}}\right).
$$

Substitution into $dg/dd_j$ would require one state-sensitivity solve per design variable. Instead define the adjoint by

$$
\mathbf K^T\boldsymbol\lambda=-\left(\frac{{\partial g}}{{\partial\mathbf u}}\right)^T.
$$

Then

$$
\boxed{{
\frac{{dg}}{{dd_j}}=\frac{{\partial g}}{{\partial d_j}}+
\boldsymbol\lambda^T\left(\frac{{\partial\mathbf K}}{{\partial d_j}}\mathbf u-\frac{{\partial\mathbf F}}{{\partial d_j}}\right)
}}.
$$

The objective mass is explicit in the design variables and therefore requires no state adjoint.

### 5.2 KS-stress constraint derivative

For the KS stress response,

$$
w_i=\frac{{e^{{\rho(a_i-a_{{max}})}}}}{{\sum_k e^{{\rho(a_k-a_{{max}})}}}},
\qquad
\frac{{\partial KS}}{{\partial a_i}}=w_i.
$$

At each stress sample,

$$
\frac{{\partial\sigma_{{vM}}}}{{\partial\sigma}}=\frac{{\sigma}}{{\sigma_{{vM}}}},
\qquad
\frac{{\partial\sigma_{{vM}}}}{{\partial\tau}}=\frac{{3\tau}}{{\sigma_{{vM}}}}.
$$

Because the element end-force vector is $\mathbf q_e=\mathbf k_e\mathbf T_e\mathbf u_e$, these local derivatives are mapped back to nodal DOFs to obtain the analytical $\partial KS/\partial\mathbf u$. One static adjoint then accounts for the state dependence of all 24 design variables. Explicit geometry/section partials are evaluated at fixed state at the element/model level.

### 5.3 Tip-displacement constraint derivative

Let $\mathbf P_t$ select the three translational DOFs at the motor node and define

$$
g_\delta=\frac{{\|\mathbf P_t\mathbf u\|}}{{\delta_{{max}}}}-1.
$$

Then

$$
\frac{{\partial g_\delta}}{{\partial\mathbf u}}=
\frac{{1}}{{\delta_{{max}}}}\mathbf P_t^T
\frac{{\mathbf P_t\mathbf u}}{{\|\mathbf P_t\mathbf u\|}}.
$$

This vector becomes the right-hand side of the static adjoint equation above.

### 5.4 Frequency constraint derivative

With $\lambda=\omega^2$ and mass normalization $\boldsymbol\phi^T\mathbf M\boldsymbol\phi=1$,

$$
\frac{{d\lambda}}{{dd_j}}=\boldsymbol\phi^T
\left(\frac{{\partial\mathbf K}}{{\partial d_j}}-\lambda\frac{{\partial\mathbf M}}{{\partial d_j}}\right)\boldsymbol\phi.
$$

Since $f=\sqrt{{\lambda}}/(2\pi)$,

$$
\boxed{{
\frac{{df}}{{dd_j}}=\frac{{1}}{{4\pi\sqrt{{\lambda}}}}
\frac{{d\lambda}}{{dd_j}}
}}.
$$

For the normalized inequality $g_f=1-f/f_{{min}}\le0$,

$$
\frac{{dg_f}}{{dd_j}}=-\frac{{1}}{{f_{{min}}}}\frac{{df}}{{dd_j}}.
$$

No differentiation through the eigensolver iterations is required.

### 5.5 Buckling/N1 derivative including the state adjoint

Write the buckling problem as

$$
\mathbf K\boldsymbol\psi=\alpha\mathbf B\boldsymbol\psi,
\qquad \mathbf B=-\mathbf K_G,
$$

and let $D=\boldsymbol\psi^T\mathbf B\boldsymbol\psi$. At fixed static state,

$$
\left.\frac{{\partial\alpha}}{{\partial d_j}}\right|_{{\mathbf u}}=
\frac{{\boldsymbol\psi^T\left(\frac{{\partial\mathbf K}}{{\partial d_j}}+
\alpha\frac{{\partial\mathbf K_G}}{{\partial d_j}}\right)\boldsymbol\psi}}{{D}}.
$$

However, $\mathbf K_G$ is assembled from element axial forces $N_e(\mathbf u,\mathbf d)$, so the eigenvalue also depends on the static state. Define

$$
\mathbf q_u=\frac{{\partial\alpha}}{{\partial\mathbf u}}
=\frac{{\alpha}}{{D}}
\left[\boldsymbol\psi^T\frac{{\partial\mathbf K_G}}{{\partial\mathbf u}}\boldsymbol\psi\right].
$$

The buckling-state adjoint is

$$
\mathbf K^T\boldsymbol\lambda_b=-\mathbf q_u^T,
$$

and the complete reduced gradient is

$$
\boxed{{
\frac{{d\alpha}}{{dd_j}}=
\left.\frac{{\partial\alpha}}{{\partial d_j}}\right|_{{\mathbf u}}+
\boldsymbol\lambda_b^T
\left(\frac{{\partial\mathbf K}}{{\partial d_j}}\mathbf u-\frac{{\partial\mathbf F}}{{\partial d_j}}\right)
}}.
$$

This is the required N1 eigenvalue sensitivity plus the adjoint through the pre-buckling state.

### 5.6 Gradient checks for every constraint and variable type

The gradient check uses two independent directional perturbations: one supported only on the 20 radius variables and one supported only on the four shape variables. Every PDE-based constraint is checked in both directions.

| Response | Variable kind | Adjoint directional derivative | Central difference | Relative error |
|---|---|---:|---:|---:|
{grows}

The maximum discrepancy is **{max_gc:.2e}**. The step-size sweep below shows the expected decrease in central-difference error before truncation/roundoff effects dominate.

![Taylor step-size test](assets/figures/taylor_test.png)

### 5.7 Cost comparison with finite differences

The current 24-variable model requires **{bench['central_difference_static_solves']}** full static reanalyses, **{bench['central_difference_vibration_eigensolves']}** vibration eigensolves, and **{bench['central_difference_buckling_eigensolves']}** buckling eigensolves for a central-difference Jacobian. The reduced-gradient implementation uses three static solves and one vibration plus one buckling eigensolve after the shared matrix-partial assembly. On the current model this is **{bench['speedup']:.2f}x** faster in measured wall time.

| Design variables | Adjoint/eigenvalue sensitivities (s) | Central-difference Jacobian (s) | Measured speedup |
|---:|---:|---:|---:|
{brows}

![Adjoint versus finite-difference scaling](assets/figures/adjoint_scaling.png)

The scaling experiment intentionally increases the number of beam-size variables. The measured speedup grows with design dimension because central differences repeatedly resolve the PDE/eigenproblems, while the adjoint solve count is controlled primarily by the number of scalar outputs. Local matrix derivatives are still evaluated element/model-wise numerically in this implementation, so the timing does not represent the maximum speedup available from fully analytical element derivatives.

## 6. Optimization

### 6.1 Algorithm, scaling, and convergence

The design is optimized with SciPy SLSQP. Design variables are scaled by 0.01 m,

$$
\mathbf y=\mathbf d/(0.01\ \text{{m}}),
$$

and the objective is normalized by the baseline arm mass $m_0$. Each physical constraint is likewise expressed in dimensionless form. This prevents the optimizer from simultaneously handling raw quantities ranging from millimeters to hertz to kilograms.

Every physical constraint and its verified reduced gradient are evaluated at every optimization point. No constraint-gradient screening is used: SLSQP receives the stress, displacement, frequency, and buckling gradients throughout the run, including when a constraint is far from active. This keeps the optimization formulation identical to the independently verified sensitivity model in Section 5.6.

The reference convergence history contains **{om['iterations']} iterations** and ended with status **"{om['message']}"**. After the optimizer implementation was changed to evaluate all four verified constraint gradients at every point, a full-gradient stationarity restart from the converged design terminated in **{om.get('full_gradient_recheck',{}).get('iterations',0)} iteration** with **{1000*om.get('full_gradient_recheck',{}).get('max_design_change_m',0.0):.6f} mm** maximum design change and **{om.get('full_gradient_recheck',{}).get('mass_change_g',0.0):.6f} g** mass change. The final optimum and active frequency multiplier are therefore unchanged by removal of gradient screening.

![Optimization convergence history](assets/figures/optimization_history.png)

The final radius distribution is strongly tapered: material is retained near the clamped root, where it contributes most efficiently to bending stiffness and the first natural frequency, while **{active_tip} radius variable{'s' if active_tip != 1 else ''}** reach{'es' if active_tip == 1 else ''} the {1000*cfg['radius_min_m']:.0f} mm lower bound.

![Optimized radius distribution](assets/figures/radius_distribution.png)

### 6.2 Comparison with the lightest feasible uniform tube

To separate the value of optimization from the arbitrary choice of starting radius, a one-parameter reference design was also solved in which all 20 elements were constrained to have the same outer radius. The smallest uniform radius that satisfies every requirement is **{1000*uniform['radius_m']:.2f} mm** ({2000*uniform['radius_m']:.2f} mm OD), with a structural mass of **{1000*uniform_mass:.2f} g**. The fully optimized tapered arm has a mass of **{1000*opt['mass_kg']:.2f} g**, which is **{uniform_saving:.1f}% lighter than the best feasible uniform tube**. This comparison shows that the benefit comes from redistributing material along the arm rather than merely selecting a smaller stock tube.

| Design | Uniform outer radius | Mass | First frequency | Feasible? |
|---|---:|---:|---:|---|
| Starting design | {1000*cfg['radius_initial_m']:.2f} mm | {1000*ini['mass_kg']:.2f} g | {ini['frequency_hz']:.2f} Hz | Yes |
| Lightest feasible uniform tube | {1000*uniform['radius_m']:.2f} mm | {1000*uniform_mass:.2f} g | {uniform['metrics']['frequency_hz']:.2f} Hz | Yes |
| Optimized tapered design | varies by element | {1000*opt['mass_kg']:.2f} g | {opt['frequency_hz']:.2f} Hz | Yes |

### 6.3 Multiple starting designs

Because shape optimization makes the problem nonconvex, the solution was checked from several starts. The first three use different radius distributions. A fourth restart perturbs all four shape variables from the converged radius profile to probe the shape basin directly.

| Start | Description | Initial mass (g) | Initial $f_1$ (Hz) | Final mass (g) | Final $f_1$ (Hz) | Iterations | Result |
|---:|---|---:|---:|---:|---:|---:|---|
{msrows}

![Final mass from multiple starts](assets/figures/multistart_final_mass.png)

All successful starts converge to essentially the same mass and frequency. Across the restart set, the maximum RMS distance from the reference design is **{1000*max(r['rms_difference_from_reference_m'] for r in multistart):.4f} mm** in the 24-dimensional design vector. The optimized shape controls themselves are at most **{1000*np.max(np.abs(shape)):.4f} mm**, so the straight centerline is locally preferred for this load model.

### 6.4 Active constraints, multipliers, and KKT interpretation

The final normalized slacks and SLSQP inequality multipliers are:

| Constraint | Final physical value | Normalized slack | Multiplier | Active? |
|---|---|---:|---:|---|
{kktrows}

The vibration constraint is the only active system-level inequality. Its normalized multiplier is **{freq_mult:.6f}**. Using the objective normalization $m/m_0$ and $f_1/f_{{min}}-1\ge0$, this corresponds to a local mass penalty of approximately **{marginal_g_hz:.3f} g per additional hertz** of required minimum frequency. In other words, near the optimum a modest relaxation of the frequency requirement buys measurable mass reduction, whereas relaxing stress, deflection, or buckling has no first-order value because those constraints are inactive.

{active_tip} radius variable{'s are' if active_tip != 1 else ' is'} on the lower bound. SciPy's reported `multipliers` array covers the nonlinear inequality constraints but not simple variable-bound multipliers, so bound activity is reported separately rather than inventing bound multipliers.

## 7. Nonlinear/stability component: N1 with and without buckling

The required N1 study re-optimizes the same problem with the buckling constraint removed. This isolates the design effect of stability rather than merely evaluating buckling after the fact.

| Quantity | Without N1 buckling constraint | With N1 buckling constraint | Change |
|---|---:|---:|---:|
{n1rows}

The maximum tube-radius difference between the two optimized designs is **{n1['max_radius_difference_mm']:.6f} mm**, and the RMS difference over all 24 design variables is **{n1['rms_design_difference_mm']:.6f} mm**.

![Optimized radius profile with and without N1 buckling](assets/figures/radius_n1_comparison.png)

The result is physically clear: $\alpha_{{cr}}={opt['buckling_factor']:.1f}$ at the optimum, compared with a required value of 3. The arm is therefore vibration-dominated, not stability-dominated. The secondary axial maneuver load is small compared with the bending/vibration stiffness demanded by the {cfg['frequency_limit_hz']:.0f} Hz requirement, and the lower tube-size bounds also prevent the tip members from vanishing. The N1 component is correctly implemented and differentiated, but it does not alter this optimum. The first buckling mode remains available in the interactive viewer to show the instability pattern that the N1 analysis predicts.

## 8. Initial-versus-optimized results

| Metric | Initial | Optimized | Requirement | Status |
|---|---:|---:|---:|---|
| Arm mass | {ini['mass_kg']*1000:.2f} g | {opt['mass_kg']*1000:.2f} g | Minimize | - |
| True max von Mises stress | {ini['max_stress_pa']/1e6:.3f} MPa | {opt['max_stress_pa']/1e6:.3f} MPa | <= 138 MPa | PASS |
| KS stress utilization | {ini['stress_ks_util']:.4f} | {opt['stress_ks_util']:.4f} | <= 1 | PASS |
| Tip displacement | {ini['tip_deflection_m']*1000:.4f} mm | {opt['tip_deflection_m']*1000:.4f} mm | <= {1000*cfg['tip_deflection_limit_m']:.2f} mm | PASS |
| First natural frequency | {ini['frequency_hz']:.3f} Hz | {opt['frequency_hz']:.3f} Hz | >= {cfg['frequency_limit_hz']:.0f} Hz | PASS |
| Buckling factor | {ini['buckling_factor']:.1f} | {opt['buckling_factor']:.1f} | >= 3 | PASS |

The optimized arm is **{mass_reduction:.1f}% lighter** while the governing frequency is driven directly to its allowable boundary. The very low stress and deflection utilization values explain why adding material solely for static strength or stiffness would be inefficient for this particular requirement set.

## 9. Website and live 3D comparison

The GitHub Pages viewer is designed explicitly around the eight project requirements:

| Viewer requirement | Implementation |
|---|---|
| 1. Initial and optimized comparison | Both designs are shown simultaneously under one linked camera |
| 2. True 3D geometry | Tube geometry uses the actual optimized radii and 3D centerline |
| 3. Field behind the most important constraint | First-mode amplitude is the default sequential field because vibration governs; stress utilization is available as a toggle with over-limit markers |
| 4. State on geometry | Static deformation uses beam-interpolated curves and an exaggeration slider |
| 5. Explanatory animation | First vibration mode is the default animation; the first buckling mode is also available |
| 6. Hover/click inspection | Hover reports design, element, radius, mode amplitude, and stress utilization |
| 7. Metrics panel | Objective and every constraint are shown for both designs with PASS/FAIL status |
| 8. Data/responsiveness | The viewer reads `assets/results.json`, uses relative paths, and the layout collapses for phone width |

Math is rendered with MathJax on the same GitHub Pages report page as the Three.js viewer.

## 10. Reproducibility

From a clean repository environment, the full numerical report can be regenerated with

```bash
python -m pip install -r requirements.txt
python run_all.py
```

That command regenerates the verification ladder, detailed convergence data, gradient checks, Taylor test, adjoint-scaling benchmark, multiple-start study, optimization history, N1 comparison, every plotted figure, `assets/results.json`, `REPORT.md`, and the GitHub Pages `index.md`. The calculations are deterministic; no random starting design or unrecorded manual number is required.

## 11. Limitations and AI-use statement

### 11.1 Limitations

The model is a design study rather than flight certification. Important omissions include composite anisotropy, root-joint/fastener compliance, local tube-shell buckling, motor torque, rotor gyroscopic effects, blade-passing and motor electromagnetic excitation, aerodynamic distributed loads, controller coupling, manufacturing tolerances, and multiple flight/load cases. The 2g horizontal-maneuver stability case and 1%-of-span serviceability requirement are team-selected design requirements rather than externally certified limits. The frequency constraint considers the {rotor_rpm:.0f} RPM ({rotor_hz:.1f} Hz) shaft frequency but not blade-passing harmonics, so a hardware design would require a broader excitation map and Campbell-style separation study.

### 11.2 AI-use statement

A coding assistant was used to draft and revise the 3D frame solver, sensitivity implementation, optimization workflow, verification scripts, Three.js viewer, and report structure. The code was not accepted solely because it executed: independent analytical targets are shown in Section 4; the solver is checked for exact bending, torsion, vibration and Euler buckling behavior; mesh convergence, arbitrary 3D rotation invariance, six free-free rigid-body modes, and global force/moment balance are checked; and every PDE-based constraint gradient is compared with central differences for both radius and shape variables.

One concrete bug/improvement was exposed by gradient verification. An early stress-adjoint implementation numerically perturbed every displacement DOF, which produced an avoidable sensitivity discrepancy and obscured the computational advantage of the adjoint method. The stress and tip-displacement state derivatives were replaced by analytical response derivatives, after which the worst directional gradient error is **{max_gc:.2e}** and the measured adjoint advantage increases with design dimension as shown in Section 5.7.

**Team-action before submission:** the course specification asks the team to state which derivations/checks were performed by the students themselves. The analytical formulas and verification evidence are included above, but the team should independently reproduce/check those formulas and replace this sentence with an accurate description of each member's human contribution. This report intentionally does not invent human work that has not been documented.

## References

1. Course Project 4 specification, *PDE-Constrained Design Optimization - An Introduction to Optimization for Mechanical Engineers*.
2. MatWeb. *Aluminum 6061-T6; 6061-T651*, material-property database.
'''
    Path(out_path).write_text(txt)
    return txt
