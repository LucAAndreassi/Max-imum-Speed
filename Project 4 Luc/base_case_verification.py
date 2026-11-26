D = 40.0 mm
t = 2.5 mm
L = 0.6096 m
E = 64 GPa
rho = 1520 kg/m³
T = 133.446648457815 N
tip mass = 0.680388555 kg
g = 9.80665 m/s²
linear mesh = 20
nonlinear verification = 80

V1 cantilever tip-load displacement/slope/root moment
V2 simply-supported uniform-load displacement
V3 2/4/8/16 element mesh convergence and fitted order
V4 FE vs NACA tip-mass frequencies
V5 equilibrium residuals and two planar rigid modes
V6 all 18 optimizer gradients + Taylor remainder
V7 40/80/160 element N2 continuum-BVP comparison

return {
    "inputs": {...},
    "ladder": {...},
    "mesh_convergence": [...],
    "matlab_comparison": {...},
    "gradient_checks": [...],
    "taylor_test": {...},
    "nonlinear_refinement": [...],
    "status": "pass"
}
