import numpy as np
from scipy.linalg import eig, eigh
from .config import CFG
from .frame3d import element_matrices, local_geometric_stiffness, section_props
from .drone_arm import (
    initial_design, solve_static, stress_ks_utilization, tip_deflection,
    natural_modes, buckling_factor,
)
from .sensitivities import (
    matrix_partials, stress_gradient, deflection_gradient,
    frequency_gradient, buckling_gradient,
)


def relerr(a, b):
    return abs(a-b)/max(abs(b), 1e-30)


def rotation_matrix_xyz(ax, ay, az):
    cx, sx = np.cos(ax), np.sin(ax)
    cy, sy = np.cos(ay), np.sin(ay)
    cz, sz = np.cos(az), np.sin(az)
    Rx = np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]], float)
    Ry = np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]], float)
    Rz = np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]], float)
    return Rz @ Ry @ Rx


def beam_model(n=20, L=0.30, r=0.012, t=0.0015, load=None, rotate=None, constrained=True):
    E, G, rho = CFG.E, CFG.G, CFG.rho
    if rotate is None:
        rotate = np.eye(3)
    xyz = np.asarray([rotate @ np.array([L*i/n, 0.0, 0.0]) for i in range(n+1)])
    nd = 6*(n+1)
    K = np.zeros((nd, nd))
    M = np.zeros((nd, nd))
    for e in range(n):
        ids = np.r_[np.arange(6*e,6*e+6), np.arange(6*(e+1),6*(e+1)+6)]
        ke, me, *_ = element_matrices(xyz[e], xyz[e+1], r, t, E, G, rho)
        K[np.ix_(ids,ids)] += ke
        M[np.ix_(ids,ids)] += me
    F = np.zeros(nd)
    if load is not None:
        F[-6:] += load
    free = np.arange(nd)[6:] if constrained else np.arange(nd)
    u = np.zeros(nd)
    if load is not None and constrained:
        u[free] = np.linalg.solve(K[np.ix_(free,free)], F[free])
    return xyz, K, M, F, u, free


def first_frequency(n, L, r, t, cfg=CFG):
    _, K, M, _, _, free = beam_model(n, L, r, t, None)
    vals, _ = eig(K[np.ix_(free,free)], M[np.ix_(free,free)])
    vals = np.real(vals[np.isfinite(vals) & (np.real(vals)>1e-8) & (np.abs(np.imag(vals))<1e-7)])
    return float(np.sqrt(np.min(vals))/(2*np.pi))


def run_verification(cfg=CFG):
    tests = []
    details = {}
    r = 0.012
    t = cfg.wall_thickness
    L = 0.30
    A, I, _, J = section_props(r, t)

    # V1 exact cantilever bending.
    P = 10.0
    tip = np.array([0,0,P,0,0,0], float)
    _, K, M, F, u, free = beam_model(24, L, r, t, tip)
    num = u[-4]
    exact_defl = P*L**3/(3*cfg.E*I)
    tests.append({'test':'V1 cantilever tip deflection','numerical':num,'exact':exact_defl,'relative_error':relerr(num,exact_defl)})

    # V2 exact torsion.
    Tq = 0.5
    tip = np.array([0,0,0,Tq,0,0], float)
    _, K, M, F, u, free = beam_model(12, L, r, t, tip)
    num = u[-3]
    exact_torsion = Tq*L/(cfg.G*J)
    tests.append({'test':'V2 cantilever torsion','numerical':num,'exact':exact_torsion,'relative_error':relerr(num,exact_torsion)})

    # V3 mesh-convergence rate for the first bending frequency.
    beta = 1.875104068711961
    fexact = beta**2/(2*np.pi*L**2)*np.sqrt(cfg.E*I/(cfg.rho*A))
    ns = np.array([1,2,4,8,16], int)
    fvals = np.array([first_frequency(int(n), L, r, t, cfg) for n in ns])
    errs = np.abs(fvals-fexact)/fexact
    hs = L/ns
    valid = errs > 1e-13
    order = float(np.polyfit(np.log(hs[valid]), np.log(errs[valid]), 1)[0])
    details['mesh_convergence'] = {
        'elements': ns.tolist(), 'element_size_m': hs.tolist(),
        'frequency_hz': fvals.tolist(), 'relative_error': errs.tolist(),
        'exact_frequency_hz': float(fexact), 'observed_order': order,
    }
    tests.append({'test':f'V3 frequency mesh convergence (p={order:.2f})',
                  'numerical':float(fvals[-1]),'exact':float(fexact),'relative_error':relerr(fvals[-1],fexact)})

    # V4 arbitrary 3D rotation invariance. Rotate both the beam and its load.
    R = rotation_matrix_xyz(np.deg2rad(23.0), np.deg2rad(-31.0), np.deg2rad(17.0))
    Flocal = np.array([0.0,0.0,P])
    Fglobal = R @ Flocal
    tiprot = np.r_[Fglobal, np.zeros(3)]
    _, Kr, Mr, Fr, ur, freer = beam_model(12, L, r, t, tiprot, rotate=R)
    defl_rot = float(np.linalg.norm(ur[-6:-3]))
    tests.append({'test':'V4 arbitrary 3D rotation invariance',
                  'numerical':defl_rot,'exact':exact_defl,'relative_error':relerr(defl_rot,exact_defl)})
    details['rotation_invariance'] = {'angles_deg':[23.0,-31.0,17.0], 'deflection_m':defl_rot, 'exact_m':float(exact_defl)}

    # V5 free-free rigid body modes: exactly six near-zero generalized eigenvalues.
    _, Kff, Mff, _, _, _ = beam_model(8, L, r, t, None, constrained=False)
    vals = eigh(Kff, Mff, eigvals_only=True)
    ref = float(vals[6]) if len(vals) > 6 and vals[6] > 0 else float(np.max(np.abs(vals)))
    tol = max(1e-4, 1e-8*ref)
    nzero = int(np.sum(np.abs(vals) < tol))
    tests.append({'test':'V5 free-free rigid-body modes','numerical':float(nzero),'exact':6.0,'relative_error':abs(nzero-6)/6.0})
    details['rigid_body_modes'] = {'zero_mode_count':nzero, 'eigenvalue_tolerance':tol, 'lowest_eigenvalues':vals[:8].tolist()}

    # V6 exact first cantilever frequency on a fine mesh.
    fnum = first_frequency(32, L, r, t, cfg)
    tests.append({'test':'V6 cantilever first frequency','numerical':fnum,'exact':fexact,'relative_error':relerr(fnum,fexact)})

    # V7 Euler buckling of a cantilever under unit compression.
    n = 24
    xyz, K, M, F, u, free = beam_model(n, L, r, t, None)
    KG = np.zeros_like(K)
    for e in range(n):
        ids = np.r_[np.arange(6*e,6*e+6), np.arange(6*(e+1),6*(e+1)+6)]
        ke, me, Tr, kl, Le, props = element_matrices(xyz[e],xyz[e+1],r,t,cfg.E,cfg.G,cfg.rho)
        kg = Tr.T @ local_geometric_stiffness(-1.0, Le) @ Tr
        KG[np.ix_(ids,ids)] += kg
    vals, _ = eig(K[np.ix_(free,free)], (-KG)[np.ix_(free,free)])
    vals = np.real(vals[np.isfinite(vals) & (np.real(vals)>1e-8) & (np.abs(np.imag(vals))<1e-7)])
    pnum = float(np.min(vals))
    pex = float(np.pi**2*cfg.E*I/(4*L**2))
    tests.append({'test':'V7 Euler cantilever buckling','numerical':pnum,'exact':pex,'relative_error':relerr(pnum,pex)})

    # V8 global force and moment equilibrium on the actual drone arm.
    d = initial_design(cfg)
    u, K, M, F, nodes, elems = solve_static(d,cfg)
    Rvec = K@u-F
    reaction_force = Rvec[:3]
    reaction_moment = Rvec[3:6]
    applied_force = np.zeros(3)
    applied_moment = np.zeros(3)
    root = nodes[0]
    for i, x in enumerate(nodes):
        fi = F[6*i:6*i+3]
        mi = F[6*i+3:6*i+6]
        applied_force += fi
        applied_moment += np.cross(x-root, fi) + mi
    force_err = float(np.linalg.norm(reaction_force + applied_force)/max(np.linalg.norm(applied_force),1e-30))
    moment_scale = max(np.linalg.norm(applied_moment), cfg.arm_length*np.linalg.norm(applied_force), 1e-30)
    moment_err = float(np.linalg.norm(reaction_moment + applied_moment)/moment_scale)
    balance = max(force_err, moment_err)
    tests.append({'test':'V8 global force and moment balance','numerical':balance,'exact':0.0,'relative_error':balance})
    details['global_balance'] = {'force_relative_error':force_err,'moment_relative_error':moment_err}

    # V9 directional gradient checks. Use a mildly curved valid arm so that
    # first-order shape sensitivities are nonzero and meaningfully checked.
    d = initial_design(cfg)
    d[-4:] = np.array([0.004,-0.003,0.003,-0.002])
    partials = matrix_partials(d,cfg)
    checks = []
    responses = [
        ('stress_ks', lambda x: (lambda uu: stress_ks_utilization(x,uu,cfg))(solve_static(x,cfg)[0]), stress_gradient),
        ('tip_deflection', lambda x: (lambda uu: tip_deflection(x,uu,cfg))(solve_static(x,cfg)[0]), deflection_gradient),
        ('frequency', lambda x: natural_modes(x,cfg,1)[0][0], frequency_gradient),
        ('buckling', lambda x: (lambda uu: buckling_factor(x,uu,cfg)[0])(solve_static(x,cfg)[0]), buckling_gradient),
    ]
    directions = []
    vr = np.zeros_like(d)
    vr[:cfg.n_elem] = np.linspace(-1,1,cfg.n_elem)
    vr /= np.linalg.norm(vr)
    directions.append(('radius',vr))
    vs = np.zeros_like(d)
    vs[-4:] = np.array([1,-0.5,0.75,-1.0])
    vs /= np.linalg.norm(vs)
    directions.append(('shape',vs))
    for name, fun, gfun in responses:
        g = gfun(d,cfg,partials)
        for kind, v in directions:
            h = 1e-4
            fd = (fun(d+h*v)-fun(d-h*v))/(2*h)
            adj = float(g@v)
            err = abs(adj-fd)/max(abs(fd),abs(adj),1e-8)
            checks.append({'response':name,'variable_kind':kind,'adjoint_directional':adj,
                           'finite_difference':float(fd),'relative_error':float(err)})
    maxerr = max(c['relative_error'] for c in checks)
    tests.append({'test':'V9 adjoint directional gradients','numerical':maxerr,'exact':0.0,'relative_error':maxerr})
    return tests, checks, details
