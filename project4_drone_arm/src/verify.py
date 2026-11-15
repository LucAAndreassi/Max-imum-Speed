import time
import numpy as np
from scipy.linalg import eig
from .config import CFG
from .frame3d import element_matrices, local_geometric_stiffness, section_props
from .drone_arm import initial_design, metrics, solve_static, free_dofs, assemble, stress_ks_utilization, tip_deflection, natural_modes, buckling_factor
from .sensitivities import matrix_partials, stress_gradient, deflection_gradient, frequency_gradient, buckling_gradient


def relerr(a,b): return abs(a-b)/max(abs(b),1e-30)


def beam_model(n=20,L=0.30,r=0.012,t=0.0015,load=None,rotate=None):
    E,G,rho=CFG.E,CFG.G,CFG.rho
    if rotate is None: rotate=np.eye(3)
    xyz=[]
    for i in range(n+1): xyz.append(rotate@np.array([L*i/n,0.,0.]))
    xyz=np.asarray(xyz); nd=6*(n+1); K=np.zeros((nd,nd)); M=np.zeros((nd,nd))
    for e in range(n):
        ids=np.r_[np.arange(6*e,6*e+6),np.arange(6*(e+1),6*(e+1)+6)]
        ke,me,*_=element_matrices(xyz[e],xyz[e+1],r,t,E,G,rho)
        K[np.ix_(ids,ids)]+=ke; M[np.ix_(ids,ids)]+=me
    F=np.zeros(nd)
    if load is not None: F[-6:]+=load
    free=np.arange(nd)[6:]
    u=np.zeros(nd)
    if load is not None: u[free]=np.linalg.solve(K[np.ix_(free,free)],F[free])
    return xyz,K,M,F,u,free


def run_verification(cfg=CFG):
    tests=[]
    r=0.012; t=cfg.wall_thickness; L=0.30; A,I,_,J=section_props(r,t)

    # V1 cantilever bending
    P=10.0; tip=np.array([0,0,P,0,0,0],float)
    _,K,M,F,u,free=beam_model(24,L,r,t,tip)
    num=u[-4]; exact=P*L**3/(3*cfg.E*I)
    tests.append({'test':'V1 cantilever tip deflection','numerical':num,'exact':exact,'relative_error':relerr(num,exact)})

    # V2 torsion
    Tq=0.5; tip=np.array([0,0,0,Tq,0,0],float)
    _,K,M,F,u,free=beam_model(12,L,r,t,tip)
    num=u[-3]; exact=Tq*L/(cfg.G*J)
    tests.append({'test':'V2 cantilever torsion','numerical':num,'exact':exact,'relative_error':relerr(num,exact)})

    # V3 first cantilever bending frequency
    _,K,M,F,u,free=beam_model(32,L,r,t,None)
    vals,_=eig(K[np.ix_(free,free)],M[np.ix_(free,free)])
    vals=np.real(vals[np.isfinite(vals)&(np.real(vals)>1e-8)&(np.abs(np.imag(vals))<1e-7)])
    fnum=np.sqrt(np.min(vals))/(2*np.pi)
    beta=1.875104068711961
    fex=beta**2/(2*np.pi*L**2)*np.sqrt(cfg.E*I/(cfg.rho*A))
    tests.append({'test':'V3 cantilever first frequency','numerical':fnum,'exact':fex,'relative_error':relerr(fnum,fex)})

    # V4 cantilever Euler buckling under unit compression
    n=24; xyz,K,M,F,u,free=beam_model(n,L,r,t,None)
    KG=np.zeros_like(K)
    for e in range(n):
        ids=np.r_[np.arange(6*e,6*e+6),np.arange(6*(e+1),6*(e+1)+6)]
        ke,me,Tr,kl,Le,props=element_matrices(xyz[e],xyz[e+1],r,t,cfg.E,cfg.G,cfg.rho)
        kg=Tr.T@local_geometric_stiffness(-1.0,Le)@Tr
        KG[np.ix_(ids,ids)]+=kg
    vals,_=eig(K[np.ix_(free,free)],(-KG)[np.ix_(free,free)])
    vals=np.real(vals[np.isfinite(vals)&(np.real(vals)>1e-8)&(np.abs(np.imag(vals))<1e-7)])
    pnum=np.min(vals); pex=np.pi**2*cfg.E*I/(4*L**2)
    tests.append({'test':'V4 Euler cantilever buckling','numerical':pnum,'exact':pex,'relative_error':relerr(pnum,pex)})

    # V5 global force balance on the actual drone arm.
    d=initial_design(cfg); u,K,M,F,nodes,elems=solve_static(d,cfg)
    R=K@u-F
    applied=np.array([F[0::6].sum(),F[1::6].sum(),F[2::6].sum()])
    reaction=np.array([R[0],R[1],R[2]])
    balance=np.linalg.norm(applied+reaction)/max(np.linalg.norm(applied),1e-30)
    tests.append({'test':'V5 global force balance','numerical':balance,'exact':0.0,'relative_error':balance})

    # V6 directional gradient checks: use a mildly curved valid design so that
    # first-order shape sensitivities are nonzero and can be meaningfully checked.
    d=initial_design(cfg)
    d[-4:]=np.array([0.004,-0.003,0.003,-0.002])
    partials=matrix_partials(d,cfg)
    checks=[]
    responses=[
        ('stress_ks',lambda x: (lambda uu: stress_ks_utilization(x,uu,cfg))(solve_static(x,cfg)[0]),stress_gradient),
        ('tip_deflection',lambda x: (lambda uu: tip_deflection(x,uu,cfg))(solve_static(x,cfg)[0]),deflection_gradient),
        ('frequency',lambda x: natural_modes(x,cfg,1)[0][0],frequency_gradient),
        ('buckling',lambda x: (lambda uu: buckling_factor(x,uu,cfg)[0])(solve_static(x,cfg)[0]),buckling_gradient),
    ]
    directions=[]
    vr=np.zeros_like(d); vr[:cfg.n_elem]=np.linspace(-1,1,cfg.n_elem); vr/=np.linalg.norm(vr); directions.append(('radius',vr))
    vs=np.zeros_like(d); vs[-4:]=np.array([1,-0.5,0.75,-1.0]); vs/=np.linalg.norm(vs); directions.append(('shape',vs))
    for name,fun,gfun in responses:
        g=gfun(d,cfg,partials)
        for kind,v in directions:
            h=1e-4
            fd=(fun(d+h*v)-fun(d-h*v))/(2*h)
            adj=float(g@v)
            err=abs(adj-fd)/max(abs(fd),abs(adj),1e-8)
            checks.append({'response':name,'variable_kind':kind,'adjoint_directional':adj,'finite_difference':float(fd),'relative_error':float(err)})
    maxerr=max(c['relative_error'] for c in checks)
    tests.append({'test':'V6 adjoint directional gradients','numerical':maxerr,'exact':0.0,'relative_error':maxerr})
    return tests,checks
