import numpy as np
from scipy.sparse.linalg import spsolve, eigsh
from .config import CFG
from .drone_arm import (
    assemble, solve_static, free_dofs, stress_ks_utilization, tip_deflection,
    buckling_eigenpair, geometric_stiffness, nodes_from_design, unpack_design, dof_map,
)
from .frame3d import element_matrices, local_geometric_stiffness


def _step(x):
    return 2e-7*max(1.0,abs(float(x)))


def matrix_partials(d,cfg=CFG):
    d=np.asarray(d,float)
    K0,M0,F0,_,_=assemble(d,cfg,True)
    out=[]
    for j in range(len(d)):
        h=_step(d[j])
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        Kp,Mp,Fp,_,_=assemble(dp,cfg,True)
        Km,Mm,Fm,_,_=assemble(dm,cfg,True)
        out.append(((Kp-Km)/(2*h),(Mp-Mm)/(2*h),(Fp-Fm)/(2*h)))
    return out


def _response_partial_d(response,d,u,cfg=CFG):
    d=np.asarray(d,float)
    g=np.zeros(len(d))
    for j in range(len(d)):
        h=_step(d[j])
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        g[j]=(response(dp,u,cfg)-response(dm,u,cfg))/(2*h)
    return g


def _static_adjoint_from_partials(d,u,K,gu,gd,partials):
    free=free_dofs(CFG) if len(u)==6*(CFG.n_elem+1) else np.arange(len(u))[6:]
    lam=np.zeros_like(u)
    lam[free]=spsolve(K[free][:,free].T.tocsr(),-gu[free])
    out=np.asarray(gd,float).copy()
    for j,(dK,dM,dF) in enumerate(partials):
        out[j]+=lam@(dK@u-dF)
    return out


def _stress_partial_u(d,u,cfg=CFG):
    """Analytical partial d(KS stress utilization)/du at fixed design."""
    radii,_=unpack_design(d,cfg)
    nodes=nodes_from_design(d,cfg)
    entries=[]
    vm_values=[]
    angles=np.linspace(0.0,2.0*np.pi,8,endpoint=False)
    for e,r in enumerate(radii):
        ids=dof_map(e,e+1)
        _,_,T,kl,L,props=element_matrices(
            nodes[e],nodes[e+1],r,cfg.wall_thickness,cfg.E,cfg.G,cfg.rho)
        A,Iy,Iz,J=props
        q=kl@(T@u[ids])
        for end in (0,1):
            base=0 if end==0 else 6
            N=q[base+0]; Tq=q[base+3]; My=q[base+4]; Mz=q[base+5]
            tau=Tq*r/J
            for th in angles:
                y=r*np.cos(th); z=r*np.sin(th)
                sig=N/A - Mz*y/Iz + My*z/Iy
                vm=np.sqrt(sig**2+3.0*tau**2+1e-24)
                dq=np.zeros(12)
                dq[base+0]=(sig/vm)/A
                dq[base+3]=(3.0*tau/vm)*(r/J)
                dq[base+4]=(sig/vm)*(z/Iy)
                dq[base+5]=(sig/vm)*(-y/Iz)
                vm_values.append(vm)
                entries.append((ids,T,kl,dq))
    a=np.asarray(vm_values)/cfg.stress_limit
    amax=float(np.max(a))
    ex=np.exp(cfg.ks_rho*(a-amax))
    weights=ex/np.sum(ex)
    gu=np.zeros_like(u)
    for w,(ids,T,kl,dq) in zip(weights,entries):
        # q = kl T u_e, and dKS/dvm = softmax_weight / stress_limit.
        dks_dq=(w/cfg.stress_limit)*dq
        gu[ids]+=T.T@kl.T@dks_dq
    return gu


def stress_gradient(d,cfg=CFG,partials=None):
    d=np.asarray(d,float)
    u,K,M,F,_,_=solve_static(d,cfg)
    if partials is None:
        partials=matrix_partials(d,cfg)
    gu=_stress_partial_u(d,u,cfg)
    gd=_response_partial_d(stress_ks_utilization,d,u,cfg)
    free=free_dofs(cfg)
    lam=np.zeros_like(u)
    lam[free]=spsolve(K[free][:,free].T.tocsr(),-gu[free])
    for j,(dK,dM,dF) in enumerate(partials):
        gd[j]+=lam@(dK@u-dF)
    return gd


def deflection_gradient(d,cfg=CFG,partials=None):
    d=np.asarray(d,float)
    u,K,M,F,_,_=solve_static(d,cfg)
    if partials is None:
        partials=matrix_partials(d,cfg)
    tip=cfg.n_elem
    v=u[6*tip:6*tip+3]
    n=np.linalg.norm(v)
    gu=np.zeros_like(u)
    if n>0:
        gu[6*tip:6*tip+3]=v/n
    gd=_response_partial_d(tip_deflection,d,u,cfg)
    free=free_dofs(cfg)
    lam=np.zeros_like(u)
    lam[free]=spsolve(K[free][:,free].T.tocsr(),-gu[free])
    for j,(dK,dM,dF) in enumerate(partials):
        gd[j]+=lam@(dK@u-dF)
    return gd


def frequency_gradient(d,cfg=CFG,partials=None):
    d=np.asarray(d,float)
    if partials is None:
        partials=matrix_partials(d,cfg)
    K,M,_,_,_=assemble(d,cfg,False)
    free=free_dofs(cfg)
    Kr=K[free][:,free].tocsr(); Mr=M[free][:,free].tocsr()
    vals,vecs=eigsh(Kr,k=1,M=Mr,sigma=0.0,which='LM')
    lam=float(vals[0]); phi=np.real(vecs[:,0])
    g=np.zeros(len(d))
    for j,(dK,dM,dF) in enumerate(partials):
        dKr=dK[free][:,free]; dMr=dM[free][:,free]
        dlam=phi@(dKr-lam*dMr)@phi
        g[j]=dlam/(4*np.pi*np.sqrt(lam))
    return g


def buckling_gradient(d,cfg=CFG,partials=None):
    # N1 buckling sensitivity: generalized-eigenvalue sensitivity plus one
    # static adjoint solve through the pre-buckling axial-force state.
    d=np.asarray(d,float)
    u,K,M,F,_,_=solve_static(d,cfg)
    alpha,phi,Klin,KG=buckling_eigenpair(d,u,cfg)
    if alpha>1e8 or phi is None:
        return np.zeros(len(d))
    free=free_dofs(cfg)
    B=(-KG[free][:,free])
    phir=phi[free]
    denom=float(phir@B@phir)
    if abs(denom)<1e-14:
        return np.zeros(len(d))
    if partials is None:
        partials=matrix_partials(d,cfg)

    # q_u = partial alpha / partial u through the axial-force dependence of KG.
    gu=np.zeros_like(u)
    radii,_=unpack_design(d,cfg)
    nodes=nodes_from_design(d,cfg)
    for e,r in enumerate(radii):
        ids=dof_map(e,e+1)
        _,_,T,kl,L,props=element_matrices(nodes[e],nodes[e+1],r,cfg.wall_thickness,cfg.E,cfg.G,cfg.rho)
        A=props[0]
        dN_du=(cfg.E*A/L)*(T[6,:]-T[0,:])
        dKGdN=T.T@local_geometric_stiffness(1.0,L)@T
        pe=phi[ids]
        coeff=alpha*float(pe@dKGdN@pe)/denom
        gu[ids]+=coeff*dN_du

    # Direct design derivative at fixed pre-buckling displacement.
    gd=np.zeros(len(d))
    for j,(dK,dM,dF) in enumerate(partials):
        h=_step(d[j])
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        KGp,_=geometric_stiffness(dp,u,cfg)
        KGm,_=geometric_stiffness(dm,u,cfg)
        dKG=(KGp-KGm)/(2*h)
        pe=phi[free]
        dKr=dK[free][:,free]
        dKGr=dKG[free][:,free]
        gd[j]=float(pe@(dKr+alpha*dKGr)@pe)/denom

    lam=np.zeros_like(u)
    lam[free]=spsolve(K[free][:,free].T.tocsr(),-gu[free])
    for j,(dK,dM,dF) in enumerate(partials):
        gd[j]+=lam@(dK@u-dF)
    return gd
