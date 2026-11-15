import numpy as np
from .config import CFG
from .drone_arm import assemble, solve_static, free_dofs, stress_ks_utilization, tip_deflection, natural_modes, buckling_factor, buckling_eigenpair, geometric_stiffness, nodes_from_design, unpack_design, dof_map


def _step(x):
    return 2e-7*max(1.0,abs(float(x)))


def matrix_partials(d,cfg=CFG):
    d=np.asarray(d,float)
    free=free_dofs(cfg)
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


def _response_partial_u(response,d,u,cfg=CFG):
    gu=np.zeros_like(u)
    free=free_dofs(cfg)
    # Response evaluation is cheap and requires no additional PDE solves.
    for k in free:
        h=1e-9*max(1.0,abs(float(u[k])))
        up=u.copy(); um=u.copy(); up[k]+=h; um[k]-=h
        gu[k]=(response(d,up,cfg)-response(d,um,cfg))/(2*h)
    return gu


def static_adjoint_gradient(response,d,cfg=CFG,partials=None):
    u,K,M,F,_,_=solve_static(d,cfg)
    free=free_dofs(cfg)
    if partials is None:
        partials=matrix_partials(d,cfg)
    gu=_response_partial_u(response,d,u,cfg)
    gd=_response_partial_d(response,d,u,cfg)
    lam=np.zeros_like(u)
    lam[free]=np.linalg.solve(K[np.ix_(free,free)].T,-gu[free])
    for j,(dK,dM,dF) in enumerate(partials):
        gd[j] += lam @ (dK@u-dF)
    return gd


def stress_gradient(d,cfg=CFG,partials=None):
    return static_adjoint_gradient(stress_ks_utilization,d,cfg,partials)


def deflection_gradient(d,cfg=CFG,partials=None):
    return static_adjoint_gradient(tip_deflection,d,cfg,partials)


def frequency_gradient(d,cfg=CFG,partials=None):
    d=np.asarray(d,float)
    if partials is None:
        partials=matrix_partials(d,cfg)
    K,M,_,_,_=assemble(d,cfg,False)
    free=free_dofs(cfg)
    Kr=K[np.ix_(free,free)]; Mr=M[np.ix_(free,free)]
    from scipy.linalg import eigh
    vals,vecs=eigh(Kr,Mr,subset_by_index=[0,0])
    lam=float(vals[0]); phi=vecs[:,0]
    # eigh normalization gives phi^T M phi = 1.
    g=np.zeros(len(d))
    for j,(dK,dM,dF) in enumerate(partials):
        dKr=dK[np.ix_(free,free)]; dMr=dM[np.ix_(free,free)]
        dlam=phi@(dKr-lam*dMr)@phi
        g[j]=dlam/(4*np.pi*np.sqrt(lam))
    return g


def buckling_gradient(d,cfg=CFG,partials=None):
    # Efficient N1 buckling sensitivity: generalized-eigenvalue sensitivity
    # plus one static adjoint solve through the pre-buckling state.
    from .frame3d import element_matrices, local_geometric_stiffness
    d=np.asarray(d,float)
    u,K,M,F,_,_=solve_static(d,cfg)
    alpha,phi,Klin,KG=buckling_eigenpair(d,u,cfg)
    if alpha>1e8 or phi is None:
        return np.zeros(len(d))
    free=free_dofs(cfg)
    B=(-KG)[np.ix_(free,free)]
    phir=phi[free]
    denom=float(phir@B@phir)
    if abs(denom)<1e-14:
        return np.zeros(len(d))
    if partials is None:
        partials=matrix_partials(d,cfg)

    # Partial derivative with respect to the state u through axial forces in KG.
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

    # Direct design partial: eigenvalue sensitivity at fixed pre-buckling u.
    gd=np.zeros(len(d))
    for j,(dK,dM,dF) in enumerate(partials):
        h=_step(d[j])
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        KGp,_=geometric_stiffness(dp,u,cfg)
        KGm,_=geometric_stiffness(dm,u,cfg)
        dKG=(KGp-KGm)/(2*h)
        pe=phi[free]
        dKr=dK[np.ix_(free,free)]
        dKGr=dKG[np.ix_(free,free)]
        gd[j]=float(pe@(dKr+alpha*dKGr)@pe)/denom

    lam=np.zeros_like(u)
    lam[free]=np.linalg.solve(K[np.ix_(free,free)].T,-gu[free])
    for j,(dK,dM,dF) in enumerate(partials):
        gd[j]+=lam@(dK@u-dF)
    return gd
