import numpy as np
from scipy.linalg import eig
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, eigsh
from .config import CFG
from .frame3d import element_matrices, local_geometric_stiffness, section_props


def control_shape(xnorm, controls):
    xp = np.array([0.0,0.2,0.4,0.6,0.8,1.0])
    zp = np.array([0.0,*controls,0.0])
    return np.interp(xnorm, xp, zp)


def unpack_design(d, cfg=CFG):
    d = np.asarray(d, dtype=float)
    radii = d[:cfg.n_elem]
    shape = d[cfg.n_elem:cfg.n_elem+4]
    return radii, shape


def initial_design(cfg=CFG):
    return np.r_[np.full(cfg.n_elem, cfg.radius_initial), np.zeros(4)]


def bounds(cfg=CFG):
    b = [(cfg.radius_min,cfg.radius_max)]*cfg.n_elem
    b += [(-cfg.shape_bound,cfg.shape_bound)]*4
    return b


def nodes_from_design(d, cfg=CFG):
    _, shape = unpack_design(d,cfg)
    x = np.linspace(0.0,cfg.arm_length,cfg.n_elem+1)
    xn = x/cfg.arm_length
    z = control_shape(xn,shape)
    y = np.zeros_like(x)
    return np.c_[x,y,z]


def dof_map(i,j):
    return np.r_[np.arange(6*i,6*i+6),np.arange(6*j,6*j+6)]


def fixed_dofs(cfg=CFG):
    return np.arange(6)


def free_dofs(cfg=CFG):
    return np.arange(6*(cfg.n_elem+1))[6:]


def assemble(d, cfg=CFG, include_load=True):
    radii,_ = unpack_design(d,cfg)
    nodes = nodes_from_design(d,cfg)
    ndof = 6*len(nodes)
    rows=[]; cols=[]; kvals=[]; mvals=[]
    F = np.zeros(ndof)
    elems = []
    for e in range(cfg.n_elem):
        ids = dof_map(e,e+1)
        kg,mg,T,kl,L,props = element_matrices(nodes[e],nodes[e+1],radii[e],cfg.wall_thickness,cfg.E,cfg.G,cfg.rho)
        rr,cc=np.meshgrid(ids,ids,indexing='ij')
        rows.extend(rr.ravel()); cols.extend(cc.ravel())
        kvals.extend(kg.ravel()); mvals.extend(mg.ravel())
        A = props[0]
        m = cfg.rho*A*L
        # Lump self-weight equally to the translational z DOF at both nodes.
        if include_load:
            F[6*e+2] -= 0.5*m*cfg.gravity
            F[6*(e+1)+2] -= 0.5*m*cfg.gravity
        elems.append((ids,T,kl,L,props))
    tip = cfg.n_elem
    td=np.arange(6*tip,6*tip+3)
    rows.extend(td); cols.extend(td); kvals.extend([0.0]*3); mvals.extend([cfg.motor_mass]*3)
    K=coo_matrix((kvals,(rows,cols)),shape=(ndof,ndof)).tocsr()
    M=coo_matrix((mvals,(rows,cols)),shape=(ndof,ndof)).tocsr()
    if include_load:
        # Maximum combined coaxial thrust is applied vertically. Counter-rotating
        # motor torque is assumed to cancel. A separate team-selected 2-g axial
        # inertial load provides the N1 pre-buckling compression case.
        F[6*tip+0] += -cfg.maneuver_axial_force
        F[6*tip+2] += cfg.rotor_thrust - cfg.motor_mass*cfg.gravity
        F[6*tip+5] += cfg.motor_torque_z
    return K,M,F,nodes,elems


def solve_static(d,cfg=CFG):
    K,M,F,nodes,elems = assemble(d,cfg,True)
    free = free_dofs(cfg)
    u = np.zeros_like(F)
    u[free] = spsolve(K[free][:,free],F[free])
    return u,K,M,F,nodes,elems


def arm_mass(d,cfg=CFG):
    radii,_ = unpack_design(d,cfg)
    nodes = nodes_from_design(d,cfg)
    total = 0.0
    for e,r in enumerate(radii):
        A,_,_,_ = section_props(r,cfg.wall_thickness)
        L = np.linalg.norm(nodes[e+1]-nodes[e])
        total += cfg.rho*A*L
    return total


def element_end_forces(d,u,cfg=CFG):
    radii,_ = unpack_design(d,cfg)
    nodes = nodes_from_design(d,cfg)
    out=[]
    for e,r in enumerate(radii):
        ids=dof_map(e,e+1)
        _,_,T,kl,L,props=element_matrices(nodes[e],nodes[e+1],r,cfg.wall_thickness,cfg.E,cfg.G,cfg.rho)
        ul=T@u[ids]
        q=kl@ul
        out.append((q,L,props,T,ids))
    return out


def stress_values(d,u,cfg=CFG):
    # Smooth von Mises sampling around the circular tube perimeter. Using
    # signed normal stress avoids non-differentiable absolute-value envelopes.
    radii,_=unpack_design(d,cfg)
    vals=[]
    angles=np.linspace(0.0,2.0*np.pi,8,endpoint=False)
    for e,(q,L,props,T,ids) in enumerate(element_end_forces(d,u,cfg)):
        A,Iy,Iz,J=props
        r=radii[e]
        for end in (0,1):
            base=0 if end==0 else 6
            N=q[base+0]
            Tq=q[base+3]
            My=q[base+4]
            Mz=q[base+5]
            tau=Tq*r/J
            for th in angles:
                y=r*np.cos(th); z=r*np.sin(th)
                sig=N/A - Mz*y/Iz + My*z/Iy
                vm=np.sqrt(sig**2+3.0*tau**2+1e-24)
                vals.append(vm)
    return np.asarray(vals)


def stress_ks_utilization(d,u,cfg=CFG):
    a=stress_values(d,u,cfg)/cfg.stress_limit
    amax=float(np.max(a))
    return amax + np.log(np.sum(np.exp(cfg.ks_rho*(a-amax))))/cfg.ks_rho


def max_stress(d,u,cfg=CFG):
    return float(np.max(stress_values(d,u,cfg)))


def tip_deflection(d,u,cfg=CFG):
    tip=cfg.n_elem
    return float(np.linalg.norm(u[6*tip:6*tip+3]))


def natural_modes(d,cfg=CFG,nmodes=6):
    K,M,F,nodes,elems=assemble(d,cfg,False)
    free=free_dofs(cfg)
    Kr=K[free][:,free].tocsr()
    Mr=M[free][:,free].tocsr()
    k=min(nmodes,Kr.shape[0]-2)
    vals,vecs=eigsh(Kr,k=k,M=Mr,sigma=0.0,which='LM')
    order=np.argsort(vals); vals=np.real(vals[order]); vecs=np.real(vecs[:,order])
    vals=np.maximum(vals,1e-16)
    freqs=np.sqrt(vals)/(2*np.pi)
    modes=[]
    for j in range(len(vals)):
        v=np.zeros(K.shape[0]); v[free]=vecs[:,j]
        # normalize max translation to one for visualization
        trans=v.reshape(-1,6)[:,:3]
        s=np.max(np.linalg.norm(trans,axis=1))
        if s>0: v/=s
        modes.append(v)
    return freqs,modes


def geometric_stiffness(d,u,cfg=CFG):
    radii,_=unpack_design(d,cfg)
    nodes=nodes_from_design(d,cfg)
    ndof=6*len(nodes)
    rows=[]; cols=[]; vals=[]
    axial=[]
    for e,r in enumerate(radii):
        ids=dof_map(e,e+1)
        _,_,T,kl,L,props=element_matrices(nodes[e],nodes[e+1],r,cfg.wall_thickness,cfg.E,cfg.G,cfg.rho)
        ul=T@u[ids]
        A=props[0]
        N=cfg.E*A/L*(ul[6]-ul[0])
        kgl=local_geometric_stiffness(N,L)
        kg=T.T@kgl@T
        rr,cc=np.meshgrid(ids,ids,indexing='ij')
        rows.extend(rr.ravel()); cols.extend(cc.ravel()); vals.extend(kg.ravel())
        axial.append(N)
    KG=coo_matrix((vals,(rows,cols)),shape=(ndof,ndof)).tocsr()
    return KG,np.asarray(axial)


def buckling_eigenpair(d,u,cfg=CFG):
    K,_,_,_,_=assemble(d,cfg,False)
    KG,_=geometric_stiffness(d,u,cfg)
    free=free_dofs(cfg)
    Kr=K[free][:,free].toarray()
    B=(-KG[free][:,free]).toarray()
    vals,vecs=eig(Kr,B)
    good=[]
    for j,val in enumerate(vals):
        if np.isfinite(val) and abs(np.imag(val))<1e-6 and np.real(val)>1e-6:
            good.append((float(np.real(val)),j))
    if not good:
        return 1e9,None,K,KG
    alpha,j=min(good,key=lambda z:z[0])
    phi=np.zeros(K.shape[0]); phi[free]=np.real(vecs[:,j])
    return alpha,phi,K,KG


def buckling_factor(d,u,cfg=CFG):
    alpha,phi,_,_=buckling_eigenpair(d,u,cfg)
    return alpha,phi


def metrics(d,cfg=CFG):
    u,K,M,F,nodes,elems=solve_static(d,cfg)
    freqs,modes=natural_modes(d,cfg,6)
    alpha,bmode=buckling_factor(d,u,cfg)
    return {
        'mass_kg':arm_mass(d,cfg),
        'max_stress_pa':max_stress(d,u,cfg),
        'stress_ks_util':stress_ks_utilization(d,u,cfg),
        'tip_deflection_m':tip_deflection(d,u,cfg),
        'frequency_hz':float(freqs[0]),
        'buckling_factor':float(alpha),
    },u,modes,bmode,nodes
