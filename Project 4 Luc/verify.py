"""Deterministic verification ladder. Synthetic tests are not drone results."""
from dataclasses import replace
import time
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import solve_bvp, quad
from scipy.optimize import brentq
from arm import Beam, Elastica, section, hermite


def run_verification():
    rows=[]
    def check(name,computed,expected,tol,unit="1"):
        error=abs(computed-expected)/abs(expected) if expected else abs(computed)
        passed=bool(np.isfinite(error) and error<=tol)
        rows.append(dict(test=name,computed=float(computed),expected=float(expected),
                         error=float(error),error_kind="relative" if expected else "absolute",
                         tolerance=tol,unit=unit,passed=passed))
    # Deliberately synthetic effective beam properties, not CFRP allowables.
    L,E,D,t,mu,P=1.,70e9,.04,.002,.4,100.
    A,I=section(D,t);EI=E*I
    b=Beam(L,E,mu,0.,P,0.,8)
    u,r,_=b.solve(D,t)
    check("V1 tip force: displacement",u[-2],P*L**3/(3*EI),1e-9,"m")
    check("V1 tip force: slope",u[-1],P*L**2/(2*EI),1e-9,"rad")
    check("V1 root moment magnitude",abs(r[1]),P*L,1e-9,"N m")
    udl=replace(b,tip_force=0.,gravity=25.)
    q=-mu*udl.gravity
    uu,_,_=udl.solve(D,t)
    check("V2 uniform load: tip",uu[-2],q*L**4/(8*EI),1e-9,"m")
    ss=replace(udl,elements=8)
    free=np.setdiff1d(np.arange(18),[0,16])
    us,_,_=ss.solve(D,t,free)
    check("V2 simply supported midspan",us[8],5*q*L**4/(384*EI),1e-9,"m")
    # The quartic exact field cannot be exactly represented by cubic Hermite FE.
    displacement_convergence=[]
    for n in [2,4,8,16]:
        bb=replace(udl,elements=n);v,_,_=bb.solve(D,t)
        gx,gw=np.polynomial.legendre.leggauss(5);err=0.
        for e in range(n):
            for a,w in zip((gx+1)/2,gw/2):
                x=(e+a)*L/n
                exact=q*x*x*(6*L*L-4*L*x+x*x)/(24*EI)
                numerical=hermite(a,L/n)@v[2*e:2*e+4]
                err+=w*L/n*(numerical-exact)**2
        displacement_convergence.append(dict(elements=n,h=L/n,error=float(np.sqrt(err))))
    slope=np.polyfit(np.log([a["h"] for a in displacement_convergence]),
                     np.log([a["error"] for a in displacement_convergence]),1)[0]
    check("V3 displacement L2 convergence order",slope,4.,.01)
    beta=brentq(lambda x:np.cos(x)*np.cosh(x)+1,1.,2.)
    exact_f=beta**2/(2*np.pi*L*L)*np.sqrt(EI/mu)
    frequency_convergence=[]
    for n in [2,4,8,16]:
        f=replace(b,elements=n).fe_modes(D,t,1)[0][0]
        frequency_convergence.append(dict(elements=n,h=L/n,computed=float(f),
                                           expected=float(exact_f),error=float(abs(f/exact_f-1))))
    check("V4 cantilever frequency",frequency_convergence[-1]["computed"],exact_f,2e-6,"Hz")
    freq_order=np.log(frequency_convergence[0]["error"]/frequency_convergence[2]["error"])/np.log(4)
    check("V4 frequency convergence order",freq_order,4.,.035)
    Mtip=.8
    ratio=Mtip/(mu*L)
    beta_tip=brentq(lambda x:1+np.cos(x)*np.cosh(x)+ratio*x*(
        np.cos(x)*np.sinh(x)-np.sin(x)*np.cosh(x)),.1,1.875)
    exact_tip=beta_tip**2/(2*np.pi*L*L)*np.sqrt(EI/mu)
    check("V4 frequency including tip mass",replace(b,tip_mass=Mtip,elements=16).fe_modes(D,t,1)[0][0],
          exact_tip,2e-6,"Hz")
    check("V4 NACA eq13 zero tip mass",b.modes(D,t,1)[0][0],exact_f,1e-12,"Hz")
    check("V4 NACA eq13 including tip mass",replace(b,tip_mass=Mtip).modes(D,t,1)[0][0],exact_tip,1e-12,"Hz")
    for mode,(analytical,fe) in enumerate(zip(replace(b,tip_mass=Mtip).modes(D,t,3)[0],
                                            replace(b,tip_mass=Mtip,elements=24).fe_modes(D,t,3)[0]),1):
        check(f"V4 NACA versus FE mode {mode}",fe,analytical,1e-5,"Hz")
    # Only the bending subspace is present: translation and rotation, two modes.
    K,M,_,_=b.matrices(D,t)
    values=eigh(K.toarray(),M.toarray(),eigvals_only=True)
    zero=int(np.sum(abs(values)<1e-8*np.max(abs(values))))
    check("V5 free-free bending rigid modes",zero,2,0.)
    rigid=np.zeros(K.shape[0]);rigid[::2]=1.
    check("V5 consistent mass translation",rigid@M@rigid,mu*L,1e-12,"kg")
    bal=replace(b,tip_mass=.45,gravity=9.80665)
    _,reac,(_,_,F,_,_)=bal.solve(D,t)
    check("V5 force balance residual",reac[0]+sum(F[::2]),0.,1e-7,"N")
    xx=np.linspace(0,L,b.elements+1)
    check("V5 moment balance residual",reac[1]+xx@F[::2]+sum(F[1::2]),0.,1e-7,"N m")
    # Every differentiated response, both diameter and wall thickness.
    stress=lambda dd,tt,v:bal.end_moments(dd,tt,v)[0,0]*(dd/2)/section(dd,tt)[1]
    quantities={"tip_displacement":(lambda d,t:bal.solve(d,t)[0][-2],bal.tip_gradient(D,t)),
                "root_stress":(lambda d,t:stress(d,t,bal.solve(d,t)[0]),bal.scalar_gradient(D,t,stress)),
                "bending_frequency":(lambda d,t:bal.modes(d,t,1)[0][0],bal.frequency_gradient(D,t))}
    sweep=[]
    rng=np.random.default_rng(20261121)
    directions=[np.array([1.,0.]),np.array([0.,1.])]
    directions.extend(v/np.linalg.norm(v) for v in rng.normal(size=(3,2)))
    scales=np.array([D,t]);x=np.array([D,t])
    for name,(fun,gradient) in quantities.items():
        for k,v in enumerate(directions):
            analytic=gradient@(scales*v)
            errors=[]
            for h in [1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5]:
                fd=(fun(*(x+h*scales*v))-fun(*(x-h*scales*v)))/(2*h)
                error=abs(fd-analytic)/max(abs(analytic),1e-12)
                sweep.append(dict(response=name,direction=k,h=h,adjoint=float(analytic),
                                  finite_difference=float(fd),error=float(error)))
                errors.append(error)
            best=int(np.argmin(errors))
            check(f"V6 {name} gradient direction {k}",
                  [a for a in sweep if a["response"]==name and a["direction"]==k][best]["finite_difference"],
                  analytic,2e-5,{"tip_displacement":"m","root_stress":"Pa","bending_frequency":"Hz"}[name])
    # Independent continuum BVP for nonlinear FE (no FE residual reused).
    nl=Elastica(1.,1.,5.,elements=80)
    theta,records=nl.solve(steps=1)
    mesh=np.linspace(0,1,101)
    sol=solve_bvp(lambda x,y:np.vstack((y[1],-5*np.cos(y[0]))),
                  lambda a,b:np.array([a[0],b[1]]),mesh,
                  np.vstack((1.2*mesh*(2-mesh),2.4*(1-mesh))),tol=1e-10,max_nodes=10000)
    if not sol.success: raise RuntimeError("Independent elastica BVP failed")
    z_exact=quad(lambda s:np.sin(sol.sol(s)[0]),0,1,epsabs=1e-12)[0]
    check("V7 nonlinear tip vs independent BVP",nl.coordinates(theta)[-1,1],z_exact,1e-4,"m")
    # A pure-moment fixture has known constant curvature and exact circle.
    alpha=1.7
    moment=Elastica(1.,1.,0.,elements=12,tip_moment=alpha)
    th,_=moment.solve()
    xy=moment.coordinates(th)[-1]
    check("V7 pure moment circle x",xy[0],np.sin(alpha)/alpha,1e-11,"m")
    check("V7 pure moment circle z",xy[1],(1-np.cos(alpha))/alpha,1e-11,"m")
    v=rng.normal(size=nl.n);v/=np.linalg.norm(v)
    tc=theta.astype(complex);tc[1:]+=1e-25j*v
    numerical=nl.residual(tc)[0].imag[1:]/1e-25
    tangent=nl.residual(theta)[1][1:,1:]@v
    check("V7 consistent tangent relative error",np.linalg.norm(numerical-tangent)/np.linalg.norm(tangent),0.,1e-12)
    h=1e-4
    def val(ei):
        obj=Elastica(1.,ei,5.,elements=80)
        return obj.coordinates(obj.solve(steps=1)[0])[-1,1]
    fd=(val(1+h)-val(1-h))/(2*h)
    check("V8 nonlinear implicit EI gradient",nl.tip_EI_gradient(theta),fd,1e-6,"m / (N m2)")
    # Follower load: independent BVP has the unknown tip angle as a parameter.
    follower=Elastica(1.,1.,-.1,q=-.05,elements=80,follower_force=1.)
    tf,_=follower.solve(steps=5)
    def ode_f(s,y,p):
        return np.vstack((y[1],-np.cos(p[0]-y[0])+(.1+.05*(1-s))*np.cos(y[0]),np.sin(y[0])))
    def bc_f(a,b,p):return np.array([a[0],b[1],a[2],b[0]-p[0]])
    sf=solve_bvp(ode_f,bc_f,mesh,np.vstack((.5*mesh*(2-mesh),1-mesh,.3*mesh)),
                 p=[.5],tol=1e-10,max_nodes=10000)
    if not sf.success:raise RuntimeError("Independent follower BVP failed")
    check("V9 follower tip vs independent BVP",follower.coordinates(tf)[-1,1],sf.sol(1)[2],1e-4,"m")
    tc=tf.astype(complex);tc[1:]+=1e-25j*v
    numerical=follower.residual(tc)[0].imag[1:]/1e-25
    Hf=follower.residual(tf)[1][1:,1:]
    check("V9 follower nonsymmetric tangent",np.linalg.norm(numerical-Hf@v)/np.linalg.norm(Hf@v),0.,1e-12)
    def follower_val(ei):
        obj=Elastica(1.,ei,-.1,q=-.05,elements=80,follower_force=1.)
        return obj.coordinates(obj.solve(steps=5)[0])[-1,1]
    fdf=(follower_val(1+h)-follower_val(1-h))/(2*h)
    check("V9 follower transposed adjoint",follower.tip_EI_gradient(tf),fdf,1e-6,"m / (N m2)")
    nonlinear_sweep=[]
    analytic_f=follower.tip_EI_gradient(tf)
    for step in [1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5]:
        finite=(follower_val(1+step)-follower_val(1-step))/(2*step)
        nonlinear_sweep.append(dict(h=step,adjoint=float(analytic_f),finite_difference=float(finite),
                                    error=float(abs(finite/analytic_f-1))))
    # Quadratic convergence is assessed before roundoff, not from final zeroes.
    rr=records[0]["residuals"]
    orders=[np.log(rr[i+1]/rr[i])/np.log(rr[i]/rr[i-1]) for i in range(1,len(rr)-1)
            if rr[i+1]>1e-11 and rr[i]<.1*rr[i-1] and rr[i-1]>.001]
    check("V7 Newton convergence order",orders[-1] if orders else 0.,2.,.2)
    t0=time.perf_counter()
    for _ in range(5): bal.tip_gradient(D,t)
    adj_ms=(time.perf_counter()-t0)*1000/5
    t0=time.perf_counter()
    for _ in range(5):
        for j in range(2):
            v=np.eye(2)[j]*scales*1e-4
            bal.solve(*(x+v));bal.solve(*(x-v))
    fd_ms=(time.perf_counter()-t0)*1000/5
    return dict(status="pass" if all(r["passed"] for r in rows) else "fail",rows=rows,
                displacement_convergence=displacement_convergence,
                frequency_convergence=frequency_convergence,gradient_sweep=sweep,
                nonlinear_newton=records,nonlinear_reference_z=z_exact,
                nonlinear_gradient_sweep=nonlinear_sweep,
                timing=dict(adjoint_tip_gradient_ms=adj_ms,central_difference_two_variables_ms=fd_ms,
                            note="Measured five repeats, synthetic 8-element beam, two geometric variables. Includes assembly; machine-dependent."))
