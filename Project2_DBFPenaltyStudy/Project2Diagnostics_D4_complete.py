"""Diagnostics for DBF Project 2: penalty-induced ill-conditioning.

Produces D1-D4 numerical evidence and saves plots/data under project2_outputs/.
D4 compares the same projected-gradient convergence diagnostic before and after
using the augmented-Lagrangian formulation.
"""
from __future__ import annotations
import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Project2Models import BOUNDS, X0, augmented_lagrangian_objective, dbf_objective, penalized_objective, state, weight_constraint

OUT = Path(__file__).resolve().parent / "project2_outputs"
OUT.mkdir(exist_ok=True)
LB = BOUNDS[:,0]; UB=BOUNDS[:,1]; SCALE=UB-LB

def to_x(z): return LB + SCALE*np.asarray(z,dtype=float)
def to_z(x): return (np.asarray(x,dtype=float)-LB)/SCALE

def fd_gradient(fun,z,h=1e-6):
    z=np.asarray(z,dtype=float); g=np.empty_like(z); f0=fun(z)
    for i in range(len(z)):
        step=h*max(1.0,abs(z[i])); zp=z.copy(); zp[i]+=step
        g[i]=(fun(zp)-f0)/step
    return g

def fd_hessian(fun,z,h=2e-4):
    z=np.asarray(z,dtype=float); n=len(z); H=np.zeros((n,n)); f0=fun(z)
    for i in range(n):
        hi=h*max(1.0,abs(z[i])); zp=z.copy(); zm=z.copy(); zp[i]+=hi; zm[i]-=hi
        H[i,i]=(fun(zp)-2*f0+fun(zm))/hi**2
        for j in range(i+1,n):
            hj=h*max(1.0,abs(z[j])); zpp=z.copy(); zpm=z.copy(); zmp=z.copy(); zmm=z.copy()
            zpp[i]+=hi; zpp[j]+=hj; zpm[i]+=hi; zpm[j]-=hj; zmp[i]-=hi; zmp[j]+=hj; zmm[i]-=hi; zmm[j]-=hj
            val=(fun(zpp)-fun(zpm)-fun(zmp)+fun(zmm))/(4*hi*hj); H[i,j]=H[j,i]=val
    return 0.5*(H+H.T)

def positive_spectrum(H,floor=1e-10):
    ev=np.linalg.eigvalsh(H); pos=ev[ev>floor]
    return (ev, math.inf) if len(pos)<2 else (ev,float(pos[-1]/pos[0]))

def jacobi_condition(H):
    d=np.diag(H).copy()
    if np.any(d<=0): return math.inf,np.full_like(H,np.nan),np.array([])
    invsqrt=1/np.sqrt(d); Hs=(invsqrt[:,None]*H)*invsqrt[None,:]
    ev,kappa=positive_spectrum(Hs); return kappa,Hs,ev

def optimize_penalty(rho,z0):
    fun=lambda z: penalized_objective(to_x(z),rho)
    res=minimize(fun,np.clip(z0,0,1),method='L-BFGS-B',bounds=[(0,1)]*6,
                 options={'ftol':1e-13,'gtol':1e-9,'maxiter':3000,'maxls':50})
    return res,to_x(res.x)

def projected_gd(fun,z0,step,tol=1e-6,maxit=200000):
    z=np.clip(np.asarray(z0,float),0,1); hist=[]
    for k in range(maxit):
        f=fun(z); grad=fd_gradient(fun,z); pg=z-np.clip(z-grad,0,1); pgn=float(np.linalg.norm(pg))
        hist.append((k,f,pgn))
        if pgn<=tol: break
        z=np.clip(z-step*grad,0,1)
    return z,np.array(hist)

def augmented_lagrangian(z0,rho0=5.0,max_outer=12):
    z=np.clip(np.asarray(z0,float),0,1); lam=0.0; rho=float(rho0); history=[]; prev_violation=math.inf
    for outer in range(max_outer):
        fun=lambda zz: augmented_lagrangian_objective(to_x(zz),lam,rho)
        res=minimize(fun,z,method='L-BFGS-B',bounds=[(0,1)]*6,
                     options={'ftol':1e-13,'gtol':1e-9,'maxiter':2000,'maxls':50})
        z=res.x; x=to_x(z); g=weight_constraint(x); violation=max(0.0,g); lam=max(0.0,lam+rho*g)
        history.append((outer,dbf_objective(x),g,rho,lam,res.nit))
        if violation<1e-6: break
        if violation>0.5*prev_violation: rho*=5.0
        prev_violation=violation
    return z,np.array(history,dtype=float),rho,lam

def main():
    print('Running smooth DBF Project 2 diagnostics...')
    rhos=[0.1,1.0,10.0,100.0,1000.0,10000.0]; zstart=to_z(X0); rows=[]; optima={}
    for rho in rhos:
        res,xstar=optimize_penalty(rho,zstart); zstar=res.x; zstart=zstar
        fun=lambda z,rr=rho: dbf_objective(to_x(z))+0.5*rr*weight_constraint(to_x(z))**2
        H=fd_hessian(fun,zstar,h=5e-5); ev,kappa=positive_spectrum(H); kj,Hs,evj=jacobi_condition(H)
        s=state(xstar); g=weight_constraint(xstar)
        rows.append({'rho':rho,'success':res.success,'iterations':res.nit,'objective':res.fun,'base_objective':dbf_objective(xstar),
                     'constraint_g':g,'m2_weight':s['m2_takeoff_weight'],'kappa':kappa,'kappa_jacobi':kj,
                     'lambda_min':float(ev[ev>1e-10][0]) if np.any(ev>1e-10) else np.nan,'lambda_max':float(ev[-1])})
        optima[rho]=(xstar,zstar,H,ev)
        print(f'rho={rho:8g} f={res.fun: .6f} g={g:+.3e} kappa={kappa:.3e} jacobi={kj:.3e} nit={res.nit}')
    with (OUT/'conditioning_vs_rho.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    rho_spec=rhos[-1]; _,_,_,ev=optima[rho_spec]
    plt.figure(figsize=(7,4.5)); plt.semilogy(np.arange(1,len(ev)+1),np.maximum(np.abs(ev),1e-14),'o-')
    plt.xlabel('Eigenvalue index'); plt.ylabel('|Hessian eigenvalue|'); plt.title(f'D1: Penalized DBF Hessian spectrum (rho={rho_spec:g})')
    plt.grid(True,alpha=.25); plt.tight_layout(); plt.savefig(OUT/'D1_hessian_spectrum.png',dpi=180); plt.close()

    rr=np.array([r['rho'] for r in rows]); kk=np.array([r['kappa'] for r in rows]); kj=np.array([r['kappa_jacobi'] for r in rows])
    plt.figure(figsize=(7,4.5)); plt.loglog(rr,kk,'o-',label='Original Hessian'); plt.loglog(rr,kj,'s--',label='After Jacobi rescaling')
    plt.xlabel('Penalty weight rho'); plt.ylabel('Condition number kappa'); plt.title('D2: Intrinsic ill-conditioning test')
    plt.grid(True,which='both',alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'D2_kappa_vs_rho.png',dpi=180); plt.close()

    gd_rhos=[0.1,1.0,10.0]; plt.figure(figsize=(7,4.5)); gd_summary=[]
    z_common=to_z(np.array([7.0,6.0,10.0,2.0,100.0,90.0])); base_center=to_z(np.array([6.0,6.0,10.0,2.0,100.0,90.0]))
    direction=z_common-base_center
    for rho in gd_rhos:
        xstar,zstar,Hstar,evstar=optima[rho]; positive=evstar[evstar>1e-10]; L,mu=positive[-1],positive[0]; step=2/(L+mu)
        fun=lambda z,rr=rho: penalized_objective(to_x(z),rr)
        z_gd0=np.clip(zstar+0.08*direction/max(np.linalg.norm(direction),1e-12),0,1)
        zend,hist=projected_gd(fun,z_gd0,step,tol=1e-5,maxit=15000); f_ref=fun(zstar); gap=np.maximum(hist[:,1]-f_ref,1e-16)
        plt.semilogy(hist[:,0],gap,label=f'rho={rho:g}'); gd_summary.append((rho,len(hist),hist[-1,2],fun(zend)-f_ref))
    plt.xlabel('Iteration'); plt.ylabel('Objective gap'); plt.title('D3: Penalty weight slows projected gradient descent')
    plt.grid(True,alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'D3_gradient_descent_convergence.png',dpi=180); plt.close()
    with (OUT/'gradient_descent_summary.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['rho','iterations','final_projected_grad','final_objective_gap']); w.writerows(gd_summary)

    # D4: augmented Lagrangian and local conditioning.
    z_seed=to_z(np.array([7.0,6.0,10.0,2.0,100.0,90.0]))
    z_al,hist_al,rho_final,lam_final=augmented_lagrangian(z_seed,rho0=5.0); x_al=to_x(z_al); g_al=weight_constraint(x_al)
    x_pen,z_pen,H_pen,ev_pen=optima[10000.0]
    fun_pen_local=lambda z: dbf_objective(to_x(z))+0.5*10000.0*weight_constraint(to_x(z))**2
    fun_al_local=lambda z: (dbf_objective(to_x(z))+0.5*rho_final*(weight_constraint(to_x(z))+lam_final/rho_final)**2-0.5*lam_final**2/rho_final)
    H_al=fd_hessian(fun_al_local,z_al,h=5e-5); ev_al,k_al=positive_spectrum(H_al); k_pen=positive_spectrum(H_pen)[1]
    with (OUT/'augmented_lagrangian_history.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['outer_iteration','base_objective','constraint_g','rho','lambda','inner_iterations']); w.writerows(hist_al)
    plt.figure(figsize=(7,4.5)); plt.semilogy(hist_al[:,0],np.maximum(np.abs(hist_al[:,2]),1e-12),'o-')
    plt.xlabel('Augmented-Lagrangian outer iteration'); plt.ylabel('|weight-constraint residual|')
    plt.title('D4: Augmented Lagrangian enforces the weight constraint'); plt.grid(True,alpha=.25); plt.tight_layout()
    plt.savefig(OUT/'D4_augmented_lagrangian_constraint.png',dpi=180); plt.close()

    # D4 required before/after convergence evidence: same projected-gradient
    # diagnostic on the local fixed-penalty and final AL objectives.
    pen_pos=ev_pen[ev_pen>1e-10]; al_pos=ev_al[ev_al>1e-10]
    step_pen=2.0/(pen_pos[-1]+pen_pos[0]); step_al=2.0/(al_pos[-1]+al_pos[0])
    d=direction/max(np.linalg.norm(direction),1e-12)
    z_pen0=np.clip(z_pen+0.08*d,0,1); z_al0=np.clip(z_al+0.08*d,0,1)
    maxit_d4=20000; tol_d4=1e-5
    zend_pen,hist_pen=projected_gd(fun_pen_local,z_pen0,step_pen,tol=tol_d4,maxit=maxit_d4)
    zend_al,hist_al_gd=projected_gd(fun_al_local,z_al0,step_al,tol=tol_d4,maxit=maxit_d4)
    fref_pen=fun_pen_local(z_pen); fref_al=fun_al_local(z_al)
    gap_pen=np.maximum(hist_pen[:,1]-fref_pen,1e-16); gap_al=np.maximum(hist_al_gd[:,1]-fref_al,1e-16)

    plt.figure(figsize=(7,4.5)); plt.semilogy(hist_pen[:,0],gap_pen,label='Fixed penalty (rho=10,000)')
    plt.semilogy(hist_al_gd[:,0],gap_al,label=f'Augmented Lagrangian (rho={rho_final:g})')
    plt.xlabel('Projected-gradient iteration'); plt.ylabel('Objective gap')
    plt.title('D4: Before/after convergence with the same first-order method')
    plt.grid(True,alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'D4_before_after_convergence.png',dpi=180); plt.close()

    d4_rows=[
        ['fixed_penalty',10000.0,k_pen,weight_constraint(x_pen),len(hist_pen),hist_pen[-1,2],float(gap_pen[-1]),step_pen],
        ['augmented_lagrangian',rho_final,k_al,g_al,len(hist_al_gd),hist_al_gd[-1,2],float(gap_al[-1]),step_al],
    ]
    with (OUT/'d4_before_after_summary.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['formulation','rho','kappa','constraint_g','iterations','final_projected_grad','final_objective_gap','step']); w.writerows(d4_rows)

    print('\nD4 comparison')
    print(f'Fixed penalty rho=10000: g={weight_constraint(x_pen):+.3e}, kappa={k_pen:.3e}, iters={len(hist_pen)}, pg={hist_pen[-1,2]:.3e}, gap={gap_pen[-1]:.3e}')
    print(f'Augmented Lagrangian:      g={g_al:+.3e}, rho_final={rho_final:g}, kappa={k_al:.3e}, iters={len(hist_al_gd)}, pg={hist_al_gd[-1,2]:.3e}, gap={gap_al[-1]:.3e}')
    print('AL design:',np.array2string(x_al,precision=5)); print('Outputs written to',OUT)

if __name__=='__main__': main()
