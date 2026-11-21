import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'

import argparse, json, time
import numpy as np
from dataclasses import replace
from .config import CFG
from .drone_arm import initial_design, solve_static, stress_ks_utilization, tip_deflection, natural_modes, buckling_factor
from .sensitivities import matrix_partials, stress_gradient, deflection_gradient, frequency_gradient, buckling_gradient


def response_values(d,cfg):
    u,*_=solve_static(d,cfg)
    return np.array([
        stress_ks_utilization(d,u,cfg),
        tip_deflection(d,u,cfg),
        natural_modes(d,cfg,1)[0][0],
        buckling_factor(d,u,cfg)[0],
    ],float)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--n',type=int,required=True); args=ap.parse_args()
    cfg=replace(CFG,n_elem=args.n)
    d=initial_design(cfg); d[-4:]=np.array([0.004,-0.003,0.003,-0.002])
    t0=time.perf_counter(); partials=matrix_partials(d,cfg)
    gs=[stress_gradient(d,cfg,partials),deflection_gradient(d,cfg,partials),frequency_gradient(d,cfg,partials),buckling_gradient(d,cfg,partials)]
    tadj=time.perf_counter()-t0
    t0=time.perf_counter(); J=np.zeros((4,len(d)))
    for j in range(len(d)):
        h=1e-4*max(1.0,abs(float(d[j])))
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        J[:,j]=(response_values(dp,cfg)-response_values(dm,cfg))/(2*h)
    tfd=time.perf_counter()-t0
    errs=[float(np.linalg.norm(g-J[i])/max(np.linalg.norm(J[i]),1e-12)) for i,g in enumerate(gs)]
    nvar=len(d)
    out={
        'n_elements':cfg.n_elem,'n_design_variables':nvar,
        'adjoint_all_constraints_s':float(tadj),'central_difference_jacobian_s':float(tfd),
        'speedup':float(tfd/max(tadj,1e-12)),'relative_gradient_errors':errs,
        'adjoint_static_solves':3,'adjoint_vibration_eigensolves':1,'adjoint_buckling_eigensolves':1,
        'central_difference_static_solves':2*nvar,
        'central_difference_vibration_eigensolves':2*nvar,
        'central_difference_buckling_eigensolves':2*nvar,
    }
    print(json.dumps(out))

if __name__=='__main__': main()
