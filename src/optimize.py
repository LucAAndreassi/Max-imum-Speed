import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from .config import CFG
from .drone_arm import initial_design, bounds, metrics, arm_mass, solve_static, stress_ks_utilization, tip_deflection, natural_modes, buckling_factor
from .sensitivities import matrix_partials, stress_gradient, deflection_gradient, frequency_gradient, buckling_gradient


def fd_explicit_grad(fun,d):
    g=np.zeros(len(d))
    for j in range(len(d)):
        h=2e-6*max(1.0,abs(float(d[j])))
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        g[j]=(fun(dp)-fun(dm))/(2*h)
    return g


class Evaluator:
    def __init__(self,cfg=CFG,include_buckling=True):
        self.cfg=cfg; self.include_buckling=include_buckling
        self.last=None; self.partials=None
        self.initial=initial_design(cfg)
        self.mass0=arm_mass(self.initial,cfg)
        self.scale=np.r_[np.full(cfg.n_elem,0.01),np.full(4,0.01)]
        self.history=[]
        self.iteration=0

    def to_d(self,y): return np.asarray(y)*self.scale
    def to_y(self,d): return np.asarray(d)/self.scale

    def ensure(self,d):
        if self.last is None or not np.array_equal(np.asarray(d),self.last):
            self.last=np.asarray(d).copy(); self.partials=None

    def get_partials(self,d):
        self.ensure(d)
        if self.partials is None:
            self.partials=matrix_partials(d,self.cfg)
        return self.partials

    def obj(self,y):
        d=self.to_d(y)
        return arm_mass(d,self.cfg)/self.mass0

    def obj_jac(self,y):
        d=self.to_d(y)
        g=fd_explicit_grad(lambda x: arm_mass(x,self.cfg),d)/self.mass0
        return g*self.scale

    def c_stress(self,y):
        d=self.to_d(y); u,*_=solve_static(d,self.cfg)
        return 1.0-stress_ks_utilization(d,u,self.cfg)
    def j_stress(self,y):
        # Active-set screening: exact adjoint is evaluated once the constraint
        # is within 25% of its limit. Far-slack constraints have zero KKT
        # multiplier and do not affect the search direction.
        if self.c_stress(y) > 0.25:
            return np.zeros_like(y)
        d=self.to_d(y); g=stress_gradient(d,self.cfg,self.get_partials(d))
        return -g*self.scale

    def c_defl(self,y):
        d=self.to_d(y); u,*_=solve_static(d,self.cfg)
        return 1.0-tip_deflection(d,u,self.cfg)/self.cfg.tip_deflection_limit
    def j_defl(self,y):
        if self.c_defl(y) > 0.25:
            return np.zeros_like(y)
        d=self.to_d(y); g=deflection_gradient(d,self.cfg,self.get_partials(d))
        return -g/self.cfg.tip_deflection_limit*self.scale

    def c_freq(self,y):
        d=self.to_d(y); f,_=natural_modes(d,self.cfg,1)
        return float(f[0]/self.cfg.frequency_limit-1.0)
    def j_freq(self,y):
        d=self.to_d(y); g=frequency_gradient(d,self.cfg,self.get_partials(d))
        return g/self.cfg.frequency_limit*self.scale

    def c_buck(self,y):
        d=self.to_d(y); u,*_=solve_static(d,self.cfg); a,_=buckling_factor(d,u,self.cfg)
        return float(a/self.cfg.buckling_factor_limit-1.0)
    def j_buck(self,y):
        if self.c_buck(y) > 1.0:
            return np.zeros_like(y)
        d=self.to_d(y); g=buckling_gradient(d,self.cfg,self.get_partials(d))
        return g/self.cfg.buckling_factor_limit*self.scale

    def callback(self,yk):
        d=self.to_d(yk)
        met=metrics(d,self.cfg)[0]
        self.iteration+=1
        self.history.append({'iteration':self.iteration,**{k:float(v) for k,v in met.items()}})


def optimize_design(include_buckling=True,cfg=CFG,maxiter=None):
    ev=Evaluator(cfg,include_buckling)
    y0=ev.to_y(ev.initial)
    b=[(lo/s,hi/s) for (lo,hi),s in zip(bounds(cfg),ev.scale)]
    cons=[
        {'type':'ineq','fun':ev.c_stress,'jac':ev.j_stress},
        {'type':'ineq','fun':ev.c_defl,'jac':ev.j_defl},
        {'type':'ineq','fun':ev.c_freq,'jac':ev.j_freq},
    ]
    if include_buckling:
        cons.append({'type':'ineq','fun':ev.c_buck,'jac':ev.j_buck})
    t0=time.perf_counter()
    res=minimize(ev.obj,y0,method='SLSQP',jac=ev.obj_jac,bounds=b,constraints=cons,
                 callback=ev.callback,options={'maxiter':maxiter or cfg.maxiter,'ftol':cfg.ftol,'disp':False})
    elapsed=time.perf_counter()-t0
    d=ev.to_d(res.x)
    return d,res,ev.history,elapsed


def save_history(history,path):
    import csv
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if not history: return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(history[0].keys())); w.writeheader(); w.writerows(history)
