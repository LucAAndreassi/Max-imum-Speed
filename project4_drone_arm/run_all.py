from pathlib import Path
import json, time
import numpy as np
import matplotlib.pyplot as plt

from src.config import CFG
from src.drone_arm import initial_design, solve_static, stress_ks_utilization, tip_deflection, natural_modes, buckling_factor, metrics
from src.sensitivities import matrix_partials, stress_gradient, deflection_gradient, frequency_gradient, buckling_gradient
from src.verify import run_verification
from src.optimize import optimize_design, save_history
from src.export_results import write_results
from src.report_builder import build_report

ROOT=Path(__file__).resolve().parent


def response_values(d):
    u,*_=solve_static(d,CFG)
    return np.array([
        stress_ks_utilization(d,u,CFG),
        tip_deflection(d,u,CFG),
        natural_modes(d,CFG,1)[0][0],
        buckling_factor(d,u,CFG)[0],
    ],float)


def benchmark(d):
    t0=time.perf_counter()
    partials=matrix_partials(d,CFG)
    gs=[stress_gradient(d,CFG,partials),deflection_gradient(d,CFG,partials),frequency_gradient(d,CFG,partials),buckling_gradient(d,CFG,partials)]
    tadj=time.perf_counter()-t0
    t0=time.perf_counter()
    J=np.zeros((4,len(d)))
    for j in range(len(d)):
        h=1e-4*max(1.0,abs(float(d[j])))
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        J[:,j]=(response_values(dp)-response_values(dm))/(2*h)
    tfd=time.perf_counter()-t0
    errs=[]
    for i,g in enumerate(gs): errs.append(float(np.linalg.norm(g-J[i])/max(np.linalg.norm(J[i]),1e-12)))
    return {'adjoint_all_constraints_s':tadj,'central_difference_jacobian_s':tfd,'speedup':tfd/max(tadj,1e-12),'relative_gradient_errors':errs}


def taylor_data(d):
    partials=matrix_partials(d,CFG)
    grads=[stress_gradient(d,CFG,partials),deflection_gradient(d,CFG,partials),frequency_gradient(d,CFG,partials),buckling_gradient(d,CFG,partials)]
    names=['stress KS','tip deflection','frequency','buckling']
    v=np.r_[np.linspace(-1,1,CFG.n_elem),np.array([0.5,-0.75,1.0,-0.4])]
    v=v/np.linalg.norm(v)
    hs=np.logspace(-3,-8,11)
    out={}
    for i,name in enumerate(names):
        true=float(grads[i]@v); errors=[]
        for h in hs:
            fp=response_values(d+h*v)[i]; fm=response_values(d-h*v)[i]
            fd=(fp-fm)/(2*h)
            errors.append(abs(fd-true))
        out[name]={'h':hs.tolist(),'error':[float(x) for x in errors],'adjoint_directional':true}
    return out


def make_plots(history,initial_d,opt_d,taylor):
    figdir=ROOT/'assets'/'figures'; figdir.mkdir(parents=True,exist_ok=True)
    if history:
        it=[h['iteration'] for h in history]; mass=[1000*h['mass_kg'] for h in history]; freq=[h['frequency_hz'] for h in history]
        fig,ax=plt.subplots(figsize=(7.2,4.2)); ax.plot(it,mass,marker='o',label='Arm mass (g)'); ax.set_xlabel('SLSQP iteration'); ax.set_ylabel('Arm mass (g)'); ax.grid(True,alpha=.25)
        ax2=ax.twinx(); ax2.plot(it,freq,marker='s',label='First frequency (Hz)'); ax2.axhline(CFG.frequency_limit,linestyle='--'); ax2.set_ylabel('First natural frequency (Hz)')
        lines=ax.get_lines()+ax2.get_lines(); ax.legend(lines,[l.get_label() for l in lines],loc='best'); fig.tight_layout(); fig.savefig(figdir/'optimization_history.png',dpi=180); plt.close(fig)
    x=np.arange(1,CFG.n_elem+1)
    fig,ax=plt.subplots(figsize=(7.2,4.2)); ax.plot(x,1000*initial_d[:CFG.n_elem],marker='o',label='Initial'); ax.plot(x,1000*opt_d[:CFG.n_elem],marker='o',label='Optimized'); ax.axhline(1000*CFG.radius_min,linestyle='--',label='Lower bound'); ax.set_xlabel('Element number (root to motor)'); ax.set_ylabel('Outer radius (mm)'); ax.grid(True,alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figdir/'radius_distribution.png',dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    for name,dct in taylor.items(): ax.loglog(dct['h'],dct['error'],marker='o',label=name)
    ax.set_xlabel('Central-difference step size'); ax.set_ylabel('|FD directional derivative - adjoint|'); ax.grid(True,which='both',alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figdir/'taylor_test.png',dpi=180); plt.close(fig)


def main():
    initial=initial_design(CFG)
    verification,gradient_checks=run_verification(CFG)
    verify_d=initial.copy(); verify_d[-4:]=np.array([0.004,-0.003,0.003,-0.002])
    bench=benchmark(verify_d)
    taylor=taylor_data(verify_d)

    no_buck,res0,h0,t0=optimize_design(False,CFG,maxiter=CFG.maxiter)
    opt,res,hist,elapsed=optimize_design(True,CFG,maxiter=CFG.maxiter)
    save_history(hist,ROOT/'results'/'optimization_history.csv')
    (ROOT/'results'/'verification.json').write_text(json.dumps(verification,indent=2))
    (ROOT/'results'/'gradient_checks.json').write_text(json.dumps(gradient_checks,indent=2))
    (ROOT/'results'/'benchmark.json').write_text(json.dumps(bench,indent=2))
    (ROOT/'results'/'taylor.json').write_text(json.dumps(taylor,indent=2))

    optimizer_meta={'success':bool(res.success),'message':str(res.message),'iterations':int(res.nit),'function_evaluations':int(res.nfev),'gradient_evaluations':int(res.njev),'elapsed_s':float(elapsed),'multipliers':[float(x) for x in getattr(res,'multipliers',[])], 'without_buckling_elapsed_s':float(t0)}
    payload=write_results(ROOT/'assets'/'results.json',initial,opt,no_buck,hist,verification,gradient_checks,optimizer_meta,CFG)
    payload['benchmark']=bench; payload['taylor']=taylor
    (ROOT/'assets'/'results.json').write_text(json.dumps(payload,indent=2))
    make_plots(hist,initial,opt,taylor)
    report=build_report(payload,ROOT/'REPORT.md')
    front='---\nlayout: default\ntitle: PDE-Constrained Optimization of a Tubular Drone Arm\n---\n\n'
    (ROOT/'index.md').write_text(front+report)
    print('Project regenerated successfully.')
    print(json.dumps({'initial':payload['initial']['metrics'],'optimized':payload['optimized']['metrics'],'benchmark':bench,'optimizer':optimizer_meta},indent=2))

if __name__=='__main__': main()
