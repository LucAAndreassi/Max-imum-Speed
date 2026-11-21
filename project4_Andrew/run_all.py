import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'

from pathlib import Path
from dataclasses import replace
import json, time, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt

from src.config import CFG
from src.drone_arm import (
    initial_design, solve_static, stress_ks_utilization, tip_deflection,
    natural_modes, buckling_factor, metrics,
)
from src.sensitivities import (
    matrix_partials, stress_gradient, deflection_gradient,
    frequency_gradient, buckling_gradient,
)
from src.verify import run_verification
from src.optimize import optimize_design, save_history
from src.export_results import write_results
from src.report_builder import build_report

ROOT=Path(__file__).resolve().parent


def response_values(d,cfg=CFG):
    u,*_=solve_static(d,cfg)
    return np.array([
        stress_ks_utilization(d,u,cfg),
        tip_deflection(d,u,cfg),
        natural_modes(d,cfg,1)[0][0],
        buckling_factor(d,u,cfg)[0],
    ],float)


def benchmark_one(cfg):
    d=initial_design(cfg)
    d[-4:]=np.array([0.004,-0.003,0.003,-0.002])
    t0=time.perf_counter()
    partials=matrix_partials(d,cfg)
    gs=[
        stress_gradient(d,cfg,partials),
        deflection_gradient(d,cfg,partials),
        frequency_gradient(d,cfg,partials),
        buckling_gradient(d,cfg,partials),
    ]
    tadj=time.perf_counter()-t0
    t0=time.perf_counter()
    J=np.zeros((4,len(d)))
    for j in range(len(d)):
        h=1e-4*max(1.0,abs(float(d[j])))
        dp=d.copy(); dm=d.copy(); dp[j]+=h; dm[j]-=h
        J[:,j]=(response_values(dp,cfg)-response_values(dm,cfg))/(2*h)
    tfd=time.perf_counter()-t0
    errs=[float(np.linalg.norm(g-J[i])/max(np.linalg.norm(J[i]),1e-12)) for i,g in enumerate(gs)]
    nvar=len(d)
    return {
        'n_elements':cfg.n_elem,
        'n_design_variables':nvar,
        'adjoint_all_constraints_s':float(tadj),
        'central_difference_jacobian_s':float(tfd),
        'speedup':float(tfd/max(tadj,1e-12)),
        'relative_gradient_errors':errs,
        'adjoint_static_solves':3,
        'adjoint_vibration_eigensolves':1,
        'adjoint_buckling_eigensolves':1,
        'central_difference_static_solves':2*nvar,
        'central_difference_vibration_eigensolves':2*nvar,
        'central_difference_buckling_eigensolves':2*nvar,
    }


def benchmark_scaling():
    env=os.environ.copy(); env['OPENBLAS_NUM_THREADS']='1'; env['OMP_NUM_THREADS']='1'; env['MKL_NUM_THREADS']='1'
    def one(n):
        proc=subprocess.run([sys.executable,'-m','src.benchmark_case','--n',str(n)],cwd=ROOT,env=env,capture_output=True,text=True,check=True)
        return json.loads(proc.stdout)
    ns=[10,20,40]
    with ThreadPoolExecutor(max_workers=3) as ex:
        rows=list(ex.map(one,ns))
    return rows


def taylor_data(d):
    partials=matrix_partials(d,CFG)
    grads=[
        stress_gradient(d,CFG,partials),
        deflection_gradient(d,CFG,partials),
        frequency_gradient(d,CFG,partials),
        buckling_gradient(d,CFG,partials),
    ]
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


def run_opt_case(case, reference_path=None):
    env=os.environ.copy()
    env['OPENBLAS_NUM_THREADS']='1'; env['OMP_NUM_THREADS']='1'; env['MKL_NUM_THREADS']='1'
    cmd=[sys.executable,'-m','src.optimize_case','--case',case]
    if reference_path is not None:
        cmd += ['--reference-design',str(reference_path)]
    proc=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True,check=True)
    return json.loads(proc.stdout)


def assemble_multistart(reference_case, uniform15_case, taper_case, shape_case):
    ref=np.asarray(reference_case['final_design'],float)
    named=[
        ('Uniform 35 mm, straight',reference_case),
        ('Uniform 32 mm, straight',uniform15_case),
        ('38-to-30 mm taper, straight',taper_case),
        ('Converged radii + curved shape perturbation',shape_case),
    ]
    rows=[]
    for name,case in named:
        dfin=np.asarray(case['final_design'],float)
        rows.append({
            'name':name,
            'initial_metrics':case['initial_metrics'],
            'final_metrics':case['final_metrics'],
            'success':case['success'],
            'message':case['message'],
            'iterations':case['iterations'],
            'elapsed_s':case['elapsed_s'],
            'multipliers':case.get('multipliers',[]),
            'final_design':case['final_design'],
            'rms_difference_from_reference_m':float(np.sqrt(np.mean((dfin-ref)**2))),
            'max_difference_from_reference_m':float(np.max(np.abs(dfin-ref))),
        })
    return rows


def best_uniform_feasible(cfg=CFG):
    def design(r):
        d=initial_design(cfg); d[:cfg.n_elem]=r; d[-4:]=0.0; return d
    def margins(r):
        m=metrics(design(r),cfg)[0]
        return {
            'stress':1.0-m['stress_ks_util'],
            'deflection':1.0-m['tip_deflection_m']/cfg.tip_deflection_limit,
            'frequency':m['frequency_hz']/cfg.frequency_limit-1.0,
            'buckling':m['buckling_factor']/cfg.buckling_factor_limit-1.0,
        },m
    def governing(r): return min(margins(r)[0].values())
    if governing(cfg.radius_min)>=0: rstar=cfg.radius_min
    else:
        if governing(cfg.radius_max)<0: raise RuntimeError('No feasible uniform tube exists inside radius bounds.')
        rstar=brentq(governing,cfg.radius_min,cfg.radius_max,xtol=1e-11,rtol=1e-10)
    d=design(rstar); mar,met=margins(rstar)
    return {'radius_m':float(rstar),'diameter_m':float(2*rstar),'metrics':{k:float(v) for k,v in met.items()},'normalized_margins':{k:float(v) for k,v in mar.items()},'design':d.tolist()}

def make_plots(history,initial_d,opt_d,no_buck_d,taylor,verification_details,bench_rows,multistart,uniform_feasible=None):
    figdir=ROOT/'assets'/'figures'; figdir.mkdir(parents=True,exist_ok=True)

    if history:
        it=[h['iteration'] for h in history]
        mass=[1000*h['mass_kg'] for h in history]
        freq=[h['frequency_hz'] for h in history]
        fig,ax=plt.subplots(figsize=(7.2,4.2))
        ax.plot(it,mass,marker='o',label='Arm mass (g)')
        ax.set_xlabel('SLSQP iteration'); ax.set_ylabel('Arm mass (g)'); ax.grid(True,alpha=.25)
        ax2=ax.twinx(); ax2.plot(it,freq,marker='s',label='First frequency (Hz)')
        ax2.axhline(CFG.frequency_limit,linestyle='--',label='Frequency limit')
        ax2.set_ylabel('First natural frequency (Hz)')
        lines=ax.get_lines()+ax2.get_lines(); ax.legend(lines,[l.get_label() for l in lines],loc='best')
        fig.tight_layout(); fig.savefig(figdir/'optimization_history.png',dpi=180); plt.close(fig)

    x=np.arange(1,CFG.n_elem+1)
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.plot(x,1000*initial_d[:CFG.n_elem],marker='o',label='Initial')
    ax.plot(x,1000*opt_d[:CFG.n_elem],marker='o',label='Optimized')
    if uniform_feasible is not None:
        ax.axhline(1000*uniform_feasible['radius_m'],linestyle='-.',label='Lightest feasible uniform tube')
    ax.axhline(1000*CFG.radius_min,linestyle='--',label='Lower bound')
    ax.set_xlabel('Element number (root to motor)'); ax.set_ylabel('Outer radius (mm)')
    ax.grid(True,alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figdir/'radius_distribution.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.plot(x,1000*opt_d[:CFG.n_elem],marker='o',label='With N1 buckling')
    ax.plot(x,1000*no_buck_d[:CFG.n_elem],marker='x',label='Without buckling')
    ax.set_xlabel('Element number (root to motor)'); ax.set_ylabel('Outer radius (mm)')
    ax.grid(True,alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figdir/'radius_n1_comparison.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    for name,dct in taylor.items(): ax.loglog(dct['h'],dct['error'],marker='o',label=name)
    ax.set_xlabel('Central-difference step size'); ax.set_ylabel('|FD directional derivative - adjoint|')
    ax.grid(True,which='both',alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figdir/'taylor_test.png',dpi=180); plt.close(fig)

    conv=verification_details['mesh_convergence']
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.loglog(conv['element_size_m'],conv['relative_error'],marker='o')
    ax.set_xlabel('Element length h (m)'); ax.set_ylabel('Relative error in first frequency')
    ax.grid(True,which='both',alpha=.25)
    ax.text(0.04,0.08,f"Observed order p = {conv['observed_order']:.2f}",transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(figdir/'verification_convergence.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    nv=[r['n_design_variables'] for r in bench_rows]
    ta=[r['adjoint_all_constraints_s'] for r in bench_rows]
    tf=[r['central_difference_jacobian_s'] for r in bench_rows]
    ax.plot(nv,ta,marker='o',label='Adjoint / eigenvalue sensitivities')
    ax.plot(nv,tf,marker='s',label='Central-difference Jacobian')
    ax.set_xlabel('Number of design variables'); ax.set_ylabel('Wall-clock time (s)')
    ax.grid(True,alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figdir/'adjoint_scaling.png',dpi=180); plt.close(fig)

    # Model/load schematic for readers who have not taken the course.
    fig,ax=plt.subplots(figsize=(8.0,3.0))
    xs=np.linspace(0,CFG.arm_length,CFG.n_elem+1)
    ax.plot(xs,np.zeros_like(xs),marker='.',linewidth=2)
    ax.plot([0,0],[-0.035,0.035],linewidth=5)
    cp=CFG.arm_length*np.array([0.2,0.4,0.6,0.8])
    ax.scatter(cp,np.zeros_like(cp),s=55,label='Shape-control locations')
    ax.scatter([CFG.arm_length],[0],s=110,label='Coaxial propulsion tip mass')
    ax.annotate('30 lbf total thrust',xy=(CFG.arm_length,0.0),xytext=(CFG.arm_length,0.075),arrowprops=dict(arrowstyle='->'),ha='center')
    ax.annotate('2g horizontal-maneuver inertia (N1)',xy=(CFG.arm_length,0.0),xytext=(CFG.arm_length-0.17,-0.07),arrowprops=dict(arrowstyle='->'),ha='center')
    ax.text(0.003,0.042,'Clamped root')
    ax.text(CFG.arm_length/2,0.018,'20 Euler-Bernoulli frame elements',ha='center')
    ax.set_xlim(-0.02,CFG.arm_length+0.03); ax.set_ylim(-0.09,0.09); ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel('Arm axis x (m)'); ax.set_yticks([]); ax.grid(False); ax.legend(loc='upper left',fontsize=8)
    fig.tight_layout(); fig.savefig(figdir/'model_schematic.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    labels=[r['name'] for r in multistart]
    masses=[1000*r['final_metrics']['mass_kg'] for r in multistart]
    ax.bar(np.arange(len(labels)),masses)
    ax.set_xticks(np.arange(len(labels))); ax.set_xticklabels(['Start 1','Start 2','Start 3','Shape restart'])
    ax.set_ylabel('Final arm mass (g)'); ax.grid(True,axis='y',alpha=.25)
    fig.tight_layout(); fig.savefig(figdir/'multistart_final_mass.png',dpi=180); plt.close(fig)


def main():
    print('step: initialize', flush=True)
    initial=initial_design(CFG)

    print('step: benchmark scaling', flush=True)
    bench_rows=benchmark_scaling()
    bench=next(r for r in bench_rows if r['n_elements']==CFG.n_elem)

    print('step: verification', flush=True)
    verification,gradient_checks,verification_details=run_verification(CFG)

    print('step: optimize reference', flush=True)
    ref_case=run_opt_case('reference')
    opt=np.asarray(ref_case['final_design'],float)
    hist=ref_case['history']
    optimizer_meta={
        'success':ref_case['success'],'message':ref_case['message'],'iterations':ref_case['iterations'],
        'function_evaluations':ref_case['function_evaluations'],'gradient_evaluations':ref_case['gradient_evaluations'],
        'elapsed_s':ref_case['elapsed_s'],'multipliers':ref_case.get('multipliers',[]),
    }
    refpath=ROOT/'results'/'reference_design.json'
    refpath.parent.mkdir(parents=True,exist_ok=True)
    refpath.write_text(json.dumps(ref_case['final_design']))

    print('step: multistart', flush=True)
    # Independent starts run in isolated single-threaded processes so numerical-library
    # state from one case cannot affect another.
    # Run the independent starts sequentially. This is slightly slower than
    # launching nested BLAS-backed subprocesses in parallel, but is much more
    # reliable across laptops, CI runners, and constrained course environments.
    shape_case=run_opt_case('shape',refpath)
    uniform15_case=run_opt_case('uniform15')
    taper_case=run_opt_case('taper')
    multistart=assemble_multistart(ref_case,uniform15_case,taper_case,shape_case)

    print('step: best uniform feasible design', flush=True)
    uniform_feasible=best_uniform_feasible(CFG)
    (ROOT/'results'/'uniform_feasible.json').write_text(json.dumps(uniform_feasible,indent=2))

    print('step: optimize without buckling', flush=True)
    nob_case=run_opt_case('no_buck')
    no_buck=np.asarray(nob_case['final_design'],float)
    optimizer_meta['without_buckling_elapsed_s']=nob_case['elapsed_s']

    verify_d=initial.copy(); verify_d[-4:]=np.array([0.004,-0.003,0.003,-0.002])
    print('step: taylor', flush=True)
    taylor=taylor_data(verify_d)

    print('step: save results', flush=True)
    save_history(hist,ROOT/'results'/'optimization_history.csv')
    (ROOT/'results'/'verification.json').write_text(json.dumps(verification,indent=2))
    (ROOT/'results'/'verification_details.json').write_text(json.dumps(verification_details,indent=2))
    (ROOT/'results'/'gradient_checks.json').write_text(json.dumps(gradient_checks,indent=2))
    (ROOT/'results'/'benchmark.json').write_text(json.dumps(bench,indent=2))
    (ROOT/'results'/'benchmark_scaling.json').write_text(json.dumps(bench_rows,indent=2))
    (ROOT/'results'/'taylor.json').write_text(json.dumps(taylor,indent=2))
    (ROOT/'results'/'multistart.json').write_text(json.dumps(multistart,indent=2))

    opt_met=metrics(opt,CFG)[0]; nob_met=metrics(no_buck,CFG)[0]
    n1={
        'mass_difference_g':float(1000*(opt_met['mass_kg']-nob_met['mass_kg'])),
        'max_radius_difference_mm':float(1000*np.max(np.abs(opt[:CFG.n_elem]-no_buck[:CFG.n_elem]))),
        'rms_design_difference_mm':float(1000*np.sqrt(np.mean((opt-no_buck)**2))),
    }
    (ROOT/'results'/'n1_comparison.json').write_text(json.dumps(n1,indent=2))

    payload=write_results(ROOT/'assets'/'results.json',initial,opt,no_buck,hist,verification,gradient_checks,optimizer_meta,CFG)
    payload['benchmark']=bench
    payload['benchmark_scaling']=bench_rows
    payload['taylor']=taylor
    payload['verification_details']=verification_details
    payload['multistart']=multistart
    payload['n1_comparison']=n1
    payload['uniform_feasible']=uniform_feasible
    (ROOT/'assets'/'results.json').write_text(json.dumps(payload,indent=2))

    print('step: finalize plots and report', flush=True)
    env=os.environ.copy(); env['OPENBLAS_NUM_THREADS']='1'; env['OMP_NUM_THREADS']='1'; env['MKL_NUM_THREADS']='1'
    subprocess.run([sys.executable,'-m','src.finalize_report'],cwd=ROOT,env=env,check=True)

    print('Project regenerated successfully.')
    print(json.dumps({
        'initial':payload['initial']['metrics'],
        'optimized':payload['optimized']['metrics'],
        'benchmark':bench,
        'optimizer':optimizer_meta,
        'n1_comparison':n1,
        'uniform_feasible':uniform_feasible,
        'multistart_final_masses_g':[1000*r['final_metrics']['mass_kg'] for r in multistart],
    },indent=2))

if __name__=='__main__':
    main()
