"""Path C residual-adaptive collocation, using PDE residuals and no truth.

Start from each completed combined 600-evaluation model. At each round add
the 64 largest residual locations from a fixed dense monitor, preserving all
existing constraints, then train 600 more evaluations. Never select on exact
field/parameter error. All rounds, including unfavorable changes, are saved.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
import oscillating_bl as bl


def refine(data,p,yc,scale,width,rounds=5,max_nfev=600):
    """This function has no access to true kappa or delta."""
    monitor=np.unique(np.r_[np.linspace(0,1,4097),np.linspace(0,min(6*scale,1),4097)])
    traces=[]
    refinement_start=time.perf_counter()
    elapsed=0.
    for step in range(1,rounds+1):
        aa,bb,_,_,a2,b2=bl.neural_fields(p,monitor,width,scale)
        kappa=np.exp(p[-1])
        magnitude=np.hypot(bb-kappa*a2,1-aa-kappa*b2)
        # Greedy high-residual locations with spatial separation. The monitor
        # is not an error field; it depends only on the current neural PDE.
        order=np.argsort(magnitude)[::-1]
        added=[]
        spacing=min(scale/400,1e-4)
        for j in order:
            if min(monitor[j],1-monitor[j]) < spacing:
                continue
            if np.min(abs(yc-monitor[j])) < spacing:
                continue
            if added and np.min(abs(np.asarray(added)-monitor[j])) < spacing:
                continue
            added.append(monitor[j])
            if len(added)==64: break
        yc=np.unique(np.r_[yc,added])
        history=[]
        start=time.perf_counter()
        def fun(z):
            r=bl.residual_jacobian(z,data,yc,width,scale)
            nd=len(data['y'])
            history.append([len(history),time.perf_counter()-start,float(r[:nd]@r[:nd]),
                            float(r[nd:-4]@r[nd:-4]),float(r[-4:]@r[-4:]),np.exp(z[-1])])
            return r
        def jac(z):
            return bl.residual_jacobian(z,data,yc,width,scale,True)[1]
        low,high=np.full(len(p),-np.inf),np.full(len(p),np.inf)
        low[-1],high[-1]=np.log(1e-8),np.log(.2)
        opt=least_squares(fun,p,jac=jac,bounds=(low,high),ftol=1e-9,xtol=1e-9,gtol=1e-9,max_nfev=max_nfev)
        elapsed+=time.perf_counter()-start
        p=opt.x
        aa,bb,_,_,a2,b2=bl.neural_fields(p,monitor,width,scale)
        mag=np.hypot(bb-np.exp(p[-1])*a2,1-aa-np.exp(p[-1])*b2)
        fit={'parameters':p.copy(),'scale':scale,'kappa':float(np.exp(p[-1])),
             'collocation':yc.copy(),'width':width,'history':np.array(history),
             'training_seconds':time.perf_counter()-refinement_start,'nfev':opt.nfev,'optimizer_status':opt.status,
             'monitor_max_residual':float(np.max(mag)),'round':step}
        traces.append(fit)
    return traces


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--results',default='final_results')
    ap.add_argument('--rounds',type=int,default=5)
    args=ap.parse_args()
    out=Path(args.results)
    env=json.loads((out/'environment.json').read_text())
    config=env['arguments']
    rows=[]
    initial={}
    with (out/'metrics.csv').open() as f:
        for r in csv.DictReader(f):
            if r['method']=='combined' and float(r['delta'])==.01 and float(r['noise'])==.02 and int(r['nt'])==16:
                initial[int(r['seed'])]=r
    for seed in config['seeds']:
        case=bl.Case(.01,.02,16,seed)
        data=bl.make_measurements(case)
        path=out/'models'/f'{case.key}_combined.npz'
        with np.load(path) as z:
            p,yc,scale=z['parameters'],z['collocation'],float(z['scale'])
        print(f'Adaptive PDE refinement seed={seed}',flush=True)
        # Truth only re-enters here in evaluation after each completed fit.
        traces=refine(data,p,yc,scale,config['width'],args.rounds,config['max_nfev'])
        for fit in traces:
            metrics,fields=bl.evaluate(case,'combined',fit,config['width'])
            row={**bl.asdict(case),'method':'adaptive','round':fit['round'],**metrics,
                 'training_seconds':float(initial[seed]['training_seconds'])+fit['training_seconds'],
                 'collocation_count':len(fit['collocation']),'monitor_max_residual':fit['monitor_max_residual'],
                 'optimizer_status':int(fit['optimizer_status']),'nfev':fit['nfev']}
            rows.append(row)
            tag=f'{case.key}_adaptive_round{fit["round"]}'
            np.savez_compressed(out/'models'/f'{tag}.npz',parameters=fit['parameters'],
                                scale=scale,width=config['width'],collocation=fit['collocation'])
            (out/'models'/f'{tag}.json').write_text(bl.json_text(row),encoding='utf-8')
            bl.write_csv(out/'histories'/f'{tag}.csv',[dict(zip(
                ['evaluation','seconds','data_mse','pde_mse','boundary_mse','kappa'],h)) for h in fit['history']])
            if seed==config['seeds'][0] and fit['round']==args.rounds:
                np.savez_compressed(out/'fields'/f'{case.key}_adaptive.npz',**fields)
            print({k:row[k] for k in ['seed','round','kappa_relative_error','shear_amplitude_relative_error','heldout_layer_pde_rms','monitor_max_residual']},flush=True)
        bl.write_csv(out/'adaptive_metrics.csv',rows)
    (out/'adaptive_environment.json').write_text(bl.json_text({'arguments':vars(args),
       'main_source_sha256':env['source_sha256'],
       'adaptive_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'selection_rule':'five fixed rounds; no oracle-error selection'}),encoding='utf-8')


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
