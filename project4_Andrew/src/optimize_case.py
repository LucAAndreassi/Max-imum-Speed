import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'

import argparse, json
from pathlib import Path
import numpy as np
from .config import CFG
from .drone_arm import initial_design, metrics
from .optimize import optimize_design


def make_start(case, reference_path=None):
    if case in ('reference','no_buck'):
        return initial_design(CFG)
    if case=='uniform15':
        d=initial_design(CFG); d[:CFG.n_elem]=0.032; return d
    if case=='taper':
        d=initial_design(CFG); d[:CFG.n_elem]=np.linspace(0.038,0.030,CFG.n_elem); return d
    if case=='shape':
        if not reference_path:
            raise ValueError('shape case requires --reference-design')
        ref=np.asarray(json.loads(Path(reference_path).read_text()),float)
        d=ref.copy(); d[-4:]=np.array([0.0002,-0.00015,0.0001,-0.00005]); return d
    raise ValueError(case)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--case',required=True,choices=['reference','no_buck','uniform15','taper','shape'])
    ap.add_argument('--reference-design')
    args=ap.parse_args()
    d0=make_start(args.case,args.reference_design)
    include_buckling=args.case!='no_buck'
    d,res,hist,elapsed=optimize_design(include_buckling,CFG,initial_d=d0)
    payload={
        'case':args.case,
        'initial_design':d0.tolist(),
        'final_design':d.tolist(),
        'initial_metrics':{k:float(v) for k,v in metrics(d0,CFG)[0].items()},
        'final_metrics':{k:float(v) for k,v in metrics(d,CFG)[0].items()},
        'history':hist,
        'success':bool(res.success),
        'message':str(res.message),
        'iterations':int(res.nit),
        'function_evaluations':int(res.nfev),
        'gradient_evaluations':int(res.njev),
        'elapsed_s':float(elapsed),
        'multipliers':[float(x) for x in getattr(res,'multipliers',[])],
    }
    print(json.dumps(payload))

if __name__=='__main__':
    main()
