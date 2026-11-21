import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'

import json
from pathlib import Path
import numpy as np
from .report_builder import build_report
from run_all import make_plots

ROOT=Path(__file__).resolve().parents[1]

def design_from_payload(dct):
    return np.r_[np.asarray(dct['radii_m'],float),np.asarray(dct['shape_controls_m'],float)]

def main():
    payload=json.loads((ROOT/'assets'/'results.json').read_text())
    initial=design_from_payload(payload['initial'])
    opt=design_from_payload(payload['optimized'])
    nob=design_from_payload(payload['optimized_without_buckling'])
    make_plots(payload['optimization_history'],initial,opt,nob,payload['taylor'],payload['verification_details'],payload['benchmark_scaling'],payload['multistart'],payload.get('uniform_feasible'))
    report=build_report(payload,ROOT/'REPORT.md')
    front='---\nlayout: default\ntitle: PDE-Constrained Optimization of a Tubular Drone Arm\n---\n\n'
    (ROOT/'index.md').write_text(front+report)
    print(f'finalized {len(report.splitlines())} report lines')

if __name__=='__main__': main()
