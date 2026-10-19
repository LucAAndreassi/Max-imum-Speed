"""Independent tenfold-budget control to separate sampling and early stopping.

Run after the full study: python budget_controls.py --results final_results
These runs are labeled separately and never replace the original-budget runs.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
from threadpoolctl import threadpool_limits
import oscillating_bl as bl


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--results',default='final_results')
    ap.add_argument('--max-nfev',type=int,default=6000)
    args=ap.parse_args()
    out=Path(args.results).resolve()
    config=json.loads((out/'environment.json').read_text())['arguments']
    rows=[]
    for i,seed in enumerate(config['seeds']):
        case=bl.Case(.01,.02,16,seed)
        data=bl.make_measurements(case)
        initial=case.kappa*(5 if i%2==0 else .2)
        for variant in ['vanilla','scaling','combined']:
            print(f'Budget control seed={seed}, method={variant}, max_nfev={args.max_nfev}',flush=True)
            fit=bl.train_pinn(data,initial,seed,variant,config['width'],config['collocation'],args.max_nfev)
            metrics,fields=bl.evaluate(case,variant,fit,config['width'])
            tag=f'{case.key}_{variant}_long'
            np.savez_compressed(out/'models'/f'{tag}.npz',parameters=fit['parameters'],
                                scale=fit['scale'],width=config['width'],collocation=fit['collocation'])
            bl.write_csv(out/'histories'/f'{tag}.csv',[dict(zip(
                ['evaluation','seconds','data_mse','pde_mse','boundary_mse','kappa'],h)) for h in fit['history']])
            r=bl.residual_jacobian(fit['parameters'],data,fit['collocation'],config['width'],fit['scale'])
            nd=len(data['y'])
            row={**bl.asdict(case),'method':variant,'max_nfev':args.max_nfev,**metrics,
                 'training_seconds':fit['training_seconds'],'nfev':fit['nfev'],
                 'optimizer_status':fit['optimizer_status'],'optimizer_message':fit['optimizer_message'],
                 'final_data_mse':float(r[:nd]@r[:nd]),'final_pde_mse':float(r[nd:-4]@r[nd:-4]),
                 'final_boundary_mse':float(r[-4:]@r[-4:])}
            rows.append(row)
            (out/'models'/f'{tag}.json').write_text(bl.json_text(row),encoding='utf-8')
            bl.write_csv(out/'budget_metrics.csv',rows)
            print({k:row[k] for k in ['kappa_relative_error','shear_amplitude_relative_error','heldout_layer_pde_rms','training_seconds']},flush=True)
    (out/'budget_environment.json').write_text(bl.json_text({
        'arguments':vars(args),'main_source_sha256':hashlib.sha256(Path(bl.__file__).read_bytes()).hexdigest(),
        'control_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),encoding='utf-8')


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
