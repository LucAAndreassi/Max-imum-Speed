"""Independent stored-field audit and supplementary diagnostic figures."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
import oscillating_bl as bl


def read_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--results',default='final_results')
    args=ap.parse_args()
    out=Path(args.results)
    rows=read_rows(out/'metrics.csv')
    env=json.loads((out/'environment.json').read_text())
    assert len(rows)==117
    assert len(list((out/'data').glob('*.npz')))==33
    assert hashlib.sha256(Path(bl.__file__).read_bytes()).hexdigest()==env['source_sha256']
    case=bl.Case(.01,.02,16,env['arguments']['seeds'][0])
    with np.load(out/'fields'/f'{case.key}_vanilla.npz') as z:
        y=z['y']
        err=z['predicted']-z['truth']
        num=np.mean(err*err,axis=1)
        den=np.mean(z['truth']*z['truth'],axis=1)
        manual=np.sqrt(np.sum((num[1:]+num[:-1])*np.diff(y))/np.sum((den[1:]+den[:-1])*np.diff(y)))
    stored=next(float(r['field_relative_l2']) for r in rows if int(r['seed'])==case.seed and float(r['delta'])==.01 and r['method']=='vanilla' and float(r['noise'])==.02 and int(r['nt'])==16)
    assert abs(manual-stored)<1e-13
    for p in out.rglob('*.json'):
        json.loads(p.read_text(),parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    report={'main_method_runs':len(rows),'measurement_datasets':33,
            'hand_reference_velocity':float(1-np.exp(-.5)*np.cos(.5)),
            'independently_recomputed_thin_vanilla_relative_l2':float(manual),
            'stored_relative_l2':stored,'absolute_difference':float(abs(manual-stored)),
            'all_json_standard':True,'main_source_hash_matches':True}

    # Short labels avoid a long unavailable-baseline-residual panel heading.
    for delta in [.1,.01]:
        c=bl.Case(delta,.02,16,case.seed)
        fig,axes=plt.subplots(3,2,figsize=(10,8.5))
        for i,method in enumerate(['baseline','vanilla','combined']):
            with np.load(out/'fields'/f'{c.key}_{method}.npz') as z:
                for j,key in enumerate(['error','residual']):
                    ax=axes[i,j]
                    ax.set(xlabel='Cycle fraction',ylabel='y/L')
                    ax.set_ylim(0,min(6*delta,1))
                    if method=='baseline' and key=='residual':
                        ax.set_facecolor('#f0f0f0')
                        ax.text(.5,.5,'Continuous residual not evaluated\nDiscrete FD equations are solved',ha='center',va='center',transform=ax.transAxes,fontsize=10)
                        ax.set_title('Baseline: residual diagnostic unavailable')
                    else:
                        vmax=max(np.max(abs(z[key])),1e-6)
                        im=ax.pcolormesh(z['t']/bl.TAU,z['y'],z[key],cmap='RdBu_r',shading='auto',vmin=-vmax,vmax=vmax)
                        fig.colorbar(im,ax=ax)
                        ax.set_title(f'{method}: {key} (own color scale)')
        fig.suptitle(f'P4: near-wall fields, delta/L={delta:g}, seed={case.seed}')
        bl.save_figure(fig,out/'figures'/f'P4_fields_d{delta:g}.png')

    colors={'vanilla':'#bd4930','scaling':'#4669bd','combined':'#16826c','adaptive':'#8b3b96'}
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    yy=np.linspace(0,.06,2401)
    for method in colors:
        if method=='adaptive':
            mp=out/'models'/f'{case.key}_adaptive_round5.npz'
        else:
            mp=out/'models'/f'{case.key}_{method}.npz'
        with np.load(mp) as m:
            p,scale,width=m['parameters'],float(m['scale']),int(m['width'])
            aa,bb,_,_,a2,b2=bl.neural_fields(p,yy,width,scale)
        amp=aa-1j*bb
        rc,rs=bb-np.exp(p[-1])*a2,1-aa-np.exp(p[-1])*b2
        axes[0].semilogy(yy,np.maximum(abs(amp-bl.exact_amplitude(yy,case.kappa)),1e-9),color=colors[method],label=method)
        axes[1].semilogy(yy,np.maximum(np.hypot(rc,rs),1e-9),color=colors[method],label=method)
    axes[0].set(xlabel='y/L',ylabel='Complex amplitude error / U0',title='Held-out spatial amplitude error')
    axes[1].set(xlabel='y/L',ylabel='Residual harmonic amplitude',title='Held-out near-wall PDE residual')
    for ax in axes: ax.grid(alpha=.25); ax.legend(fontsize=8)
    bl.save_figure(fig,out/'figures'/'BC_near_wall_profiles.png')

    adaptive=read_rows(out/'adaptive_metrics.csv')
    budget=read_rows(out/'budget_metrics.csv')
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,key,label in zip(axes,['kappa_relative_error','shear_amplitude_relative_error','heldout_layer_pde_rms'],
                            ['Viscosity error (%)','Shear amplitude error (%)','Held-out layer PDE RMS']):
        for method in ['vanilla','scaling','combined']:
            initial=[float(r[key]) for r in rows if int(r['seed']) in env['arguments']['seeds'] and r['method']==method and float(r['delta'])==.01 and float(r['noise'])==.02 and int(r['nt'])==16]
            longer=[float(r[key]) for r in budget if r['method']==method]
            factor=100 if 'relative_error' in key else 1
            ax.plot([600,6000],[factor*np.median(initial),factor*np.median(longer)],'o-',color=colors[method],label=method)
        final=[float(r[key]) for r in adaptive if int(r['round'])==5]
        ax.scatter([3600],[factor*np.median(final)],marker='*',s=110,color=colors['adaptive'],label='adaptive (extra points)')
        ax.set(xlabel='Maximum objective evaluations',ylabel=label)
        ax.set_yscale('log'); ax.grid(alpha=.25)
    axes[0].legend(fontsize=7)
    fig.suptitle('Longer fixed sampling vs residual-adaptive collocation; seed medians')
    bl.save_figure(fig,out/'figures'/'BC_budget_and_adaptive.png')

    fig,axes=plt.subplots(1,2,figsize=(11,4))
    offset=0
    names=['data_mse','pde_mse','boundary_mse']
    for step in range(6):
        suffix='combined' if step==0 else f'adaptive_round{step}'
        h=np.genfromtxt(out/'histories'/f'{case.key}_{suffix}.csv',delimiter=',',names=True)
        for key,color in zip(names,['#4669bd','#bd4930','#16826c']):
            axes[0].semilogy(h['evaluation']+offset,h[key],color=color,label=key if step==0 else None)
        if step: axes[0].axvline(offset,color='gray',lw=.5,alpha=.5)
        offset+=len(h)
    axes[0].set(xlabel='Cumulative objective evaluations',ylabel='Empirical loss',title='P2: initial training and adaptive refinement')
    axes[0].legend(fontsize=8); axes[0].grid(alpha=.25)
    for seed in env['arguments']['seeds']:
        rr=[r for r in adaptive if int(r['seed'])==seed]
        axes[1].semilogy([int(r['round']) for r in rr],[float(r['heldout_layer_pde_rms']) for r in rr],'o-',label=f'seed {seed}')
    axes[1].set(xlabel='Refinement round (64 added points each)',ylabel='Held-out layer PDE RMS',title='P4: all adaptive rounds, including reversals')
    axes[1].legend(fontsize=8); axes[1].grid(alpha=.25)
    bl.save_figure(fig,out/'figures'/'BC_adaptive_losses.png')
    report['budget_control_runs']=len(budget)
    report['adaptive_round_fits']=len(adaptive)
    dense=[]
    for seed in env['arguments']['seeds']:
        c=bl.Case(.01,.02,16,seed)
        for suffix in ['vanilla','scaling','combined','vanilla_long','scaling_long','combined_long','adaptive_round5']:
            with np.load(out/'models'/f'{c.key}_{suffix}.npz') as z:
                yy=np.linspace(0,.06,24001)
                aa,bb,_,_,a2,b2=bl.neural_fields(z['parameters'],yy,int(z['width']),float(z['scale']))
                k=np.exp(z['parameters'][-1])
                rms=np.sqrt(np.trapz(((bb-k*a2)**2+(1-aa-k*b2)**2)/2,yy)/.06)
            dense.append({'seed':seed,'method':suffix,'dense_layer_pde_rms':float(rms),
                          'layer_evaluation_points':len(yy)})
    bl.write_csv(out/'dense_residual_audit.csv',dense)
    numerator=np.median([r['dense_layer_pde_rms'] for r in dense if r['method']=='vanilla_long'])
    denominator=np.median([r['dense_layer_pde_rms'] for r in dense if r['method']=='adaptive_round5'])
    report['dense_residual_adaptive_improvement_factor']=float(numerator/denominator)
    (out/'verification_audit.json').write_text(bl.json_text(report),encoding='utf-8')
    print(bl.json_text(report))


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
