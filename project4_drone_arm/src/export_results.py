import json
import numpy as np
from pathlib import Path
from .config import CFG
from .drone_arm import (metrics, unpack_design, nodes_from_design, element_end_forces,
                        stress_values, natural_modes, buckling_factor)
from .frame3d import local_axes


def element_stress_utilization(d,u,cfg=CFG):
    vals=stress_values(d,u,cfg).reshape(cfg.n_elem,16)
    return vals.max(axis=1)/cfg.stress_limit


def sampled_curve(d,u,nsamp=9,cfg=CFG):
    radii,_=unpack_design(d,cfg)
    nodes=nodes_from_design(d,cfg)
    curves=[]
    for e in range(cfg.n_elem):
        x1=nodes[e]; x2=nodes[e+1]
        L,R=local_axes(x1,x2)
        ids=np.r_[np.arange(6*e,6*e+6),np.arange(6*(e+1),6*(e+1)+6)]
        ul=np.zeros(12) if u is None else np.block([R@u[ids[:3]],R@u[ids[3:6]],R@u[ids[6:9]],R@u[ids[9:12]]])
        pts=[]; disp=[]
        for xi in np.linspace(0,1,nsamp):
            H1=1-3*xi**2+2*xi**3
            H2=xi-2*xi**2+xi**3
            H3=3*xi**2-2*xi**3
            H4=-xi**2+xi**3
            ua=(1-xi)*ul[0]+xi*ul[6]
            v=H1*ul[1]+H2*L*ul[5]+H3*ul[7]+H4*L*ul[11]
            w=H1*ul[2]-H2*L*ul[4]+H3*ul[8]-H4*L*ul[10]
            dg=R.T@np.array([ua,v,w])
            p=x1+xi*(x2-x1)
            pts.append(p.tolist()); disp.append(dg.tolist())
        curves.append({'element':e,'radius_m':float(radii[e]),'points':pts,'displacement':disp})
    return curves


def design_payload(d,cfg=CFG):
    met,u,modes,bmode,nodes=metrics(d,cfg)
    util=element_stress_utilization(d,u,cfg)
    mode1=modes[0] if modes else np.zeros_like(u)
    if bmode is None: bmode=np.zeros_like(u)
    return {
        'metrics':{k:float(v) for k,v in met.items()},
        'nodes':nodes.tolist(),
        'radii_m':unpack_design(d,cfg)[0].tolist(),
        'shape_controls_m':unpack_design(d,cfg)[1].tolist(),
        'stress_utilization':util.tolist(),
        'static_curves':sampled_curve(d,u,9,cfg),
        'mode1_curves':sampled_curve(d,mode1,9,cfg),
        'buckling_curves':sampled_curve(d,bmode,9,cfg),
    }


def write_results(path,initial_d,optimized_d,no_buckling_d,history,verification,gradient_checks,optimizer_meta,cfg=CFG):
    payload={
        'configuration':{
            'arm_length_m':cfg.arm_length,'n_elements':cfg.n_elem,
            'wall_thickness_m':cfg.wall_thickness,'stress_limit_pa':cfg.stress_limit,
            'tip_deflection_limit_m':cfg.tip_deflection_limit,'frequency_limit_hz':cfg.frequency_limit,
            'buckling_factor_limit':cfg.buckling_factor_limit,'motor_mass_kg':cfg.motor_mass,
            'rotor_thrust_n':cfg.rotor_thrust,'thrust_tilt_deg':cfg.thrust_tilt_deg,
        },
        'initial':design_payload(initial_d,cfg),
        'optimized':design_payload(optimized_d,cfg),
        'optimized_without_buckling':design_payload(no_buckling_d,cfg),
        'optimization_history':history,
        'verification':verification,
        'gradient_checks':gradient_checks,
        'optimizer':optimizer_meta,
    }
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(payload,indent=2))
    return payload
