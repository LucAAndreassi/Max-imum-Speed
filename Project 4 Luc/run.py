"""One command regenerates verified results, report tables, and figures."""
import csv
import hashlib
import json
import os
import platform
import copy
import shutil
from pathlib import Path
import sys
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).parent/".mpl-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import scipy
from scipy.integrate import solve_bvp
from arm import Beam,Elastica,section,hermite
from verify import run_verification

ROOT=Path(__file__).resolve().parent


def continuum_elastica(L,EI,P,q,follower=0.):
    """Independent reference BVP including integrated distributed dead load."""
    x=np.linspace(0,L,101)
    def ode(s,y,p):
        return np.vstack((y[1],-((P+q*(L-s))*np.cos(y[0])+follower*np.cos(p[0]-y[0]))/EI,
                          np.cos(y[0]),np.sin(y[0])))
    def bc(a,b,p): return np.array([a[0],b[1],a[2],a[3],b[0]-p[0]])
    total=P+follower
    guess=np.vstack((total*(L*x-x*x/2)/EI,total*(L-x)/EI,x,total*x*x*(3*L-x)/(6*EI)))
    sol=solve_bvp(ode,bc,x,guess,p=[total*L*L/(2*EI)],tol=1e-10,max_nodes=10000)
    if not sol.success: raise RuntimeError("Continuum elastica reference failed")
    return sol


def candidate(c,inputs,load_factors=None):
    L=inputs["length_m"]; E=inputs["material"]["E_bending_Pa"]
    g=inputs["gravity_m_s2"]; mt=inputs["tip_mass_kg"]
    D,t,mu=[c[k] for k in ["outer_diameter_m","wall_m","mass_per_length_kg_m"]]
    A,I=section(D,t); P=inputs["thrust_N"]-mt*g; q=-mu*g
    b=Beam(L,E,mu,mt,inputs["thrust_N"],g,inputs["elements"])
    u,r,_=b.solve(D,t)
    freq,modes=b.modes(D,t,4)
    sample_x=np.linspace(0,L,241)
    # Exact resultant polynomial includes possible interior stress maximum.
    check_x=[0.,L]
    if q!=0 and 0<L+P/q<L: check_x.append(L+P/q)
    moments=[P*(L-x)+q*(L-x)**2/2 for x in check_x]
    peak_M=max(abs(v) for v in moments)
    sigma=peak_M*(D/2)/I
    reference=continuum_elastica(L,E*I,-mt*g,q,inputs["thrust_N"])
    nl_ref=reference.sol(L)[3]
    refinement=[]
    for n in [40,80,160]:
        nl=Elastica(L,E*I,-mt*g,q,elements=n,follower_force=inputs["thrust_N"])
        theta,records=nl.solve(steps=1)
        z=nl.coordinates(theta)[-1,1]
        refinement.append(dict(elements=n,tip_z_m=float(z),relative_error_vs_BVP=float(abs(z/nl_ref-1))))
    curve=[]
    for factor in ([0.,.25,.5,.75,1.,2.,5.,10.] if load_factors is None else load_factors):
        if factor:
            obj=Elastica(L,E*I,-mt*g*factor,q*factor,elements=80,follower_force=inputs["thrust_N"]*factor)
            th,_=obj.solve(steps=5)
            zz=obj.coordinates(th)[-1,1]
        else: zz=0.
        curve.append(dict(load_factor=factor,linear_tip_z_m=float(u[-2]*factor),nonlinear_tip_z_m=float(zz)))
    shapes=[]; mode_shape=[]
    for x in sample_x:
        e=min(int(x/L*b.elements),b.elements-1);xi=x/L*b.elements-e
        shapes.append([float(x),0.,float(hermite(xi,L/b.elements)@u[2*e:2*e+4])])
        mode_shape.append(float(hermite(xi,L/b.elements)@modes[2*e:2*e+4,0]))
    mode_shape=np.array(mode_shape);mode_shape/=max(abs(mode_shape))
    exact_tip=P*L**3/(3*E*I)+q*L**4/(8*E*I)
    # This shear stress is a homogeneous-section screen, not a laminate allowable.
    R=D/2;ri=R-t
    shear=max(abs(P),abs(P+q*L))*(R**3-ri**3)/(3*I*t)
    Mfield=P*(L-sample_x)+q*(L-sample_x)**2/2
    typical=inputs["material"]["typical_tensile_strength_Pa"]
    return dict(**c,area_m2=float(A),I_m4=float(I),tube_mass_kg=mu*L,
                total_modeled_mass_kg=mt+mu*L,net_tip_force_N=P,distributed_load_N_m=q,
                root_moment_Nm=float(abs(r[1])),peak_bending_stress_Pa=float(sigma),
                homogeneous_shear_screen_Pa=float(shear),
                typical_tensile_ratio=float(typical/sigma),
                comparison_compression_ratio=float(inputs["material"]["comparison_only"]["compressive_strength_Pa"]/sigma),
                required_axial_strength_for_sf2_Pa=float(inputs["required_safety_factor"]*sigma),
                required_homogeneous_shear_strength_for_sf2_Pa=float(inputs["required_safety_factor"]*shear),
                overall_safety_factor=None,safety_status="UNRESOLVED: compression, shear, local buckling, joint data absent",
                tip_deflection_m=float(u[-2]),exact_linear_tip_m=exact_tip,
                linear_tip_relative_error=float(abs(u[-2]/exact_tip-1)),
                first_bending_frequency_Hz=float(freq[0]),bending_frequencies_Hz=freq.tolist(),
                vibration=b.vibration_summary(D,t,4),
                nonlinear_tip_m=float(nl_ref),nonlinear_FEM_tip_m=float(refinement[-1]["tip_z_m"]),
                nonlinear_change_percent=float(100*(nl_ref/u[-2]-1)),
                nonlinear_refinement=refinement,load_path=curve,
                linear_shape=shapes,normalized_mode=mode_shape.tolist(),
                nonlinear_shape=[[float(x),0.,float(z)] for x,z in nl.coordinates_at(theta,sample_x)],
                nonlinear_reference_shape=[[float(x),0.,float(z)] for x,z in reference.sol(sample_x)[[2,3]].T],
                sample_bending_stress_Pa=(abs(Mfield)*R/I).tolist(),
                raw_displacements=u.tolist(),
                checks=dict(strength="UNRESOLVED",deflection="LIMIT NOT SPECIFIED",frequency="LIMIT NOT SPECIFIED",
                            nonlinear_design_comparison="PENDING OPTIMIZATION"))


def plots(data):
    out=ROOT/"assets";out.mkdir(exist_ok=True)
    v=data["verification"]
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})
    fig,ax=plt.subplots(1,2,figsize=(10,3.5))
    a=v["displacement_convergence"]
    ax[0].loglog([z["h"] for z in a],[z["error"] for z in a],"o-",color="#176678")
    ax[0].set(xlabel="Element length h (m)",ylabel="L2 displacement error",title="Synthetic uniform-load beam")
    a=v["frequency_convergence"]
    ax[1].loglog([z["h"] for z in a],[z["error"] for z in a],"o-",color="#176678")
    ax[1].set(xlabel="Element length h (m)",ylabel="Relative first-frequency error",title="Synthetic cantilever frequency")
    for aa in ax:
        aa.grid(alpha=.2);aa.set_xticks([.0625,.125,.25,.5],labels=['1/16','1/8','1/4','1/2'])
        aa.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout();fig.savefig(out/"mesh-convergence.png");plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(12,3.5))
    for aa,response in zip(ax,["tip_displacement","root_stress","bending_frequency"]):
        for d in range(5):
            a=[z for z in v["gradient_sweep"] if z["response"]==response and z["direction"]==d]
            aa.loglog([z["h"] for z in a],[max(z["error"],1e-16) for z in a],".-",label=f"direction {d}")
        aa.set(title=response.replace("_"," "),xlabel="Scaled central-difference step",ylabel="Relative directional error")
        aa.grid(alpha=.2)
    ax[0].legend(fontsize=7);fig.tight_layout();fig.savefig(out/"gradient-check.png");plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,3.5))
    a=v["nonlinear_gradient_sweep"]
    ax.loglog([z["h"] for z in a],[z["error"] for z in a],'o-',color='#176678')
    ax.set(xlabel='Central-difference step in EI (synthetic units)',ylabel='Relative directional error',
           title='Follower-thrust nonlinear adjoint')
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/"nonlinear-gradient.png");plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,3.6))
    for c in data["candidates"]:
        a=c["load_path"]
        ax.plot([z["load_factor"] for z in a],[1000*z["nonlinear_tip_z_m"] for z in a],"o-",label=c["name"]+" elastica")
        ax.plot([z["load_factor"] for z in a],[1000*z["linear_tip_z_m"] for z in a],"--",color=ax.lines[-1].get_color())
    ax.set(xlabel="Multiplier on net tip force and tube weight",ylabel="Upward tip deflection (mm)",title="Exploratory load path; dashed = linear")
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/"nonlinear-path.png");plt.close(fig)


def optimization_report(data):
    opt=data["optimization"]
    if opt.get("status")!="converged_preliminary_model":
        return "No current optimizer result. Run `python optimize.py --update-report`. Status: "+opt.get("status","missing"), "Optimization has not been regenerated for the current inputs."
    initial=opt["initial_design"];linear=opt["without_nonlinear"]["optimum"];best=opt["with_nonlinear"]["optimum"]
    rounded=opt.get("rounded_design")
    designs=[initial,linear,best]+([rounded] if rounded else [])
    names=["Starting design","Linear constraints","Linear + follower constraints"]+(["Rounded candidate"] if rounded else [])
    table="| Quantity | "+" | ".join(names)+" |\n|---|"+"---:|"*len(names)+"\n"
    fields=[("OD (mm)","outer_diameter_mm",".5f"),("Wall (mm)","wall_mm",".5f"),
            ("Tube mass (g)","tube_mass_g",".4f"),("Minimum modeled SF","minimum_modeled_safety_factor",".6f"),
            ("Linear tip deflection (mm)","linear_tip_deflection_mm",".3f"),
            ("Nonlinear FE tip deflection (mm)","nonlinear_tip_deflection_mm",".3f"),
            ("First unloaded bending frequency (Hz)","first_bending_frequency_Hz",".3f")]
    for label,key,fmt in fields:
        table+="| "+label+" | "+" | ".join("Not evaluated" if d[key] is None else format(d[key],fmt) for d in designs)+" |\n"
    headline=(f"The preliminary optimizer selects a **{best['outer_diameter_mm']:.3f} mm OD, {best['wall_mm']:.2f} mm wall** uniform tube "
              f"weighing **{best['tube_mass_g']:.2f} g**. "
              f"Its minimum modeled nominal strength SF is **{best['minimum_modeled_safety_factor']:.5f}**. "
              f"It has **{best['linear_tip_deflection_mm']:.1f} mm linear tip deflection** and a **{best['first_bending_frequency_Hz']:.2f} Hz** unloaded bending frequency.")
    fmin=opt.get('resolved_minimum_bending_frequency_Hz')
    if fmin is not None:
        headline+=f" The configured minimum bending frequency is **{fmin:.2f} Hz**; this does not establish operating dynamic stability."
    text="The constrained sizing search in `optimize_arm.py` uses the existing `optimize.TubeModel`, `arm.Beam` and `arm.Elastica` FEM and adjoints. Edit `optimizer_inputs.json` and run `python optimize_arm.py --update-report`; `optimize.py` remains a compatible entry point.\n\n"
    text+="The user now requires **wall thickness strictly greater than 2 mm**. The numerical lower bound is 2.01 mm; a 0.01 mm difference is not a manufacturing allowance. The OD interval is 10-60 mm and the wall upper bound is 3 mm, both explicit provisional study bounds. The reference 40 mm OD / 2.5 mm wall tube is an analyst-selected starting design under the same mass/material model. The original 1.5 mm stock candidates are excluded by the updated wall requirement.\n\n"
    text+="The effective material profile uses $E=64$ GPa, $\\rho=1520$ kg/m3, $X_t=650$ MPa, $X_c=910.1$ MPa and $S=103.4$ MPa. The modulus/tensile value comes from CFT-WF typical data; density/compression are comparison estimates. The torsional strength is used as an **unverified longitudinal-shear proxy**, not an interlaminar allowable. These are explicit cross-product estimates authorized for this preliminary study, not measured properties of a new laminate. No interpolation endpoints were invented: the default profile is constant; an optional property table supports piecewise-linear thickness interpolation inside its data range only.\n\n"
    text+="[CFT-WF source](https://media.easycomposites.co.uk/datasheets/EC-TDS-Carbon-Fibre-Woven-Finish-Tube.pdf), [Graphtek comparison source](https://www.graphitestore.com/core/media/media.nl?c=4343521&h=22DqRe5XcXtOFJUqMtPHOuVeC4o8PE9P-gnJuKBYWPrJ3COD&id=9929).\n\n"
    text+="The objective is tube mass $m=\\rho A L$. The configured tip assembly mass is constant within each sizing run, so minimizing total modeled mass gives the same geometry. SLSQP minimizes mass divided by reference mass, with design variables scaled by reference OD and wall. At 21 fixed arc-length stations, both outer-fiber normal stresses $N/A\\pm MD/(2I)$ and nominal section shear are evaluated. Each failure family uses the smooth upper bound $\\mathrm{KS}_{80}(2\\sigma/X)-1\\le0$ (and the corresponding shear ratio). Thus a feasible aggregate implies every sampled nominal stress has SF at least 2. Linear stresses recover the exact uniform-load moment field from FE reactions; nonlinear resultants come from equilibrium on the deformed FE centerline, including tip-following thrust and vertical gravity.\n\n"
    text+="The study with N2 retains the three linear strength constraints and adds three follower-load strength constraints. A minimum bore of 1 mm keeps a hollow section. The reported SF is the minimum unaggregated reserve factor across linear/nonlinear tension, compression and shear.\n\n"
    if fmin is not None:
        text+=f"The supplied script's frequency requirement is retained: the configured minimum is **{fmin:.6f} Hz**. The default rotor multiple is 1, so this is recalculated from the operating RPM divided by 60. It is a preliminary sizing criterion; rotor harmonics, gyroscopic coupling and flutter are not assessed.\n\n"
    else:
        text+="No minimum frequency limit is enabled for this study.\n\n"
    if opt['study'].get('deflection_limit_mm') is None:
        text+="No deflection limit is configured. The supplied file's stiffness/travel defaults were in missing dependency modules and have not been guessed.\n\n"
    if not initial['feasible']:
        text+="The starting 40 mm / 2.5 mm tube fails an enabled constraint; use the feasible fixed-wall reference below for a like-for-like mass saving.\n\n"
    text+=table+"\n"
    text+="The rounded geometry is a checked custom-tube candidate, not a verified stock product or a manufacturing tolerance analysis. Its nominal reserve must still cover any real geometric tolerance and property knockdown.\n\n"
    text+="### Optimizer verification and convergence\n\n"
    text+=f"All **{len(opt['gradient_checks'])}** mass/constraint partial-and-adjoint checks passed against central differences, including the density-dependent self-weight. Optimizer unit tests separately check the closed-form static solution, optional deflection/frequency gradients, tabulated properties, the follower tangent, an independent exact strength-only lower bound, infeasibility handling, the strict wall rule, rejection of nonstationary solver successes, recovery from failed starts, frequency-reference crossings, and shadow prices against reoptimization.\n\n"
    text+="| Response | Variable | Adjoint/direct analytic | Central difference | Relative error |\n|---|---|---:|---:|---:|\n"
    for r in opt["gradient_checks"]:
        text+=f"| {r['response']} | {r['variable']} | {r['adjoint']:.8g} | {r['finite_difference']:.8g} | {r['relative_error']:.3g} |\n"
    text+="\n| Case | Start | Solver success | Accepted (incl. KKT) | Tube mass (g) | Solver message |\n|---|---:|---|---|---:|---|\n"
    for case,key in [("Linear","without_nonlinear"),("Linear + follower","with_nonlinear")]:
        for r in opt[key]["starts"]:
            mass_text=f"{r['design']['tube_mass_g']:.5f}" if r['design'] else "Unavailable"
            text+=f"| {case} | {r['start']+1} | {r['success']} | {r['accepted']} | {mass_text} | {r['message']} |\n"
    text+="\nConfigured starts, both bounding-box corners, and a feasible fixed-wall reference are searched, with duplicate starts removed. Only solver-successful, feasible runs passing primal, dual, stationarity and complementarity tolerances are eligible. A failed start is recorded and the remaining starts continue.\n\n![Optimization convergence](assets/optimization-convergence.png)\n\n"
    reference=opt['with_nonlinear'].get('fixed_wall_reference')
    if reference and reference['point']:
        p=reference['point']
        text+=f"The fixed-wall nominal reference uses **{p['outer_diameter_mm']:.3f} mm OD / {p['wall_mm']:.3f} mm wall**, at **{p['mass_kg']*1000:.3f} g**. Its OD is screened on a {reference['grid_mm']:g} mm grid and it passes all enabled constraints.\n\n"
        savings=opt.get('fixed_wall_mass_reduction_percent')
        if savings is not None:
            text+=f"Optimizing both dimensions saves **{savings:.3f}%** tube mass relative to this feasible fixed-wall reference.\n\n"
    k=opt["with_nonlinear"]["kkt"]
    text+="### KKT interpretation\n\n| Active requirement | Scaled multiplier | Residual |\n|---|---:|---:|\n"
    for a in k["active"]:
        text+=f"| {a['name']} | {a['scaled_multiplier']:.8g} | {a['residual']:.3g} |\n"
    text+=f"\nStationarity infinity norm: **{k['stationarity_inf_norm']:.3g}**; complementarity infinity norm: **{k['complementarity_inf_norm']:.3g}**; maximum violation: **{k['maximum_violation']:.3g}**. Multipliers are reconstructed by nonnegative least squares using mass/reference-mass and dimensions/reference-dimensions. The table converts bound multipliers back to the mm convention; JSON also records kg/mm, active-set rank and all four acceptance checks. This is local optimality evidence, not proof of a global optimum.\n\n"
    active={a["name"]:a["scaled_multiplier"] for a in k["active"]}
    if "wall_lower_bound" in active:
        sensitivity=initial["tube_mass_g"]*active["wall_lower_bound"]
        text+=f"The wall lower bound is active: increasing it by 0.1 mm has a local predicted mass cost of approximately **{sensitivity*.1:.3f} g**, while the active set remains unchanged.\n\n"
    shadows=opt['with_nonlinear'].get('shadow_checks',[])
    if shadows:
        text+="| Tightened active requirement | Step | Predicted mass change (g) | Reoptimized mass change (g) |\n|---|---:|---:|---:|\n"
        for c in shadows:
            measured="Unavailable" if c['change_kg'] is None else f"{1000*c['change_kg']:.6f}"
            text+=f"| {c['name']} | {c['step']:g} {c['unit']} | {1000*c['predicted_kg']:.6f} | {measured} |\n"
        text+="\nThese perturbations check local multiplier predictions; their sizes are numerical steps, not manufacturing tolerances.\n\n"
    text+="### Nonlinear and mesh comparison\n\n"
    text+=f"Adding follower-load constraints changes the optimized mass by **{opt['nonlinear_mass_change_percent']:.3g}%**. The selected active requirements are **{', '.join(active)}**. This comparison applies to the configured load, frequency and strength limits; it does not establish dynamic stability.\n\n"
    text+="| Linear elements | Nonlinear elements | Minimum modeled SF | Nonlinear tip (mm) | First frequency (Hz) | Max constraint residual |\n|---:|---:|---:|---:|---:|---:|\n"
    for r in opt["mesh_refinement"]:
        text+=f"| {r['linear_elements']} | {r['nonlinear_elements']} | {r['minimum_modeled_safety_factor']:.8f} | {r['nonlinear_tip_deflection_mm']:.6f} | {r['first_bending_frequency_Hz']:.6f} | {r['max_constraint']:.3g} |\n"
    text+="\nThe continuous optimum passes the refined nominal strength checks. Section shear remains a homogeneous-section screen. Shell buckling, joints, laminate interaction, damage and dynamic follower-force stability are unassessed. Therefore the **overall structural SF remains unresolved**, even though the specified nominal material-strength model meets SF 2.\n"
    return text,headline


def render_report(data):
    c=data["candidates"];v=data["verification"]
    table="| Quantity | "+" | ".join(x["name"] for x in c)+" |\n|---|---:|---:|\n"
    fields=[("Tube mass (g)","tube_mass_kg",1000,".2f"),
            ("Tube + coaxial assembly mass (g)","total_modeled_mass_kg",1000,".2f"),
            ("Peak nominal bending stress (MPa)","peak_bending_stress_Pa",1e-6,".3f"),
            ("Homogeneous-section shear screen (MPa)","homogeneous_shear_screen_Pa",1e-6,".3f"),
            ("Typical tensile strength / bending stress (NOT overall SF)","typical_tensile_ratio",1,".2f"),
            ("Comparable-product compression / bending stress (estimate only)","comparison_compression_ratio",1,".2f"),
            ("Minimum axial tensile AND compressive strength needed for SF 2 (MPa)","required_axial_strength_for_sf2_Pa",1e-6,".3f"),
            ("Minimum homogeneous shear strength needed for SF 2 (MPa)","required_homogeneous_shear_strength_for_sf2_Pa",1e-6,".3f"),
            ("Linear upward tip deflection (mm)","tip_deflection_m",1000,".4f"),
            ("First bending frequency with tip mass (Hz)","first_bending_frequency_Hz",1,".3f"),
            ("Nonlinear continuum reference tip deflection (mm)","nonlinear_tip_m",1000,".4f"),
            ("Nonlinear change from linear (%)","nonlinear_change_percent",1,".5f")]
    for label,key,scale,fmt in fields:
        table+="| "+label+" | "+" | ".join(format(x[key]*scale,fmt) for x in c)+" |\n"
    table+="| Overall SF >= 2 | UNRESOLVED | UNRESOLVED |\n| Deflection / frequency acceptance | Limits pending | Limits pending |\n"
    ladder="| Check | Unit | Computed | Reference | Error | Acceptance |\n|---|---|---:|---:|---:|---|\n"
    for x in v["rows"]:
        ladder+=f"| {x['test']} | {x['unit']} | {x['computed']:.8g} | {x['expected']:.8g} | {x['error']:.3g} ({x['error_kind']}) | {'PASS' if x['passed'] else 'FAIL'}; <= {x['tolerance']:.3g} |\n"
    refinements="| Tube | Elements | Nonlinear FE tip (mm) | Relative difference from continuum BVP |\n|---|---:|---:|---:|\n"
    for a in c:
        for z in a["nonlinear_refinement"]:
            refinements+=f"| {a['name']} | {z['elements']} | {1000*z['tip_z_m']:.6f} | {z['relative_error_vs_BVP']:.3g} |\n"
    opt_text,opt_headline=optimization_report(data)
    replacements={"@@RESULTS@@":table,"@@VERIFICATION@@":ladder,"@@REFINEMENT@@":refinements,
                  "@@OPTIMIZATION@@":opt_text,"@@OPTIMIZATION_HEADLINE@@":opt_headline,
                  "@@COUNTS@@":f"{sum(x['passed'] for x in v['rows'])}/{len(v['rows'])}",
                  "@@ADJOINT_MS@@":f"{v['timing']['adjoint_tip_gradient_ms']:.3f}",
                  "@@FD_MS@@":f"{v['timing']['central_difference_two_variables_ms']:.3f}",
                  "@@ROOT_MOMENT@@":f"{data['known_loads']['thrust_only_root_moment_Nm']:.6f}",
                  "@@TIP_MASS_LB@@":f"{data['inputs']['tip_mass_kg']/.45359237:g}",
                  "@@TIP_MASS_KG@@":f"{data['inputs']['tip_mass_kg']:.9f}",
                  "@@SHAFT@@":f"{data['known_loads']['shaft_frequency_Hz']:.6f}",
                  "@@BPF@@":f"{data['known_loads']['blade_passing_frequency_Hz']:.6f}",
                  "@@MASS_CHANGE@@":f"{100*(1-c[1]['tube_mass_kg']/c[0]['tube_mass_kg']):.2f}"}
    report=(ROOT/"report_template.md").read_text(encoding="utf-8")
    for a,b in replacements.items(): report=report.replace(a,b)
    if "@@" in report: raise RuntimeError("Unexpanded report token")
    (ROOT/"index.md").write_text(report,encoding="utf-8")
    # Plain-browser preview shares exact generated Markdown content with Pages.
    (ROOT/"assets"/"report.md").write_text(report.split("---",2)[2].strip(),encoding="utf-8")


def main():
    os.chdir(ROOT)
    inp=json.loads((ROOT/"inputs.json").read_text())
    if inp["applied_torque_Nm"]!=0: raise ValueError("Reduced solver requires zero applied torque")
    v=run_verification()
    if v["status"]!="pass":
        print(json.dumps([r for r in v["rows"] if not r["passed"]],indent=2))
        raise RuntimeError("Verification failed; report not regenerated")
    data=dict(schema_version=1,status="PRELIMINARY - NO VERIFIED SAFE OPTIMUM",inputs=inp,
              known_loads=dict(shaft_frequency_Hz=inp["rpm"]/60,
                               blade_passing_frequency_Hz=inp["rpm"]*inp["blades_per_rotor"]/60,
                               thrust_only_root_moment_Nm=inp["thrust_N"]*inp["length_m"]),
              provenance=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                              source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                              for p in ["inputs.json","arm.py","vibration.py","verify.py","run.py","report_template.md"]}),
              verification=v,candidates=[candidate(c,inp) for c in inp["candidates"]],
              optimization=dict(status="blocked_by_missing_design_inputs",optimized_design=None,
                                reason="Compression/shear/local buckling/joint allowables and stiffness limits absent; stock products are discrete; uniform two-variable scope needs confirmation."))
    opt_path=ROOT/"output"/"optimization"/"results.json"
    if opt_path.exists():
        from optimize import fingerprint
        opt=json.loads(opt_path.read_text())
        study=json.loads((ROOT/"optimizer_inputs.json").read_text())
        if opt.get("provenance")==fingerprint(inp,study) and opt.get("status")=="converged_preliminary_model":
            data["optimization"]=opt
            data["status"]="PRELIMINARY OPTIMUM - MODELED SF >= 2 WITH ESTIMATED MATERIAL"
            viewer=[]
            for name,summary in [("Starting tube - estimated material",opt["initial_design"]),
                                 ("Optimized tube - estimated material",opt["with_nonlinear"]["optimum"])]:
                props=summary["material_properties"];view_inputs=copy.deepcopy(inp)
                view_inputs["material"]["E_bending_Pa"]=props["E_Pa"]
                view_inputs["material"]["typical_tensile_strength_Pa"]=props["tensile_strength_Pa"]
                D,t=summary["outer_diameter_mm"]/1000,summary["wall_mm"]/1000
                A,_=section(D,t)
                c=dict(name=name,outer_diameter_m=D,wall_m=t,mass_per_length_kg_m=props["density_kg_m3"]*A,
                       source="optimizer_inputs.json estimated profile")
                vdesign=candidate(c,view_inputs,load_factors=[0.,.25,.5,.75,1.])
                vdesign["optimization_summary"]=summary
                viewer.append(vdesign)
            data["optimization"]["viewer_designs"]=viewer
            shutil.copyfile(ROOT/"output"/"optimization"/"convergence.png",ROOT/"assets"/"optimization-convergence.png")
        else:
            data["optimization"]={"status":"stale_result","optimized_design":None,
                                  "reason":"Run python optimize.py --update-report for the current sources and inputs"}
    (ROOT/"assets").mkdir(exist_ok=True)
    (ROOT/"assets"/"results.json").write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    with (ROOT/"assets"/"verification.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(v["rows"][0]));writer.writeheader();writer.writerows(v["rows"])
    plots(data);render_report(data)
    print(f"Verification: {len(v['rows'])} checks passed. Regenerated index.md, assets/results.json, tables, figures.")
    for c in data["candidates"]:
        print(f"{c['name']}: tube {c['tube_mass_kg']*1000:.2f} g; stress {c['peak_bending_stress_Pa']/1e6:.3f} MPa; tip {c['tip_deflection_m']*1000:.3f} mm; f1 {c['first_bending_frequency_Hz']:.3f} Hz; SF unresolved")


if __name__=="__main__":main()
