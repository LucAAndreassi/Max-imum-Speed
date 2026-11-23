r"""Uniform-tube mass minimization using own FE solves and implicit gradients.

Run: .venv\Scripts\python.exe optimize.py --update-report
All SF statements apply only to the explicitly recorded effective-material model.
"""
import argparse
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import eigh

from arm import Beam, Elastica, section
from vibration import first_frequency_partials, modal_summary, tip_mass_frequencies, METHOD as VIBRATION_METHOD

ROOT = Path(__file__).resolve().parent
PROPERTY_KEYS = ("E_Pa", "density_kg_m3", "tensile_strength_Pa",
                 "compressive_strength_Pa", "shear_strength_Pa")


class Material:
    """Complex-safe constants or genuine tabulated thickness interpolation."""
    def __init__(self, config):
        self.config = copy.deepcopy(config)
        self.rows = config.get("thickness_samples", [])
        if self.rows:
            self.knots = np.array([r["wall_mm"] for r in self.rows], float)
            if len(self.rows) < 2 or not np.all(np.isfinite(self.knots)) or np.any(np.diff(self.knots) <= 0):
                raise ValueError("Material thickness samples need at least two strictly increasing finite wall_mm values")
            values = [r[k] for r in self.rows for k in PROPERTY_KEYS]
        else:
            values = [config[k] for k in PROPERTY_KEYS]
        if any(v is None or not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError("All material properties must be finite and positive")

    def at(self, wall_mm):
        if not self.rows:
            return {k: self.config[k] for k in PROPERTY_KEYS}
        t = float(np.real(wall_mm))
        if t < self.knots[0]-1e-12 or t > self.knots[-1]+1e-12:
            raise ValueError("Wall thickness is outside the material table: extrapolation is disabled")
        j = min(max(np.searchsorted(self.knots, t, side="right")-1, 0), len(self.knots)-2)
        a, b = self.rows[j], self.rows[j+1]
        fraction = (wall_mm-self.knots[j])/(self.knots[j+1]-self.knots[j])
        return {k: a[k]+fraction*(b[k]-a[k]) for k in PROPERTY_KEYS}


def smooth_max(values, rho):
    shift = float(np.max(np.real(values)))
    return shift + np.log(np.sum(np.exp(rho*(values-shift))))/rho


def minimum_frequency(loads,study):
    """Resolve either a fixed Hz limit or the supplied optimizer's rotor margin."""
    fixed=study.get("minimum_bending_frequency_Hz")
    multiple=study.get("minimum_bending_frequency_rotor_multiple")
    if multiple is None:
        return fixed
    if fixed is not None:
        raise ValueError("Specify either a fixed frequency or a rotor multiple, not both")
    if not np.isfinite(multiple) or multiple<=0 or not np.isfinite(loads["rpm"]) or loads["rpm"]<=0:
        raise ValueError("Rotor frequency multiple and operating RPM must be positive and finite")
    return multiple*loads["rpm"]/60


class TubeModel:
    """d=[OD_mm,wall_mm]. All physical equations and FE matrices use SI."""
    def __init__(self, loads, study, nonlinear=True):
        self.loads, self.study = copy.deepcopy(loads), copy.deepcopy(study)
        self.nonlinear = nonlinear
        self.material = Material(study["material"])
        self.L, self.g = loads["length_m"], loads["gravity_m_s2"]
        self.T, self.mt = loads["thrust_N"], loads["tip_mass_kg"]
        self.sf, self.rho = study["required_safety_factor"], study["ks_parameter"]
        self.n, self.nl_n = study["linear_elements"], study["nonlinear_elements"]
        self.minimum_frequency = minimum_frequency(loads,study)
        self.stations = np.linspace(0,self.L,study["evaluation_stations"])
        self.names = ["linear_tension", "linear_compression", "linear_shear"]
        if nonlinear:
            self.names += ["nonlinear_tension", "nonlinear_compression", "nonlinear_shear"]
        if study.get("deflection_limit_mm") is not None:
            self.names += ["linear_deflection"]
            if nonlinear: self.names += ["nonlinear_deflection"]
        if self.minimum_frequency is not None:
            self.names += ["bending_frequency"]
        self.names += ["inner_diameter"]
        self._key = None
        self._value = None

    def properties(self, d):
        p = self.material.at(d[1])
        D, t = d/1000
        A, I = section(D,t)
        return p, D, t, A, I, p["density_kg_m3"]*A

    def linear_parts(self,d):
        p,D,t,A,I,mu = self.properties(d)
        beam = Beam(self.L,p["E_Pa"],mu,self.mt,self.T,self.g,self.n)
        return beam, beam.matrices(D,t)

    def nonlinear_model(self,d,n=None):
        p,D,t,A,I,mu = self.properties(d)
        return Elastica(self.L,p["E_Pa"]*I,-self.mt*self.g,-mu*self.g,
                        elements=self.nl_n if n is None else n,follower_force=self.T)

    def linear_fields(self,d,u):
        p,D,t,A,I,mu = self.properties(d)
        _,(K,_,F,_) = self.linear_parts(d)
        reactions = K@u-F
        s = self.stations
        q = -mu*self.g
        moment = -reactions[1]+reactions[0]*s+q*s*s/2
        shear = reactions[0]+q*s
        axial = np.zeros_like(moment)
        return axial, moment, shear, u[-2]

    def nonlinear_fields(self,d,theta):
        p,D,t,A,I,mu = self.properties(d)
        n = len(theta)-1
        xi,w = np.polynomial.legendre.leggauss(4)
        xi,w = (xi+1)/2,w/2
        l = self.L/n
        angle = theta[:-1,None]*(1-xi)+theta[1:,None]*xi
        sgauss = (np.arange(n)[:,None]+xi)*l
        dx = l*(np.cos(angle)@w)
        dz = l*(np.sin(angle)@w)
        x = np.r_[0.,np.cumsum(dx)]
        z = np.r_[0.,np.cumsum(dz)]
        # Integral_s^L (x(v)-x(s)) dv = Integral_s^L (L-v) cos(theta(v)) dv.
        gravity_integral = l*((self.L-sgauss)*np.cos(angle))@w
        tail = np.r_[np.cumsum(gravity_integral[::-1])[::-1],0.]
        Fx = -self.T*np.sin(theta[-1])
        Fz = self.T*np.cos(theta[-1])-self.mt*self.g
        s = np.linspace(0,self.L,n+1)
        q = -mu*self.g
        moment = (x[-1]-x)*Fz-(z[-1]-z)*Fx+q*tail
        vertical = Fz+q*(self.L-s)
        axial = Fx*np.cos(theta)+vertical*np.sin(theta)
        shear = -Fx*np.sin(theta)+vertical*np.cos(theta)
        ids = np.rint(np.linspace(0,n,len(self.stations))).astype(int)
        return axial[ids],moment[ids],shear[ids],z[-1]

    def utilizations(self,d,fields):
        p,D,t,A,I,mu = self.properties(d)
        axial,moment,shear,_ = fields
        stresses = np.r_[axial/A+moment*D/(2*I),axial/A-moment*D/(2*I)]
        R,ri = D/2,D/2-t
        # Smooth force magnitude only removes the derivative singularity at zero.
        tau = np.sqrt(shear**2+1e-20)*(R**3-ri**3)/(3*I*t)
        true = [self.sf*stresses/p["tensile_strength_Pa"],
                -self.sf*stresses/p["compressive_strength_Pa"],
                self.sf*tau/p["shear_strength_Pa"]]
        return true

    def state_outputs(self,d,state,nonlinear):
        fields = self.nonlinear_fields(d,state) if nonlinear else self.linear_fields(d,state)
        out = [smooth_max(v,self.rho)-1 for v in self.utilizations(d,fields)]
        if self.study.get("deflection_limit_mm") is not None:
            out.append(np.sqrt(fields[3]**2+1e-24)/(self.study["deflection_limit_mm"]/1000)-1)
        return np.array(out)

    def mass(self,d):
        return self.properties(d)[-1]*self.L

    def state_gradient(self,d,state,nonlinear,parts):
        """Partial complex step + an actual transposed tangent/stiffness solve."""
        h = 1e-25
        active = np.arange(1,len(state)) if nonlinear else parts[3]
        values = self.state_outputs(d,state,nonlinear)
        qu = np.zeros((len(values),len(active)))
        for j,i in enumerate(active):
            sc = state.astype(complex); sc[i] += 1j*h
            qu[:,j] = self.state_outputs(d,sc,nonlinear).imag/h
        tangent = parts if nonlinear else parts[0][active,:][:,active].toarray()
        if nonlinear:
            tangent = tangent[1:,1:]
        adjoint = np.linalg.solve(tangent.T,qu.T)
        grad = np.zeros((len(values),2))
        for j in range(2):
            dc = d.astype(complex); dc[j] += 1j*h
            qd = self.state_outputs(dc,state,nonlinear).imag/h
            if nonlinear:
                rd = self.nonlinear_model(dc).residual(state)[0].imag[active]/h
            else:
                _,(K,_,F,_) = self.linear_parts(dc)
                rd = (K.imag@state-F.imag)[active]/h
            grad[:,j] = qd-adjoint.T@rd
        return grad

    def evaluate(self,d,gradients=True):
        d = np.asarray(d,float)
        key = (tuple(d),gradients)
        if self._key == key: return self._value
        beam,(K,M,F,free) = self.linear_parts(d)
        p,D,t,A,I,mu = self.properties(d)
        u,reactions,parts = beam.solve(D,t)
        values = self.state_outputs(d,u,False)
        jac = self.state_gradient(d,u,False,(K,M,F,free)) if gradients else None
        constraints = list(values[:3]); rows = list(jac[:3]) if gradients else []
        fields = {"linear":self.linear_fields(d,u)}
        theta = None
        if self.nonlinear:
            nl = self.nonlinear_model(d)
            theta,history = nl.solve(steps=4)
            nv = self.state_outputs(d,theta,True)
            nj = self.state_gradient(d,theta,True,nl.residual(theta)[1]) if gradients else None
            constraints.extend(nv[:3])
            if gradients: rows.extend(nj[:3])
            fields["nonlinear"] = self.nonlinear_fields(d,theta)
        if self.study.get("deflection_limit_mm") is not None:
            constraints.append(values[3])
            if gradients: rows.append(jac[3])
            if self.nonlinear:
                constraints.append(nv[3])
                if gradients: rows.append(nj[3])
        f = tip_mass_frequencies(self.L,p["E_Pa"]*I,mu,self.mt,1)[0]
        if self.minimum_frequency is not None:
            fmin = self.minimum_frequency
            constraints.append(fmin/f-1)
            if gradients:
                df=[]
                _,fp=first_frequency_partials(self.L,p["E_Pa"]*I,mu,self.mt)
                for j in range(2):
                    dc=d.astype(complex);dc[j]+=1e-25j
                    pc,_,_,_,ic,muc=self.properties(dc)
                    derivative=fp["EI"]*(pc["E_Pa"]*ic).imag/1e-25+fp["mass_per_length"]*muc.imag/1e-25
                    df.append(-fmin/f**2*derivative)
                rows.append(df)
        constraints.append((self.study["minimum_inner_diameter_mm"]-(d[0]-2*d[1]))/10)
        if gradients: rows.append([-.1,.2])
        mass_grad = np.zeros(2)
        if gradients:
            for j in range(2):
                dc=d.astype(complex);dc[j]+=1e-25j
                mass_grad[j]=self.mass(dc).imag/1e-25
        result=dict(d=d,mass=float(mu*self.L),g=np.array(constraints,float),
                    jac=np.array(rows,float) if gradients else None,mass_grad=mass_grad,
                    u=u,theta=theta,fields=fields,frequency=float(f))
        self._key,self._value=key,result
        return result

    def summary(self,d):
        a=self.evaluate(d,False)
        p,D,t,A,I,mu=self.properties(np.array(d))
        checks={}
        reserve=[]
        for state,fields in a["fields"].items():
            util=self.utilizations(np.array(d),fields)
            checks[state]={}
            for kind,v in zip(["tension","compression","shear"],util):
                demand=max(float(np.max(np.real(v))),0.)
                sf=self.sf/demand if demand else None
                checks[state][kind]={"safety_factor":sf,"status":"PASS" if sf is None or sf>=self.sf-1e-6 else "FAIL"}
                if sf is not None:reserve.append(sf)
        return dict(outer_diameter_mm=float(d[0]),wall_mm=float(d[1]),inner_diameter_mm=float(d[0]-2*d[1]),
                    tube_mass_g=a["mass"]*1000,total_modeled_mass_g=(a["mass"]+self.mt)*1000,
                    minimum_modeled_safety_factor=min(reserve),strength_checks=checks,
                    linear_tip_deflection_mm=float(a["fields"]["linear"][3]*1000),
                    nonlinear_tip_deflection_mm=float(a["fields"]["nonlinear"][3]*1000) if self.nonlinear else None,
                    first_bending_frequency_Hz=a["frequency"],
                    vibration_method=VIBRATION_METHOD,
                    vibration=modal_summary(self.L,p["E_Pa"]*I,mu,self.mt,1),
                    constraint_residuals=dict(zip(self.names,a["g"].tolist())),
                    feasible=bool(np.max(a["g"])<=2e-6),
                    material_properties=p,overall_structure_safety_factor=None,
                    overall_structure_status="UNRESOLVED: shell buckling, laminate interaction, joints and dynamic stability excluded")


def validate(loads,study):
    if loads.get("applied_torque_Nm",0)!=0:raise ValueError("This model requires zero applied torque")
    for key in ["outer_diameter_bounds_mm","wall_bounds_mm"]:
        b=study[key]
        if len(b)!=2 or not all(np.isfinite(v) and v>0 for v in b) or b[1]<=b[0]:
            raise ValueError(f"Invalid {key}")
    if study["wall_bounds_mm"][0]<=2:
        raise ValueError("The user requires wall thickness strictly greater than 2 mm")
    if study["outer_diameter_bounds_mm"][0] <= 2*study["wall_bounds_mm"][1]:
        raise ValueError("Bounds must keep every evaluated section hollow (OD_min > 2*wall_max)")
    if study["required_safety_factor"]<2:raise ValueError("Required safety factor must be at least 2")
    if study["linear_elements"]<2 or study["nonlinear_elements"]<20:
        raise ValueError("Use at least 2 linear and 20 nonlinear elements")
    if study["evaluation_stations"]<3 or study["nonlinear_elements"]%(study["evaluation_stations"]-1):
        raise ValueError("Nonlinear element count must be a multiple of evaluation_stations - 1")
    for key in ["deflection_limit_mm","minimum_bending_frequency_Hz"]:
        if study.get(key) is not None and (not np.isfinite(study[key]) or study[key]<=0):
            raise ValueError(f"{key} must be positive or null")
    minimum_frequency(loads,study)
    if study["ks_parameter"]<=0:raise ValueError("KS parameter must be positive")
    material=Material(study["material"])
    material.at(study["wall_bounds_mm"][0]);material.at(study["wall_bounds_mm"][1])


def gradient_check(model,d):
    d=np.asarray(d,float)
    a=model.evaluate(d)
    rows=[]
    for j in range(2):
        h=1e-4*max(abs(d[j]),1.)
        delta=np.eye(2)[j]*h
        plus,minus=model.evaluate(d+delta,False),model.evaluate(d-delta,False)
        fd=(plus["g"]-minus["g"])/(2*h)
        for i,name in enumerate(model.names):
            error=abs(fd[i]-a["jac"][i,j])/max(abs(fd[i]),abs(a["jac"][i,j]),1e-7)
            rows.append(dict(response=name,variable=["OD_mm","wall_mm"][j],adjoint=float(a["jac"][i,j]),
                             finite_difference=float(fd[i]),relative_error=float(error),passed=bool(error<3e-4)))
        massfd=(plus["mass"]-minus["mass"])/(2*h)
        error=abs(massfd-a["mass_grad"][j])/max(abs(massfd),1e-12)
        rows.append(dict(response="mass",variable=["OD_mm","wall_mm"][j],adjoint=float(a["mass_grad"][j]),
                         finite_difference=float(massfd),relative_error=float(error),passed=bool(error<1e-6)))
    if not all(r["passed"] for r in rows):
        raise RuntimeError("Optimizer gradient check failed: "+str([r for r in rows if not r["passed"]]))
    return rows


def solve_case(loads,study,nonlinear,progress=True):
    """Use the integrated scaled multistart optimizer for both static models."""
    from optimize_arm import solve
    validate(loads,study)
    return solve(TubeModel(loads,study,nonlinear),progress=progress)


def fingerprint(loads,study):
    return dict(loads_sha256=hashlib.sha256(json.dumps(loads,sort_keys=True).encode()).hexdigest(),
                study_sha256=hashlib.sha256(json.dumps(study,sort_keys=True).encode()).hexdigest(),
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ["arm.py","vibration.py","optimize.py","optimize_arm.py"]})


def run_study(loads,study,progress=True):
    validate(loads,study)
    t0=time.perf_counter()
    check_point=np.mean([study["outer_diameter_bounds_mm"],study["wall_bounds_mm"]],axis=1)
    check_point[0]=max(check_point[0],30.)
    check_point[0]=min(check_point[0],study["outer_diameter_bounds_mm"][1]*.95)
    model=TubeModel(loads,study,True)
    checks=gradient_check(model,check_point)
    if progress:print(f"Optimizer preflight: {len(checks)} mass/constraint derivatives passed",flush=True)
    without=solve_case(loads,study,False,progress)
    with_nl=solve_case(loads,study,True,progress)
    from optimize_arm import check_shadow_prices
    with_nl["shadow_checks"]=check_shadow_prices(model,with_nl)
    best=with_nl["optimum"]
    d=np.array([best["outer_diameter_mm"],best["wall_mm"]])
    refinement=[]
    for factor in [1,2,4]:
        refined=copy.deepcopy(study)
        refined["linear_elements"]*=factor;refined["nonlinear_elements"]*=factor
        r=TubeModel(loads,refined,True).summary(d)
        refinement.append(dict(linear_elements=refined["linear_elements"],nonlinear_elements=refined["nonlinear_elements"],
                               first_bending_frequency_Hz=r["first_bending_frequency_Hz"],
                               minimum_modeled_safety_factor=r["minimum_modeled_safety_factor"],
                               nonlinear_tip_deflection_mm=r["nonlinear_tip_deflection_mm"],
                               max_constraint=max(r["constraint_residuals"].values()),feasible=r["feasible"]))
    if not all(r["feasible"] and r["minimum_modeled_safety_factor"]>=study["required_safety_factor"]-1e-6 for r in refinement):
        raise RuntimeError("Optimal design failed refined-mesh feasibility check")
    rounding=np.array([study["rounding_mm"]["outer_diameter"],study["rounding_mm"]["wall"]])
    dr=np.ceil((d-1e-10)/rounding)*rounding
    rounded=None
    if np.all(dr<=np.array([study["outer_diameter_bounds_mm"][1],study["wall_bounds_mm"][1]])):
        rounded=model.summary(dr)
    initial=model.summary(np.array(study["reference_design_mm"],float))
    reference=with_nl["fixed_wall_reference"]["point"]
    return dict(status="converged_preliminary_model",provenance=fingerprint(loads,study),
                loads=loads,study=study,gradient_checks=checks,initial_design=initial,
                resolved_minimum_bending_frequency_Hz=model.minimum_frequency,
                fixed_wall_mass_reduction_percent=(100*(1-best["tube_mass_g"]/(reference["mass_kg"]*1000))
                                                   if reference else None),
                without_nonlinear=without,with_nonlinear=with_nl,rounded_design=rounded,
                mesh_refinement=refinement,
                mass_reduction_percent=100*(1-best["tube_mass_g"]/initial["tube_mass_g"]),
                nonlinear_mass_change_percent=100*(best["tube_mass_g"]/without["optimum"]["tube_mass_g"]-1),
                elapsed_s=time.perf_counter()-t0,
                conclusion="SF >= 2 is satisfied only for the modeled nominal linear and nonlinear tensile/compressive/shear checks using the recorded estimated material. Overall structural safety is unresolved.")


def export(result,out):
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    rows=result["with_nonlinear"]["history"]
    with (out/"history.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Directly runnable analyze.py input at the continuous optimum.
    loads=copy.deepcopy(result["loads"]);best=result["with_nonlinear"]["optimum"];p=best["material_properties"]
    loads["required_safety_factor"]=result["study"]["required_safety_factor"]
    loads["material"]["name"]=result["study"]["material"]["name"]
    loads["material"].update(E_bending_Pa=p["E_Pa"],design_tensile_strength_Pa=p["tensile_strength_Pa"],
                            design_compressive_strength_Pa=p["compressive_strength_Pa"],design_shear_strength_Pa=p["shear_strength_Pa"])
    loads["material"]["optimizer_estimate_basis"]=result["study"]["material"]
    A,_=section(best["outer_diameter_mm"]/1000,best["wall_mm"]/1000)
    loads["candidates"]=[dict(name="Continuous optimized tube - estimated material",
                            outer_diameter_m=best["outer_diameter_mm"]/1000,wall_m=best["wall_mm"]/1000,
                            mass_per_length_kg_m=float(p["density_kg_m3"]*A),source="optimizer_inputs.json estimated material profile")]
    (out/"optimal_inputs.json").write_text(json.dumps(loads,indent=2)+"\n",encoding="utf-8")
    os.environ.setdefault("MPLCONFIGDIR",str(ROOT/".mpl-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    for start in sorted(set(r["start"] for r in rows)):
        r=[v for v in rows if v["start"]==start]
        axes[0].plot([v["iteration"] for v in r],[v["mass_g"] for v in r],'.-',label=f"Start {start+1}")
        axes[1].plot([v["iteration"] for v in r],[max(v["max_constraint"],1e-12) for v in r],'.-')
    axes[0].set(xlabel="SLSQP iteration",ylabel="Tube mass (g)",title="Uniform tube mass optimization")
    axes[1].set(xlabel="SLSQP iteration",ylabel="Positive constraint violation",yscale="log",title="All enabled FE constraints")
    axes[0].legend();fig.tight_layout();fig.savefig(out/"convergence.png",dpi=150);plt.close(fig)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs",type=Path,default=ROOT/"inputs.json")
    p.add_argument("--study",type=Path,default=ROOT/"optimizer_inputs.json")
    p.add_argument("--output","--output-dir",type=Path,default=ROOT/"output"/"optimization")
    p.add_argument("--update-report",action="store_true",help="Regenerate report from default study/result paths after optimizing")
    args=p.parse_args(argv)
    loads=json.loads(args.inputs.read_text());study=json.loads(args.study.read_text())
    try:
        result=run_study(loads,study)
    except (ValueError,RuntimeError) as e:
        p.exit(2,str(e)+"\n")
    export(result,args.output.resolve())
    best=result["with_nonlinear"]["optimum"]
    print(f"Optimum: OD {best['outer_diameter_mm']:.6f} mm; wall {best['wall_mm']:.6f} mm; mass {best['tube_mass_g']:.4f} g")
    print(f"Minimum modeled SF: {best['minimum_modeled_safety_factor']:.6f}; mass change vs starting design {-result['mass_reduction_percent']:+.3f}%")
    if result['fixed_wall_mass_reduction_percent'] is not None:
        print(f"Mass saved vs feasible fixed-wall reference: {result['fixed_wall_mass_reduction_percent']:.3f}%")
    print(result["conclusion"])
    if args.update_report:
        if args.output.resolve()!=ROOT/"output"/"optimization" or args.study.resolve()!=ROOT/"optimizer_inputs.json" or args.inputs.resolve()!=ROOT/"inputs.json":
            p.exit(2,"Report integration uses default input, study and output paths. Custom run was saved successfully.\n")
        subprocess.run([sys.executable,str(ROOT/"run.py")],cwd=ROOT,check=True)
    return result


if __name__=="__main__":main()
