"""Run one parameterized beam FE model without rebuilding the report.

Example:
  python analyze.py --od-mm 30 --wall-mm 1.5 --line-mass-kg-m 0.1938 --nonlinear
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import save_npz

from arm import Beam, Elastica, section
from check_design import evaluate_strength

ROOT = Path(__file__).resolve().parent
FT_TO_M = 0.3048
LBF_TO_N = 4.4482216152605
LB_TO_KG = 0.45359237


def positive(value):
    x = float(value)
    if not np.isfinite(x) or x <= 0:
        raise argparse.ArgumentTypeError("Must be finite and greater than zero")
    return x


def nonnegative(value):
    x = float(value)
    if not np.isfinite(x) or x < 0:
        raise argparse.ArgumentTypeError("Must be finite and nonnegative")
    return x


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", type=Path, default=ROOT / "inputs.json", help="Base material/load JSON")
    p.add_argument("--candidate", type=int, default=0, help="Zero-based candidate from JSON; default 0")
    p.add_argument("--length-ft", type=positive, help="Clamp-to-tip length in ft")
    p.add_argument("--thrust-lbf", type=nonnegative, help="Combined upward thrust in lbf")
    p.add_argument("--tip-mass-lb", type=nonnegative, help="Total carried assembly mass in lb")
    p.add_argument("--od-mm", type=positive, help="Uniform outer diameter in mm")
    p.add_argument("--wall-mm", type=positive, help="Uniform wall thickness in mm")
    p.add_argument("--E-gpa", type=positive, help="Effective longitudinal/bending modulus in GPa")
    mass = p.add_mutually_exclusive_group()
    mass.add_argument("--line-mass-kg-m", type=positive, help="Known tube mass per length in kg/m")
    mass.add_argument("--density-kg-m3", type=positive, help="User-specified density; calculates line mass as rho*A")
    p.add_argument("--elements", type=int, help="Number of linear bending FE elements")
    p.add_argument("--nonlinear", action="store_true", help="Also solve FE large rotation with tip-normal follower thrust")
    p.add_argument("--vibration-method",choices=["naca_tn_3459","fem"],default="naca_tn_3459",
                   help="Default: NACA TN 3459 eqs. 13-14; fem retains the matrix eigenproblem for comparison")
    p.add_argument("--nl-elements", type=int, default=80, help="Nonlinear angle FE elements; default 80")
    p.add_argument("--tensile-mpa", type=positive, help="Applicable ultimate tensile strength in MPa")
    p.add_argument("--compression-mpa", type=positive, help="Applicable ultimate compressive strength in MPa")
    p.add_argument("--shear-mpa", type=positive, help="Applicable ultimate shear strength in MPa")
    p.add_argument("--output", type=Path, default=ROOT / "output" / "latest", help="Run output directory, overwritten on rerun")
    args = p.parse_args(argv)
    return p, args


def run(argv=None):
    parser, args = parse_args(argv)
    config = json.loads(args.inputs.read_text(encoding="utf-8"))
    if config.get("applied_torque_Nm", 0) != 0:
        parser.error("This bending-only model requires zero applied torque")
    if not 0 <= args.candidate < len(config["candidates"]):
        parser.error("Candidate index is outside the inputs.json candidate list")
    c = config["candidates"][args.candidate]
    D = args.od_mm / 1000 if args.od_mm is not None else c["outer_diameter_m"]
    t = args.wall_mm / 1000 if args.wall_mm is not None else c["wall_m"]
    if not np.isfinite(D) or not np.isfinite(t) or D <= 2*t or t <= 0:
        parser.error("Require outer diameter > 2 * wall thickness > 0")
    A, I = section(D, t)
    changed = not (np.isclose(D,c["outer_diameter_m"],rtol=1e-12,atol=0) and
                   np.isclose(t,c["wall_m"],rtol=1e-12,atol=0))
    if changed and args.density_kg_m3 is None and args.line_mass_kg_m is None:
        parser.error("Changed geometry needs --density-kg-m3 or --line-mass-kg-m; stock tube mass cannot be reused for a different section")
    if args.density_kg_m3 is not None:
        mu = args.density_kg_m3 * A
        mass_basis = "User density multiplied by exact annular area"
    elif args.line_mass_kg_m is not None:
        mu = args.line_mass_kg_m
        mass_basis = "User-specified line mass"
    else:
        mu = c["mass_per_length_kg_m"]
        mass_basis = "Configured line mass for selected tube (see candidate source)"
    L = args.length_ft * FT_TO_M if args.length_ft is not None else config["length_m"]
    thrust = args.thrust_lbf * LBF_TO_N if args.thrust_lbf is not None else config["thrust_N"]
    mt = args.tip_mass_lb * LB_TO_KG if args.tip_mass_lb is not None else config["tip_mass_kg"]
    E = args.E_gpa * 1e9 if args.E_gpa is not None else config["material"]["E_bending_Pa"]
    n = args.elements if args.elements is not None else config["elements"]
    g = config["gravity_m_s2"]
    for name,value in [("length",L),("modulus",E),("line mass",mu)]:
        if not np.isfinite(value) or value <= 0: parser.error(f"{name} must be positive and finite")
    if not isinstance(n,int) or n < 1 or args.nl_elements < 1:
        parser.error("Element counts must be positive integers")
    if not all(np.isfinite(value) and value >= 0 for value in [mt,thrust,g]):
        parser.error("Mass, thrust and gravity must be finite and nonnegative")
    material = dict(config["material"])
    for arg,key in [(args.tensile_mpa,"design_tensile_strength_Pa"),
                    (args.compression_mpa,"design_compressive_strength_Pa"),
                    (args.shear_mpa,"design_shear_strength_Pa")]:
        if arg is not None: material[key] = arg * 1e6
    beam = Beam(L,E,mu,mt,thrust,g,n)
    u,reactions,(K,M,F,free,_) = beam.solve(D,t)
    mode_count=min(4,len(free))
    frequencies,modes = beam.modes(D,t,mode_count,method=args.vibration_method)
    fe_frequencies,fe_modes = beam.fe_modes(D,t,mode_count)
    analytical = beam.vibration_summary(D,t,mode_count)
    x = np.linspace(0,L,n+1)
    q = -mu*g
    # Recover the uniform-load moment field from FE support reactions and equilibrium.
    critical = [0.,L]
    if q != 0 and 0 < -reactions[0]/q < L: critical.append(-reactions[0]/q)
    def moment(xx): return -reactions[1]+reactions[0]*xx+q*xx**2/2
    peak_M = max(abs(moment(xx)) for xx in critical)
    sigma = peak_M*D/(2*I)
    R,ri = D/2,D/2-t
    shear = max(abs(reactions[0]),abs(reactions[0]+q*L))*(R**3-ri**3)/(3*I*t)
    residual = float(np.linalg.norm((K@u-F)[free],np.inf))
    result = {
        "model":"Uniform clamped Euler-Bernoulli bending FE; 2 DOF per node [w,slope]",
        "status":"Analysis run; not an optimized or certified design",
        "parameters_SI":dict(length_m=L,thrust_N=thrust,tip_mass_kg=mt,E_Pa=E,
                              outer_diameter_m=D,wall_m=t,mass_per_length_kg_m=float(mu),gravity_m_s2=g),
        "mass_basis":mass_basis,"selected_base_candidate":c["name"],
        "mesh":dict(elements=n,nodes=n+1,total_dofs=len(u),free_dofs=len(free)),
        "linear":dict(tip_deflection_m=float(u[-2]),tip_slope_rad=float(u[-1]),
                      root_reaction_force_N=float(reactions[0]),root_reaction_moment_Nm=float(reactions[1]),
                      peak_bending_stress_Pa=float(sigma),homogeneous_shear_screen_Pa=float(shear),
                      bending_frequencies_Hz=frequencies.tolist(),free_dof_residual_N_or_Nm=residual),
        "tube_mass_kg":float(mu*L),"total_modeled_mass_kg":float(mu*L+mt),
        "material_for_strength_screen":material,
        "strength_checks":evaluate_strength(sigma,shear,material,config["required_safety_factor"]),
        "vibration":dict(method=analytical["method"] if args.vibration_method=="naca_tn_3459" else "consistent_mass_fem",
                         analytical=analytical,fe_frequencies_Hz=fe_frequencies.tolist(),
                         fe_relative_frequency_error=(fe_frequencies/np.array(analytical["frequencies_Hz"])-1).tolist(),
                         modes_file="modes.npy",fe_modes_file="fe_modes.npy",
                         mode_normalization="Continuum mass integral plus point tip mass" if args.vibration_method=="naca_tn_3459" else "FE mass matrix",
                         note="Analytical nodal samples are not discrete K,M eigenvectors; fe_modes.npy contains the matrix eigenvectors"),
        "nonlinear":None,
        "limitations":["Effective beam properties, not a shell or ply-level laminate model",
                       "Frequencies are unloaded bending modes with point tip mass; rotor inertia and flutter omitted",
                       "Homogeneous transverse-shear screen and nominal maximum-stress checks do not cover local shell or joint failure"]
    }
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    with (out/"nodes.csv").open("w",newline="",encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node","x_m","w_m","slope_rad","recovered_moment_Nm","nominal_bending_stress_Pa"])
        for i,xx in enumerate(x):
            mm = moment(xx)
            writer.writerow([i,float(xx),float(u[2*i]),float(u[2*i+1]),float(mm),float(abs(mm)*D/(2*I))])
    mesh = dict(nodes_xyz_m=[[float(xx),0.,0.] for xx in x],
                connectivity=[[i,i+1] for i in range(n)],
                dof_order="[w0,slope0,w1,slope1,...]",fixed_dofs=[0,1],
                free_dofs=free.tolist(),outer_diameter_m=D,wall_m=t)
    (out/"mesh.json").write_text(json.dumps(mesh,indent=2)+"\n",encoding="utf-8")
    save_npz(out/"K.npz",K); save_npz(out/"M.npz",M)
    np.save(out/"F.npy",F); np.save(out/"u.npy",u); np.save(out/"modes.npy",modes)
    np.save(out/"fe_modes.npy",fe_modes);np.save(out/"fe_frequencies.npy",fe_frequencies)
    if args.nonlinear:
        nl = Elastica(L,E*I,-mt*g,q,elements=args.nl_elements,follower_force=thrust)
        theta,history = nl.solve()
        points = nl.coordinates(theta)
        nr = float(np.linalg.norm(nl.residual(theta)[0][1:],np.inf))
        result["nonlinear"] = dict(elements=args.nl_elements,tip_x_m=float(points[-1,0]),
                                   tip_z_m=float(points[-1,1]),tip_angle_rad=float(theta[-1]),
                                   residual_Nm=nr,newton_history=history)
        with (out/"nonlinear_nodes.csv").open("w",newline="",encoding="utf-8") as f:
            writer=csv.writer(f);writer.writerow(["node","s_m","x_m","z_m","theta_rad"])
            for i,((xx,zz),th) in enumerate(zip(points,theta)):
                writer.writerow([i,i*L/args.nl_elements,float(xx),float(zz),float(th)])
    elif (out/"nonlinear_nodes.csv").exists():
        # Avoid leaving stale nonlinear output from a prior run in the same directory.
        (out/"nonlinear_nodes.csv").write_text("Not computed for this run; rerun with --nonlinear.\n",encoding="utf-8")
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(f"FE model: {n} elements, {n+1} nodes, {len(free)} free DOF")
    print(f"Tube OD {D*1000:g} mm, wall {t*1000:g} mm, length {L:g} m")
    print(f"Tube mass: {mu*L*1000:.3f} g")
    print(f"Linear tip deflection: {u[-2]*1000:.6f} mm; nominal bending stress: {sigma/1e6:.6f} MPa")
    print(f"First bending frequency ({args.vibration_method}): {frequencies[0]:.6f} Hz; FE comparison {fe_frequencies[0]:.6f} Hz")
    if result["nonlinear"]:
        print(f"Follower-thrust nonlinear FE tip deflection: {result['nonlinear']['tip_z_m']*1000:.6f} mm")
    print(f"Nominal strength screen: {result['strength_checks']['nominal_maximum_stress_screen']['status']}")
    print(f"Saved mesh, matrices, loads, displacements, modes and results to: {out}")
    return result


if __name__ == "__main__":
    run()
