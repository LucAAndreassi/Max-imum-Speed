"""Run the supplied 3D frame solver at the current optimized design.

python compare_drone_arm.py
The supplied drone_arm.py imports a missing arm_inputs module. This runner
provides its defaults in memory from the existing project cards, then passes
every physical Arm input explicitly. It does not modify the supplied solver.
Poisson's ratio is an auxiliary assumption for torsion, checked independently
of the bending comparison; it is not added to the material card.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

import numpy as np
from scipy.sparse import save_npz

from arm import Beam
from vibration import tip_mass_frequencies

ROOT = Path(__file__).resolve().parent


def load_solver(parameters, seed=20261123):
    """Supply only the missing defaults module while loading the actual file."""
    defaults = ModuleType("arm_inputs")
    mapping = dict(LENGTH_IN=parameters["length"]/.0254,
                   INITIAL_OD_IN=parameters["outer_diameter"]/.0254,
                   INITIAL_WALL_IN=parameters["wall"]/.0254,
                   YOUNG_MODULUS_GPA=parameters["young"]/1e9,
                   POISSON_RATIO=parameters["poisson"], DENSITY_KG_M3=parameters["density"],
                   MOTOR_MASS_COMBINED_LBM=parameters["tip_mass"]/.45359237,
                   THRUST_PER_MOTOR_LBF=parameters["thrust_per_motor"]/(.45359237*9.80665),
                   ROTOR_RPM=parameters["rpm"], BEAM_ELEMENTS=parameters["elements"],
                   GRAVITY_M_S2=parameters["gravity"], EIGENSOLVER_SEED=seed)
    defaults.__dict__.update(mapping)
    previous = sys.modules.get("arm_inputs")
    sys.modules["arm_inputs"] = defaults
    name = "drone_arm_comparison_solver"
    spec = importlib.util.spec_from_file_location(name, ROOT/"drone_arm.py")
    solver = importlib.util.module_from_spec(spec)
    sys.modules[name] = solver
    try:
        spec.loader.exec_module(solver)
    finally:
        if previous is None:
            del sys.modules["arm_inputs"]
        else:
            sys.modules["arm_inputs"] = previous
    return solver, mapping


def run_case(solver, arm):
    model = solver.build_model(arm)
    u, reaction = solver.statics(model)
    frequencies, modes = solver.modes(model, count=8)
    free = slice(6, None)
    K, M = model.stiffness[free, free], model.mass[free, free]
    phi = modes[free]
    modal_residual = K@phi-(M@phi)*(2*np.pi*frequencies)**2
    relative_residual = np.linalg.norm(modal_residual, axis=0)/np.linalg.norm(K@phi, axis=0)
    # Circular, isotropic effective sections have two degenerate bending modes;
    # ARPACK can return any mixture of their two bending planes.
    families = dict(axial=[0], torsion=[3], lateral_bending=[1, 5], vertical_bending=[2, 4])
    fractions = []
    for j in range(len(frequencies)):
        energy = {}
        for family, local in families.items():
            v = np.zeros_like(u)
            dofs = np.array([6*n+d for n in range(len(model.x)) for d in local])
            v[dofs] = modes[dofs, j]
            energy[family] = float(v@(model.stiffness@v))
        total = sum(energy.values())
        fractions.append({k:v/total for k,v in energy.items()})
    result = dict(elements=arm.elements, nodes=len(model.x), total_dofs=len(u), free_dofs=len(u)-6,
                  tube_mass_kg=float(arm.tube_mass), tip_z_mm=float(u[-4]*1000),
                  tip_slope_rad=float(-u[-2]), peak_bending_stress_MPa=float(solver.bending_stress(model,u).max()/1e6),
                  frequencies_Hz=frequencies.tolist(), modal_strain_energy_fractions=fractions,
                  modal_relative_residuals=relative_residual.tolist(),
                  equilibrium_residual_inf=float(np.linalg.norm(reaction[6:],np.inf)),
                  vertical_force_balance_N=float(reaction[2]+model.load[2::6].sum()),
                  vertical_tip_stiffness_N_m=float(solver.tip_stiffness(model)))
    if not np.isfinite(u).all() or np.max(relative_residual)>1e-6:
        raise RuntimeError("Nonfinite response or failed frame eigenpair residual check")
    np.testing.assert_allclose(phi.T@(M@phi),np.eye(len(frequencies)),atol=1e-9)
    return result, model, u, modes


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poisson",type=float,default=.3,
                        help="Assumed auxiliary Poisson ratio for the new solver's torsional block")
    parser.add_argument("--output",type=Path,default=ROOT/"output/drone-arm-comparison")
    args=parser.parse_args(argv)
    study_result=json.loads((ROOT/"output/optimization/results.json").read_text())
    loads=study_result["loads"]
    best=study_result["with_nonlinear"]["optimum"]
    material=best["material_properties"]
    parameters=dict(length=loads["length_m"],outer_diameter=best["outer_diameter_mm"]/1000,
                    wall=best["wall_mm"]/1000,young=material["E_Pa"],poisson=args.poisson,
                    density=material["density_kg_m3"],tip_mass=loads["tip_mass_kg"],
                    thrust_per_motor=loads["thrust_N"]/2,rpm=loads["rpm"],
                    elements=study_result["study"]["linear_elements"],gravity=loads["gravity_m_s2"])
    solver, defaults=load_solver(parameters)
    arm=solver.Arm(**parameters)
    current,model,u,modes=run_case(solver,arm)
    old=Beam(arm.length,arm.young,arm.density*arm.area,arm.tip_mass,
             2*arm.thrust_per_motor,arm.gravity,arm.elements)
    old_u,_,_=old.solve(arm.outer_diameter,arm.wall)
    old_fe=old.fe_modes(arm.outer_diameter,arm.wall,4)[0]
    naca=tip_mass_frequencies(arm.length,arm.young*arm.inertia,arm.density*arm.area,arm.tip_mass,4)
    # Independent closed-form linear static check, including tube self-weight.
    P=2*arm.thrust_per_motor-arm.tip_mass*arm.gravity
    w=arm.density*arm.area*arm.gravity
    exact_tip=(P*arm.length**3/3-w*arm.length**4/8)/(arm.young*arm.inertia)*1000
    refinement=[]
    for factor in [1,2,4]:
        case=current if factor==1 else run_case(solver,replace(arm,elements=arm.elements*factor))[0]
        refinement.append(case)
    poisson_sweep=[]
    for nu in sorted(set([0.,args.poisson,.49])):
        case=current if nu==arm.poisson else run_case(solver,replace(arm,poisson=nu))[0]
        poisson_sweep.append(dict(poisson=nu,tip_z_mm=case["tip_z_mm"],first_frequency_Hz=case["frequencies_Hz"][0],
                                  first_mode_bending_fraction=sum(case["modal_strain_energy_fractions"][0][k]
                                                                  for k in ["lateral_bending","vertical_bending"])))
    previous=json.loads((ROOT/"output/optimized-fe/results.json").read_text())
    delta_tip=current["tip_z_mm"]-old_u[-2]*1000
    delta_f=current["frequencies_Hz"][0]-naca[0]
    np.testing.assert_allclose(current["tip_z_mm"],exact_tip,rtol=1e-7)
    np.testing.assert_allclose(current["frequencies_Hz"][:2],naca[0],rtol=2e-6)
    np.testing.assert_allclose([v["tip_z_mm"] for v in poisson_sweep],current["tip_z_mm"],rtol=1e-10)
    np.testing.assert_allclose([v["first_frequency_Hz"] for v in poisson_sweep],current["frequencies_Hz"][0],rtol=1e-9)
    assert all(v["first_mode_bending_fraction"]>1-1e-10 for v in poisson_sweep)
    paths=["drone_arm.py","compare_drone_arm.py","arm.py","vibration.py","output/optimization/results.json"]
    result=dict(parameters_si=asdict(arm),defaults_dependency=defaults,
                notes=["Supplied drone_arm.py is unmodified; its missing arm_inputs defaults are supplied in memory",
                       "Gross thrust is split over two motors to avoid doubling the existing total load",
                       "Poisson ratio is an unverified auxiliary assumption used only by the uncoupled torsion block",
                       "New code is a linear 3D Euler-Bernoulli frame, not the finite-rotation follower-thrust solver",
                       "First two frequencies form the expected two-plane bending pair"],
                reference=dict(linear_tip_z_mm=float(old_u[-2]*1000),fe_bending_frequencies_Hz=old_fe.tolist(),
                               naca_bending_frequencies_Hz=naca.tolist(),closed_form_tip_z_mm=float(exact_tip),
                               prior_export_linear_tip_z_mm=previous["linear"]["tip_deflection_m"]*1000,
                               prior_export_nonlinear_tip_z_mm=previous["nonlinear"]["tip_z_m"]*1000),
                new_frame=current,mesh_refinement=refinement,poisson_sensitivity=poisson_sweep,
                differences=dict(linear_tip_mm=float(delta_tip),linear_tip_percent=float(delta_tip/(old_u[-2]*1000)*100),
                                 first_frequency_Hz=float(delta_f),first_frequency_percent=float(delta_f/naca[0]*100)),
                source_sha256={str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    np.save(args.output/"displacements.npy",u);np.save(args.output/"modes.npy",modes)
    save_npz(args.output/"K.npz",model.stiffness);save_npz(args.output/"M.npz",model.mass)
    np.save(args.output/"F.npy",model.load)
    np.savetxt(args.output/"nodes.csv",np.column_stack((model.x,u.reshape(-1,6))),delimiter=",",
               header="x_m,ux_m,uy_m,uz_m,rx_rad,ry_rad,rz_rad",comments="")
    print(f"Same geometry: {arm.outer_diameter*1000:.9f} mm OD x {arm.wall*1000:.9f} mm wall")
    print(f"E={arm.young/1e9:g} GPa, density={arm.density:g} kg/m3, total thrust={2*arm.thrust_per_motor:.9f} N")
    print(f"New frame: {current['nodes']} nodes, {current['free_dofs']} free DOFs")
    print(f"Linear tip: existing {old_u[-2]*1000:.12f} mm; new {current['tip_z_mm']:.12f} mm ({result['differences']['linear_tip_percent']:+.6g}%)")
    print(f"First frequency: NACA {naca[0]:.12f} Hz; new {current['frequencies_Hz'][0]:.12f} Hz ({result['differences']['first_frequency_percent']:+.6g}%)")
    print("New frame modes (Hz): "+", ".join(f"{f:.6f}" for f in current["frequencies_Hz"]))
    print("Poisson sensitivity, independent static solution, modal residuals and mesh refinement checked.")
    print(f"Saved comparison: {args.output/'results.json'}")
    return result


if __name__=="__main__":
    main()
