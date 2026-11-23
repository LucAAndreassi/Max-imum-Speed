"""Independent checks of the new optimizer and density-dependent gradients."""
import contextlib
import copy
import io
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from scipy.optimize import brentq

from optimize import Material, TubeModel, fingerprint, gradient_check, minimum_frequency, solve_case, validate
from optimize_arm import check_shadow_prices, solve

ROOT=Path(__file__).resolve().parent


class OptimizerTests(unittest.TestCase):
    def setUp(self):
        self.loads=json.loads((ROOT/"inputs.json").read_text())
        self.study=json.loads((ROOT/"optimizer_inputs.json").read_text())
        # Original strength-only fixtures remain independent of the new
        # default rotor-margin requirement; frequency has dedicated tests.
        self.study["minimum_bending_frequency_rotor_multiple"]=None

    def test_rotor_frequency_requirement_tracks_operating_rpm(self):
        self.study["minimum_bending_frequency_rotor_multiple"]=1.
        self.assertAlmostEqual(minimum_frequency(self.loads,self.study),3500/60)
        self.loads["rpm"]=4200.
        self.assertEqual(TubeModel(self.loads,self.study,False).minimum_frequency,70.)
        self.study["minimum_bending_frequency_Hz"]=50.
        with self.assertRaisesRegex(ValueError,"not both"):
            validate(self.loads,self.study)

    def test_actual_density_dependent_static_solution(self):
        d=np.array([28.,2.4]); model=TubeModel(self.loads,self.study,False)
        a=model.evaluate(d)
        D,t=d/1000;di=D-2*t
        # Independent closed-form section/deflection, not calls to FE element code.
        A=np.pi*(D*D-di*di)/4;I=np.pi*(D**4-di**4)/64
        p=model.material.at(d[1]);mu=p["density_kg_m3"]*A
        L,g=model.L,model.g;P=model.T-model.mt*g
        expected=P*L**3/(3*p["E_Pa"]*I)-mu*g*L**4/(8*p["E_Pa"]*I)
        self.assertAlmostEqual(a["u"][-2]/expected,1,places=8)
        self.assertAlmostEqual(a["mass"]/(mu*L),1,places=12)

    def test_complete_adjoint_with_density_weight_and_optional_limits(self):
        self.study["deflection_limit_mm"]=20.
        self.study["minimum_bending_frequency_Hz"]=40.
        rows=gradient_check(TubeModel(self.loads,self.study,True),[32.,2.5])
        self.assertTrue(all(r["passed"] for r in rows))
        self.assertIn("bending_frequency",{r["response"] for r in rows})
        self.assertIn("nonlinear_deflection",{r["response"] for r in rows})

    def test_table_interpolation_and_its_state_gradient(self):
        # Synthetic tables test interpolation mechanics; not report material data.
        self.study["minimum_bending_frequency_Hz"]=40.
        base=self.study["material"]
        keys=["E_Pa","density_kg_m3","tensile_strength_Pa","compressive_strength_Pa","shear_strength_Pa"]
        base["thickness_samples"]=[dict(wall_mm=t,**{k:base[k]*scale for k in keys})
                                    for t,scale in [(2.,.9),(3.2,1.1)]]
        m=Material(base)
        for key in keys:self.assertAlmostEqual(m.at(2.6)[key]/base[key],1.,places=12)
        with self.assertRaises(ValueError):m.at(3.3)
        self.assertTrue(all(r["passed"] for r in gradient_check(TubeModel(self.loads,self.study,True),[32.,2.5])))

    def test_nonconservative_follower_tangent_is_nonsymmetric(self):
        model=TubeModel(self.loads,self.study,True)
        nl=model.nonlinear_model(np.array([18.,2.1]))
        theta,_=nl.solve()
        H=nl.residual(theta)[1][1:,1:]
        self.assertGreater(np.linalg.norm(H-H.T),1e-5)

    def test_linear_optimum_against_independent_strength_bound(self):
        # Constant-property, bending-dominated fixture. The independent optimum
        # lies at the minimum wall. KS aggregation may increase mass slightly.
        self.study["starting_designs_mm"]=[[40.,2.5]]
        r=solve_case(self.loads,self.study,False,False)
        opt=r["optimum"];p=self.study["material"]
        t=self.study["wall_bounds_mm"][0]/1000;L=self.loads["length_m"];g=self.loads["gravity_m_s2"]
        P=self.loads["thrust_N"]-self.loads["tip_mass_kg"]*g
        def mass_and_util(dmm):
            D=dmm/1000;di=D-2*t
            A=np.pi*(D*D-di*di)/4;I=np.pi*(D**4-di**4)/64
            mu=p["density_kg_m3"]*A
            M=P*L-mu*g*L*L/2
            return mu*L,self.study["required_safety_factor"]*M*D/(2*I*p["tensile_strength_Pa"])
        exact_D=brentq(lambda d:mass_and_util(d)[1]-1,*self.study["outer_diameter_bounds_mm"])
        exact_mass=mass_and_util(exact_D)[0]*1000
        self.assertAlmostEqual(opt["wall_mm"],t*1000,places=6)
        self.assertGreaterEqual(opt["minimum_modeled_safety_factor"],2.-1e-6)
        self.assertGreaterEqual(opt["tube_mass_g"],exact_mass-1e-5)
        self.assertLess(opt["tube_mass_g"]/exact_mass,1.01)
        self.assertLess(r["kkt"]["stationarity_inf_norm"],1e-5)
        self.assertEqual(r["settings"]["engine"],"optimize_arm")
        self.assertTrue(r["kkt"]["passed"])
        self.assertTrue(r["starts"][r["best_run"]]["accepted"])
        self.assertIn("optimize_arm.py",fingerprint(self.loads,self.study)["source_sha256"])

    def test_solver_success_without_kkt_stationarity_is_rejected(self):
        # A feasible oversized tube is not a stationary mass minimum. A solver
        # success flag alone must never make it eligible for design selection.
        fake=SimpleNamespace(x=np.array([1.,1.]),success=True,status=0,nit=1,message="synthetic success")
        with patch("optimize_arm.minimize",return_value=fake):
            with self.assertRaisesRegex(RuntimeError,"KKT-passing"):
                solve(TubeModel(self.loads,self.study,False),progress=False)

    def test_failed_start_is_recorded_and_other_starts_continue(self):
        model=TubeModel(self.loads,self.study,False)
        actual=model.evaluate
        def evaluate(d,gradients=True):
            if np.allclose(d,[10.,2.01],rtol=0,atol=1e-12):
                raise RuntimeError("synthetic FEM failure at one start")
            return actual(d,gradients)
        with patch.object(model,"evaluate",side_effect=evaluate):
            result=solve(model,progress=False)
        failed=next(r for r in result["starts"] if r["label"]=="Small tube")
        self.assertFalse(failed["accepted"])
        self.assertIn("synthetic FEM failure",failed["message"])
        self.assertTrue(result["optimum"]["feasible"])

    def test_frequency_reference_and_shadow_price_against_resolve(self):
        # Exercise the supplied script's reference-crossing logic and frequency
        # price using the current density-dependent modal FEM and adjoint.
        self.study["minimum_bending_frequency_Hz"]=60.
        model=TubeModel(self.loads,self.study,False)
        result=solve(model,progress=False)
        reference=result["fixed_wall_reference"]
        self.assertEqual(reference["threshold_kind"],"crossing")
        self.assertAlmostEqual(reference["frequency_threshold"]["constraint_residuals"]["bending_frequency"],0.,places=6)
        self.assertLessEqual(max(reference["point"]["constraint_residuals"].values()),1e-10)
        d=np.array([result["optimum"]["outer_diameter_mm"],result["optimum"]["wall_mm"]])
        self.assertGreaterEqual(model.evaluate(d,False)["frequency"],60.-1e-5)
        checks=check_shadow_prices(model,result)
        self.assertIn("bending_frequency",{c["name"] for c in checks})
        for check in checks:
            self.assertIsNotNone(check["change_kg"],check["unavailable_reason"])
            self.assertAlmostEqual(check["predicted_kg"]/check["change_kg"],1.,delta=.015)

    def test_infeasible_problem_cannot_be_reported_as_success(self):
        self.study["material"]["tensile_strength_Pa"]=1e5
        self.study["starting_designs_mm"]=[[40.,2.5]]
        with self.assertRaises(RuntimeError):solve_case(self.loads,self.study,False,False)

    def test_strict_wall_requirement_and_missing_material_validation(self):
        self.study["wall_bounds_mm"][0]=2.
        with self.assertRaises(ValueError):validate(self.loads,self.study)
        self.study["wall_bounds_mm"][0]=2.01
        self.study["material"]["compressive_strength_Pa"]=None
        with self.assertRaises(ValueError):validate(self.loads,self.study)


if __name__=="__main__":unittest.main(verbosity=2)
