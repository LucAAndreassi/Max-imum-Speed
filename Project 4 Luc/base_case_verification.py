"""Physical base-case verification for the Project 4 Luc drone-arm model.

This module generates the numerical verification data reported in Section 4.
It intentionally complements verify.py: verify.py contains synthetic unit-style
fixtures, while this file evaluates the actual 40 mm OD x 2.5 mm wall base case
using the project inputs and optimizer material model.

The public entry point is::

    run_base_case_verification(loads, study)

It returns a JSON-serializable dictionary.  Running this file directly reads
inputs.json and optimizer_inputs.json and prints the dictionary as JSON.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import eigh

from arm import Beam, Elastica, hermite, section
from optimize import TubeModel, gradient_check
from vibration import tip_mass_frequencies

ROOT = Path(__file__).resolve().parent


def _rel_error(computed: float, reference: float) -> float:
    """Relative error, using absolute error when the reference is zero."""
    if reference == 0.0:
        return abs(float(computed))
    return abs(float(computed) - float(reference)) / abs(float(reference))


def _exact_linear_field(x, L, EI, P, q):
    """Exact EB displacement for a cantilever with tip force P and UDL q."""
    x = np.asarray(x)
    return (
        P * x**2 * (3.0 * L - x) / (6.0 * EI)
        + q * x**2 * (6.0 * L**2 - 4.0 * L * x + x**2) / (24.0 * EI)
    )


def _displacement_l2_error(beam: Beam, D: float, t: float, P: float, q: float) -> float:
    """Integrated L2 error against the exact tip-load + UDL cantilever field."""
    u, _, _ = beam.solve(D, t)
    _, I = section(D, t)
    EI = beam.E * I
    n = beam.elements
    le = beam.length / n
    gx, gw = np.polynomial.legendre.leggauss(5)
    error2 = 0.0
    for e in range(n):
        for xi, weight in zip((gx + 1.0) / 2.0, gw / 2.0):
            x = (e + xi) * le
            exact = _exact_linear_field(x, beam.length, EI, P, q)
            numerical = hermite(xi, le) @ u[2 * e : 2 * e + 4]
            error2 += weight * le * float((numerical - exact) ** 2)
    return float(np.sqrt(error2))


def _continuum_follower_bvp(L, EI, thrust, tip_mass, mu, g, initial_tip_angle=0.01):
    """Independent continuum BVP for the physical follower-thrust problem.

    State y=[theta, kappa, z].  The unknown parameter is theta(L), which fixes
    the follower-force direction.  This does not reuse the nonlinear FE residual.
    """
    mesh = np.linspace(0.0, L, 121)

    def ode(s, y, p):
        theta = y[0]
        theta_L = p[0]
        dead_down = tip_mass * g + mu * g * (L - s)
        theta_ss = (-thrust * np.cos(theta_L - theta) + dead_down * np.cos(theta)) / EI
        return np.vstack((y[1], theta_ss, np.sin(theta)))

    def bc(a, b, p):
        return np.array([a[0], b[1], a[2], b[0] - p[0]])

    # Small-rotation cantilever-shaped initial guess.  It is only an initial
    # iterate for the independent continuum equations, not an FE result.
    sbar = mesh / L
    theta_guess = initial_tip_angle * sbar * (2.0 - sbar)
    kappa_guess = (2.0 * initial_tip_angle / L) * (1.0 - sbar)
    z_guess = initial_tip_angle * L * (sbar**2 - sbar**3 / 3.0)
    guess = np.vstack((theta_guess, kappa_guess, z_guess))

    sol = solve_bvp(
        ode,
        bc,
        mesh,
        guess,
        p=np.array([initial_tip_angle]),
        tol=1e-10,
        max_nodes=20000,
        verbose=0,
    )
    if not sol.success:
        raise RuntimeError(f"Independent physical follower BVP failed: {sol.message}")
    return sol


def _taylor_test(model: TubeModel, d0: np.ndarray):
    """First-order Taylor remainder for [mass, constraints] at the base case."""
    base = model.evaluate(d0, gradients=True)
    f0 = np.r_[base["mass"], base["g"]]
    jac = np.vstack((base["mass_grad"], base["jac"]))

    direction = np.array([0.6, 0.8])
    scaled_direction = np.diag(d0) @ direction
    df = jac @ scaled_direction
    normalization = np.maximum(np.abs(f0), 1.0)

    hs = np.array([1e-2, 3e-3, 1e-3, 3e-4], float)
    rows = []
    for h in hs:
        trial = model.evaluate(d0 + h * scaled_direction, gradients=False)
        f1 = np.r_[trial["mass"], trial["g"]]
        remainder = (f1 - f0 - h * df) / normalization
        norm = float(np.linalg.norm(remainder))
        rows.append({"h": float(h), "normalized_remainder": norm})

    order = float(np.polyfit(np.log(hs), np.log([r["normalized_remainder"] for r in rows]), 1)[0])
    return {"direction": direction.tolist(), "rows": rows, "measured_order": order}


def run_base_case_verification(loads, study):
    """Generate the physical seven-rung verification dataset.

    Parameters
    ----------
    loads : dict
        Parsed inputs.json.
    study : dict
        Parsed optimizer_inputs.json.

    Returns
    -------
    dict
        JSON-serializable verification data for the 40 mm x 2.5 mm base case.
    """
    loads = copy.deepcopy(loads)
    study = copy.deepcopy(study)

    d0_mm = np.asarray(study.get("reference_design_mm", [40.0, 2.5]), dtype=float)
    if not np.allclose(d0_mm, [40.0, 2.5], rtol=0.0, atol=1e-12):
        raise ValueError("Section 4 base case requires reference_design_mm = [40.0, 2.5]")

    D = d0_mm[0] / 1000.0
    t = d0_mm[1] / 1000.0
    L = float(loads["length_m"])
    g = float(loads["gravity_m_s2"])
    T = float(loads["thrust_N"])
    mt = float(loads["tip_mass_kg"])
    E = float(study["material"]["E_Pa"])
    rho = float(study["material"]["density_kg_m3"])
    linear_elements = int(study["linear_elements"])

    A, I = section(D, t)
    EI = E * I
    mu = rho * A
    tube_mass = mu * L
    tube_weight = tube_mass * g
    P = T - mt * g
    q = -mu * g

    tolerances = {
        "V1": 1e-9,
        "V2": 1e-9,
        "V3_order_relative": 0.01,
        "V4": 2e-6,
        "V5_force_N": 1e-7,
        "V5_moment_Nm": 1e-7,
        "V6_constraint": 3e-4,
        "V6_mass": 1e-6,
        "V7": 1e-4,
    }

    ladder = []

    # ------------------------------------------------------------------ V1
    # Isolate the net tip load.  Distributed self-weight is suppressed.
    v1_beam = Beam(L, E, mu, 0.0, P, gravity=0.0, elements=linear_elements)
    u1, r1, _ = v1_beam.solve(D, t)
    v1_exact_disp = P * L**3 / (3.0 * EI)
    v1_exact_slope = P * L**2 / (2.0 * EI)
    v1_exact_moment = P * L
    v1 = {
        "tip_displacement_m": float(u1[-2]),
        "exact_tip_displacement_m": float(v1_exact_disp),
        "tip_displacement_relative_error": _rel_error(u1[-2], v1_exact_disp),
        "tip_slope_rad": float(u1[-1]),
        "exact_tip_slope_rad": float(v1_exact_slope),
        "tip_slope_relative_error": _rel_error(u1[-1], v1_exact_slope),
        "root_moment_Nm": float(abs(r1[1])),
        "exact_root_moment_Nm": float(v1_exact_moment),
        "root_moment_relative_error": _rel_error(abs(r1[1]), v1_exact_moment),
    }
    v1["passed"] = bool(
        max(
            v1["tip_displacement_relative_error"],
            v1["tip_slope_relative_error"],
            v1["root_moment_relative_error"],
        )
        <= tolerances["V1"]
    )
    ladder.append(
        {
            "test": "V1",
            "name": "Isolated tip-load cantilever",
            "computed": v1["tip_displacement_m"],
            "reference": v1["exact_tip_displacement_m"],
            "error": v1["tip_displacement_relative_error"],
            "error_kind": "relative",
            "tolerance": tolerances["V1"],
            "unit": "m",
            "passed": v1["passed"],
        }
    )

    # ------------------------------------------------------------------ V2
    # Same section/span, only self-weight, simply supported translations.
    v2_beam = Beam(L, E, mu, 0.0, 0.0, gravity=g, elements=linear_elements)
    nd = 2 * (linear_elements + 1)
    free = np.setdiff1d(np.arange(nd), [0, nd - 2])
    u2, _, _ = v2_beam.solve(D, t, free=free)
    mid_node = linear_elements // 2
    v2_num = float(u2[2 * mid_node])
    v2_exact = float(5.0 * q * L**4 / (384.0 * EI))
    v2_error = _rel_error(v2_num, v2_exact)
    v2 = {
        "midspan_displacement_m": v2_num,
        "exact_midspan_displacement_m": v2_exact,
        "relative_error": v2_error,
        "passed": bool(v2_error <= tolerances["V2"]),
    }
    ladder.append(
        {
            "test": "V2",
            "name": "Simply supported self-weight beam",
            "computed": v2_num,
            "reference": v2_exact,
            "error": v2_error,
            "error_kind": "relative",
            "tolerance": tolerances["V2"],
            "unit": "m",
            "passed": v2["passed"],
        }
    )

    # ------------------------------------------------------------------ V3
    mesh_convergence = []
    exact_f = float(tip_mass_frequencies(L, EI, mu, mt, 1)[0])
    for n in [2, 4, 8, 16]:
        beam = Beam(L, E, mu, mt, T, gravity=g, elements=n)
        l2 = _displacement_l2_error(beam, D, t, P, q)
        fe_f = float(beam.fe_modes(D, t, 1)[0][0])
        mesh_convergence.append(
            {
                "linear_elements": n,
                "h_m": float(L / n),
                "displacement_L2_error_m32": l2,
                "first_frequency_Hz": fe_f,
                "frequency_relative_error": _rel_error(fe_f, exact_f),
            }
        )
    disp_order = float(
        np.polyfit(
            np.log([r["h_m"] for r in mesh_convergence]),
            np.log([r["displacement_L2_error_m32"] for r in mesh_convergence]),
            1,
        )[0]
    )
    freq_order = float(
        np.log(mesh_convergence[0]["frequency_relative_error"] / mesh_convergence[2]["frequency_relative_error"])
        / np.log(4.0)
    )
    v3_error = abs(disp_order - 4.0) / 4.0
    v3 = {
        "displacement_order": disp_order,
        "frequency_order": freq_order,
        "relative_error_from_order_4": v3_error,
        "passed": bool(v3_error <= tolerances["V3_order_relative"]),
    }
    ladder.append(
        {
            "test": "V3",
            "name": "Displacement mesh convergence order",
            "computed": disp_order,
            "reference": 4.0,
            "error": v3_error,
            "error_kind": "relative",
            "tolerance": tolerances["V3_order_relative"],
            "unit": "1",
            "passed": v3["passed"],
        }
    )

    # ------------------------------------------------------------------ V4
    actual_beam = Beam(L, E, mu, mt, T, gravity=g, elements=linear_elements)
    fe_frequencies, _ = actual_beam.fe_modes(D, t, 5)
    exact_frequencies = tip_mass_frequencies(L, EI, mu, mt, 5)
    frequency_rows = []
    for mode, (f_fe, f_exact) in enumerate(zip(fe_frequencies, exact_frequencies), 1):
        frequency_rows.append(
            {
                "mode": mode,
                "fe_Hz": float(f_fe),
                "reference_Hz": float(f_exact),
                "relative_error": _rel_error(f_fe, f_exact),
                "percent_difference": 100.0 * _rel_error(f_fe, f_exact),
            }
        )
    v4_error = frequency_rows[0]["relative_error"]
    v4 = {
        "first_frequency_FE_Hz": frequency_rows[0]["fe_Hz"],
        "first_frequency_reference_Hz": frequency_rows[0]["reference_Hz"],
        "first_frequency_relative_error": v4_error,
        "modes": frequency_rows,
        "passed": bool(v4_error <= tolerances["V4"]),
    }
    ladder.append(
        {
            "test": "V4",
            "name": "First frequency with concentrated tip mass",
            "computed": v4["first_frequency_FE_Hz"],
            "reference": v4["first_frequency_reference_Hz"],
            "error": v4_error,
            "error_kind": "relative",
            "tolerance": tolerances["V4"],
            "unit": "Hz",
            "passed": v4["passed"],
        }
    )

    # ------------------------------------------------------------------ V5
    u5, reactions5, (K5, M5, F5, _, _) = actual_beam.solve(D, t)
    xnodes = np.linspace(0.0, L, linear_elements + 1)
    force_balance = float(reactions5[0] + np.sum(F5[::2]))
    moment_balance = float(reactions5[1] + xnodes @ F5[::2] + np.sum(F5[1::2]))
    eigvals = eigh(K5.toarray(), M5.toarray(), eigvals_only=True)
    rigid_modes = int(np.sum(np.abs(eigvals) < 1e-8 * np.max(np.abs(eigvals))))
    v5 = {
        "force_balance_residual_N": abs(force_balance),
        "moment_balance_residual_Nm": abs(moment_balance),
        "free_free_planar_rigid_modes": rigid_modes,
        "passed": bool(
            abs(force_balance) <= tolerances["V5_force_N"]
            and abs(moment_balance) <= tolerances["V5_moment_Nm"]
            and rigid_modes == 2
        ),
    }
    ladder.append(
        {
            "test": "V5",
            "name": "Global force balance",
            "computed": abs(force_balance),
            "reference": 0.0,
            "error": abs(force_balance),
            "error_kind": "absolute",
            "tolerance": tolerances["V5_force_N"],
            "unit": "N",
            "passed": v5["passed"],
        }
    )

    # ------------------------------------------------------------------ V6
    base_model = TubeModel(loads, study, nonlinear=True)
    gradient_rows = gradient_check(base_model, d0_mm)
    max_constraint_error = max(r["relative_error"] for r in gradient_rows if r["response"] != "mass")
    max_mass_error = max(r["relative_error"] for r in gradient_rows if r["response"] == "mass")
    taylor = _taylor_test(base_model, d0_mm)
    v6 = {
        "gradient_rows": gradient_rows,
        "maximum_constraint_gradient_error": float(max_constraint_error),
        "maximum_mass_gradient_error": float(max_mass_error),
        "taylor_test": taylor,
        "passed": bool(
            max_constraint_error <= tolerances["V6_constraint"]
            and max_mass_error <= tolerances["V6_mass"]
            and 1.8 <= taylor["measured_order"] <= 2.2
        ),
    }
    ladder.append(
        {
            "test": "V6",
            "name": "Maximum normalized derivative discrepancy",
            "computed": float(max_constraint_error),
            "reference": 0.0,
            "error": float(max_constraint_error),
            "error_kind": "absolute",
            "tolerance": tolerances["V6_constraint"],
            "unit": "1",
            "passed": v6["passed"],
        }
    )

    # ------------------------------------------------------------------ V7
    bvp = _continuum_follower_bvp(L, EI, T, mt, mu, g)
    continuum_tip = float(bvp.sol(L)[2])
    nonlinear_refinement = []
    records_80 = None
    theta_80 = None
    nl_80 = None
    for n in [40, 80, 160]:
        nl = Elastica(L, EI, -mt * g, q=q, elements=n, follower_force=T)
        theta, records = nl.solve(steps=4)
        tip = float(nl.coordinates(theta)[-1, 1])
        nonlinear_refinement.append(
            {
                "nonlinear_elements": n,
                "tip_displacement_m": tip,
                "continuum_tip_displacement_m": continuum_tip,
                "relative_error": _rel_error(tip, continuum_tip),
            }
        )
        if n == 80:
            records_80 = records
            theta_80 = theta
            nl_80 = nl

    # Consistent tangent check at the accepted 80-element physical case.
    rng = np.random.default_rng(20261125)
    direction = rng.normal(size=nl_80.n)
    direction /= np.linalg.norm(direction)
    theta_c = theta_80.astype(complex)
    theta_c[1:] += 1e-25j * direction
    complex_step = nl_80.residual(theta_c)[0].imag[1:] / 1e-25
    tangent = nl_80.residual(theta_80)[1][1:, 1:] @ direction
    tangent_error = float(np.linalg.norm(complex_step - tangent) / np.linalg.norm(tangent))

    accepted = next(r for r in nonlinear_refinement if r["nonlinear_elements"] == 80)
    v7 = {
        "continuum_tip_displacement_m": continuum_tip,
        "refinement": nonlinear_refinement,
        "accepted_80_element_tip_displacement_m": accepted["tip_displacement_m"],
        "accepted_80_element_relative_error": accepted["relative_error"],
        "tangent_relative_error": tangent_error,
        "newton_records_80": records_80,
        "passed": bool(accepted["relative_error"] <= tolerances["V7"] and tangent_error <= 1e-12),
    }
    ladder.append(
        {
            "test": "V7",
            "name": "N2 follower-load tip displacement",
            "computed": accepted["tip_displacement_m"],
            "reference": continuum_tip,
            "error": accepted["relative_error"],
            "error_kind": "relative",
            "tolerance": tolerances["V7"],
            "unit": "m",
            "passed": v7["passed"],
        }
    )

    # Actual arm analytical values used for the independent hand-calculation table.
    actual_tip = P * L**3 / (3.0 * EI) + q * L**4 / (8.0 * EI)
    actual_slope = P * L**2 / (2.0 * EI) + q * L**3 / (6.0 * EI)
    actual_root_moment = P * L + q * L**2 / 2.0
    root_stress = abs(actual_root_moment) * D / (2.0 * I)
    actual_u, _, _ = actual_beam.solve(D, t)

    hand_calculation_comparison = {
        "section_area_m2": float(A),
        "second_moment_m4": float(I),
        "tube_mass_kg": float(tube_mass),
        "tube_weight_N": float(tube_weight),
        "net_upward_tip_force_N": float(P),
        "distributed_self_weight_N_per_m": float(q),
        "linear_tip_displacement_formula_m": float(actual_tip),
        "linear_tip_displacement_FE_m": float(actual_u[-2]),
        "linear_tip_displacement_relative_error": _rel_error(actual_u[-2], actual_tip),
        "linear_tip_slope_formula_rad": float(actual_slope),
        "root_moment_formula_Nm": float(abs(actual_root_moment)),
        "root_bending_stress_formula_Pa": float(root_stress),
        "frequency_modes": frequency_rows,
    }

    status = "pass" if all(r["passed"] for r in ladder) else "fail"
    return {
        "status": status,
        "inputs": {
            "outer_diameter_mm": float(d0_mm[0]),
            "wall_mm": float(d0_mm[1]),
            "inner_diameter_mm": float(d0_mm[0] - 2.0 * d0_mm[1]),
            "length_m": L,
            "E_Pa": E,
            "density_kg_m3": rho,
            "thrust_N": T,
            "tip_mass_kg": mt,
            "gravity_m_s2": g,
            "net_upward_tip_force_N": float(P),
            "distributed_self_weight_N_per_m": float(q),
            "linear_elements": linear_elements,
            "nonlinear_verification_elements": 80,
        },
        "section": {
            "area_m2": float(A),
            "second_moment_m4": float(I),
            "EI_Nm2": float(EI),
            "mass_per_length_kg_m": float(mu),
        },
        "tolerances": tolerances,
        "ladder": ladder,
        "V1": v1,
        "V2": v2,
        "V3": v3,
        "V4": v4,
        "V5": v5,
        "V6": v6,
        "V7": v7,
        "mesh_convergence": mesh_convergence,
        "gradient_checks": gradient_rows,
        "taylor_test": taylor,
        "nonlinear_refinement": nonlinear_refinement,
        "hand_calculation_comparison": hand_calculation_comparison,
    }


def main():
    loads = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))
    study = json.loads((ROOT / "optimizer_inputs.json").read_text(encoding="utf-8"))
    result = run_base_case_verification(loads, study)
    print(json.dumps(result, indent=2, allow_nan=False))
    if result["status"] != "pass":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
