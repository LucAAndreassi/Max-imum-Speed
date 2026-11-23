"""Multistart drone-arm mass optimization connected to the existing FEM.

Run python optimize_arm.py --update-report. Physical inputs, material data,
constraints and FE adjoints come from optimize.TubeModel and the JSON cards.
"""
from dataclasses import dataclass
from functools import lru_cache
from math import ceil, floor

import numpy as np
from scipy.optimize import brentq, minimize, nnls

BOUND_NAMES = ("OD_lower_bound", "OD_upper_bound", "wall_lower_bound", "wall_upper_bound")
BOUND_JACOBIAN = np.array([[-1., 0.], [1., 0.], [0., -1.], [0., 1.]])
TOLERANCES = dict(active=2e-5, primal=2e-6, dual=1e-10,
                  stationarity=1e-5, complementarity=2e-6, slsqp_ftol=1e-11)
SLSQP_MAX_ITERATIONS = 150


@dataclass
class Evaluation:
    point: dict
    objective: float
    gradient: np.ndarray | None
    residuals: np.ndarray  # normalized g <= 0, physical rows then bounds
    jacobian: np.ndarray | None
    names: tuple


def kkt_check(ev, mass_scale, engineering_scales):
    """Reconstruct nonnegative multipliers on the numerical active set.

    F = mass/mass_scale and x = dimensions/reference dimensions. FE
    constraints are dimensionless utilization residuals. Bound multipliers
    are also converted to kg/mm for engineering interpretation.
    """
    active = np.abs(ev.residuals) <= TOLERANCES["active"]
    multipliers = np.zeros_like(ev.residuals)
    active_matrix = ev.jacobian[active]
    rank = 0
    if np.any(active):
        multipliers[active], _ = nnls(active_matrix.T, -ev.gradient)
        rank = int(np.linalg.matrix_rank(active_matrix))
    lagrangian_gradient = ev.gradient + ev.jacobian.T @ multipliers
    checks = dict(primal=ev.point["max_scaled_violation"],
                  dual=float(max(0, -multipliers.min())),
                  stationarity=float(np.linalg.norm(lagrangian_gradient, ord=np.inf)),
                  complementarity=float(np.max(np.abs(multipliers * ev.residuals))))
    rows = []
    for i, name in enumerate(ev.names):
        scale, unit = engineering_scales[i]
        rows.append(dict(name=name, active=bool(active[i]),
                         normalized_residual=float(ev.residuals[i]),
                         scaled_multiplier=float(multipliers[i]),
                         physical_multiplier=float(multipliers[i] * mass_scale / scale), unit=unit))
    return dict(**checks, passed=all(value <= TOLERANCES[k] for k, value in checks.items()),
                active_rank=rank, licq=bool(rank == np.count_nonzero(active)),
                lagrangian_gradient=lagrangian_gradient.tolist(), rows=rows)


def fixed_wall_reference(model, point_at):
    """Screen nondecreasing nominal ODs at the configured reference wall.

    Use the project's OD rounding increment (mm). Locate a frequency
    crossing only if a minimum frequency was specified. A reference must
    satisfy every enabled physical constraint.
    """
    study = model.study
    step = study["rounding_mm"]["outer_diameter"]
    od, wall = study["reference_design_mm"]
    low, high = max(od, study["outer_diameter_bounds_mm"][0]), study["outer_diameter_bounds_mm"][1]
    result = dict(grid_mm=step, point=None, frequency_threshold=None, threshold_kind=None,
                  previous_grid_point=None, reason=None)
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Reference OD rounding increment must be positive and finite")
    if not study["wall_bounds_mm"][0] <= wall <= study["wall_bounds_mm"][1] or low > high:
        result["reason"] = "Sizing bounds exclude nondecreasing OD at the reference wall."
        return result
    if "bending_frequency" in model.names:
        def frequency_residual(d):
            return point_at(d)["constraint_residuals"]["bending_frequency"]
        if frequency_residual(low) <= 0:
            result.update(frequency_threshold=point_at(low), threshold_kind="already_meets")
        elif frequency_residual(high) <= 0:
            crossing = brentq(frequency_residual, low, high, xtol=1e-10)
            result.update(frequency_threshold=point_at(crossing), threshold_kind="crossing")
        else:
            result["reason"] = "Maximum OD at the reference wall misses the specified frequency."
            return result
    for i in range(ceil(low / step - 1e-12), floor(high / step + 1e-12) + 1):
        candidate = point_at(i * step)
        if max(candidate["constraint_residuals"].values()) <= TOLERANCES["slsqp_ftol"]:
            result["point"] = candidate
            return result
        result["previous_grid_point"] = candidate
    result["reason"] = "No nominal OD meets all specified limits at the reference wall."
    return result


def solve(model, *, progress=True, include_reference=True):
    """Use current FE responses/adjoints with scaled variables and multistart.

    Only solver-successful, feasible, KKT-passing runs can be selected. The
    result extends the existing report/export schema with full diagnostics.
    """
    from optimize import validate
    validate(model.loads, model.study)
    study = model.study
    dimensions = np.asarray(study["reference_design_mm"], float)
    box = np.array([study["outer_diameter_bounds_mm"], study["wall_bounds_mm"]], float)
    if dimensions.shape != (2,) or not np.isfinite(dimensions).all() or np.any(dimensions <= 0):
        raise ValueError("Reference dimensions must be positive and finite")
    lower, upper = (box / dimensions[:, None]).T
    mass_scale = float(model.mass(dimensions))
    scales = [(1., "kg/unit residual")] * len(model.names)
    scales += [(d, "kg/mm") for d in np.repeat(dimensions, 2)]

    @lru_cache(maxsize=8)
    def evaluate_at(x_d, x_t, gradients=True):
        d = dimensions * [x_d, x_t]
        a = model.evaluate(d, gradients=gradients)
        bounds = [lower[0]-x_d, x_d-upper[0], lower[1]-x_t, x_t-upper[1]]
        residuals = np.r_[a["g"], bounds]
        if not np.isfinite(np.r_[a["mass"], residuals]).all():
            raise ValueError("The FEM returned a nonfinite response")
        point = dict(outer_diameter_mm=float(d[0]), wall_mm=float(d[1]),
                     mass_kg=a["mass"], constraint_residuals=dict(zip(model.names, a["g"].tolist())),
                     max_scaled_violation=float(max(0., residuals.max())))
        gradient = a["mass_grad"] * dimensions / mass_scale if gradients else None
        jacobian = np.vstack((a["jac"] * dimensions, BOUND_JACOBIAN)) if gradients else None
        if gradients and not (np.isfinite(gradient).all() and np.isfinite(jacobian).all()):
            raise ValueError("The FEM returned a nonfinite derivative")
        return Evaluation(point, a["mass"] / mass_scale, gradient, residuals, jacobian,
                          tuple(model.names) + BOUND_NAMES)

    reference = fixed_wall_reference(
        model, lambda od: evaluate_at(od/dimensions[0], 1., False).point) if include_reference else None
    starts = [(f"Configured {i+1}", np.asarray(s, float)/dimensions)
              for i, s in enumerate(study["starting_designs_mm"])]
    starts += [("Small tube", lower), ("Large tube", upper)]
    if reference and reference["point"]:
        starts.append(("Fixed-wall reference", [reference["point"]["outer_diameter_mm"]/dimensions[0], 1.]))
    runs, records, seen = [], [], set()
    for label, start in starts:
        start = np.clip(start, lower, upper)
        if tuple(start) in seen:
            continue
        seen.add(tuple(start))
        history = []
        run = dict(label=label, start=len(runs), success=False, accepted=False, status=None,
                   message="", iterations=0, design=None, kkt=None, history=history)
        runs.append(run)

        def record(x):
            ev = evaluate_at(*x)
            history.append(dict(iteration=len(history), **ev.point))
            records.append(dict(start=run["start"], iteration=len(history)-1,
                                outer_diameter_mm=ev.point["outer_diameter_mm"], wall_mm=ev.point["wall_mm"],
                                mass_g=ev.point["mass_kg"]*1000, max_constraint=float(np.max(ev.residuals))))

        try:
            record(start)
            opt = minimize(lambda x: evaluate_at(*x).objective, start, method="SLSQP",
                           jac=lambda x: evaluate_at(*x).gradient,
                           bounds=list(zip(lower, upper)),
                           constraints=[dict(type="ineq", fun=lambda x: -evaluate_at(*x).residuals[:-4],
                                             jac=lambda x: -evaluate_at(*x).jacobian[:-4])],
                           callback=record, options=dict(ftol=TOLERANCES["slsqp_ftol"],
                                                         maxiter=SLSQP_MAX_ITERATIONS))
            ev = evaluate_at(*opt.x)
            kkt = kkt_check(ev, mass_scale, scales)
            run.update(design=model.summary(dimensions*opt.x), iterations=int(opt.nit),
                       success=bool(opt.success), status=int(opt.status), message=str(opt.message),
                       kkt=kkt, accepted=bool(opt.success and kkt["passed"]))
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError) as exc:
            run.update(iterations=max(0, len(history)-1), message=str(exc))
        if progress:
            print(f"{'Linear + follower' if model.nonlinear else 'Linear'} {label}: "
                  f"{run['iterations']} iterations; accepted={run['accepted']}; {run['message']}", flush=True)

    accepted = [i for i, run in enumerate(runs) if run["accepted"]]
    if not accepted:
        raise RuntimeError("No converged feasible KKT-passing design: " +
                           "; ".join(f"{r['label']}: {r['message']}" for r in runs))
    lightest = min(runs[i]["design"]["tube_mass_g"] for i in accepted)
    tied = [i for i in accepted if runs[i]["design"]["tube_mass_g"]-lightest <= mass_scale*1000*TOLERANCES["slsqp_ftol"]]
    chosen = min(tied, key=lambda i: (runs[i]["kkt"]["primal"], runs[i]["design"]["tube_mass_g"]))
    best = runs[chosen]
    # Convert bound multipliers back to the report's mm convention.
    rows = best["kkt"]["rows"]
    active = [dict(name=r["name"], residual=r["normalized_residual"]*scales[i][0],
                   scaled_multiplier=r["physical_multiplier"]/mass_scale)
              for i, r in enumerate(rows) if r["active"]]
    legacy_kkt = dict(active=active, stationarity_inf_norm=best["kkt"]["stationarity"],
                      complementarity_inf_norm=best["kkt"]["complementarity"],
                      maximum_violation=best["kkt"]["primal"], passed=best["kkt"]["passed"],
                      diagnostics=best["kkt"],
                      note="KKT norms use mass/reference mass and x=d/reference dimensions. "
                           "Active multipliers use mm bounds; diagnostics include kg/mm. Local evidence only.")
    return dict(optimum=best["design"], starts=runs, history=records, kkt=legacy_kkt,
                best_run=chosen, fixed_wall_reference=reference,
                settings=dict(engine="optimize_arm", dimension_scales_mm=dimensions.tolist(),
                              mass_scale_kg=mass_scale, tolerances=TOLERANCES.copy(),
                              max_iterations=SLSQP_MAX_ITERATIONS))


def check_shadow_prices(model, result):
    """Validate active bound/frequency mass prices with independent re-solves."""
    import copy

    rows = result["kkt"]["diagnostics"]["rows"]
    trials = []
    for variable, key in (("OD", "outer_diameter_bounds_mm"), ("wall", "wall_bounds_mm")):
        low, high = model.study[key]
        step = min(.001, (high-low)/2)  # mm numerical perturbation, not manufacturing tolerance
        for name, bounds in ((variable+"_lower_bound", [low+step, high]),
                             (variable+"_upper_bound", [low, high-step])):
            row = next(r for r in rows if r["name"] == name)
            if row["active"]:
                study = copy.deepcopy(model.study)
                study[key] = bounds
                trials.append((row, step, study, "mm", row["physical_multiplier"]*step))
    for row in rows:
        if row["name"] == "bending_frequency" and row["active"]:
            study = copy.deepcopy(model.study)
            target = model.minimum_frequency
            step = .01
            study["minimum_bending_frequency_rotor_multiple"] = None
            study["minimum_bending_frequency_Hz"] = target + step
            trials.append((row, step, study, "Hz", row["physical_multiplier"]*step/target))
    checks = []
    for row, step, study, unit, predicted in trials:
        change, reason = None, None
        try:
            trial = solve(type(model)(model.loads, study, model.nonlinear), progress=False, include_reference=False)
            change = (trial["optimum"]["tube_mass_g"]-result["optimum"]["tube_mass_g"])/1000
        except (ValueError, RuntimeError) as exc:
            reason = str(exc)
        checks.append(dict(name=row["name"], step=step, unit=unit, predicted_kg=predicted,
                           change_kg=change, unavailable_reason=reason))
    return checks


def main(argv=None):
    # Both command names share one validated FEM/report/export pipeline.
    from optimize import main as run
    return run(argv)


if __name__ == "__main__":
    main()
