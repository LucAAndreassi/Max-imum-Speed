"""Production-model adjoint versus central-difference timing benchmark.

This benchmark uses the same TubeModel and discretization as the Project 4
optimizer.  It evaluates the full mass/constraint Jacobian at the optimizer's
preflight gradient-check point and compares analytic/adjoint derivatives with
central finite differences in the two physical design variables.

Run directly:
    python benchmark_production.py

Or import from run.py:
    from benchmark_production import run_production_benchmark
"""

from __future__ import annotations

import copy
import json
import time
from pathlib import Path

import numpy as np

from optimize import TubeModel, validate

ROOT = Path(__file__).resolve().parent


def _check_point(study: dict) -> np.ndarray:
    """Return the exact preflight point used by optimize.run_study()."""
    d = np.mean(
        [study["outer_diameter_bounds_mm"], study["wall_bounds_mm"]],
        axis=1,
    ).astype(float)
    d[0] = max(d[0], 30.0)
    d[0] = min(d[0], study["outer_diameter_bounds_mm"][1] * 0.95)
    return d


def _central_difference(loads: dict, study: dict, d: np.ndarray):
    """Compute the complete constraint Jacobian and mass gradient by CD."""
    model = TubeModel(loads, study, nonlinear=True)
    n_constraints = len(model.names)
    jac = np.zeros((n_constraints, 2), dtype=float)
    mass_grad = np.zeros(2, dtype=float)

    for j in range(2):
        h = 1.0e-4 * max(abs(float(d[j])), 1.0)
        delta = np.zeros(2, dtype=float)
        delta[j] = h
        plus = model.evaluate(d + delta, gradients=False)
        minus = model.evaluate(d - delta, gradients=False)
        jac[:, j] = (plus["g"] - minus["g"]) / (2.0 * h)
        mass_grad[j] = (plus["mass"] - minus["mass"]) / (2.0 * h)

    return jac, mass_grad


def _relative_error(a: np.ndarray, b: np.ndarray, floor: float) -> np.ndarray:
    denom = np.maximum(np.maximum(np.abs(a), np.abs(b)), floor)
    return np.abs(a - b) / denom


def run_production_benchmark(
    loads: dict,
    study: dict,
    repeats: int = 10,
) -> dict:
    """Benchmark the production analytic/adjoint Jacobian against CD.

    A fresh TubeModel is created for every timed trial so the model's single-
    entry evaluation cache cannot make either method artificially faster.
    Reported times are medians across ``repeats`` trials.
    """
    if repeats < 3:
        raise ValueError("Use at least 3 repeats for a stable median benchmark")

    loads = copy.deepcopy(loads)
    study = copy.deepcopy(study)
    validate(loads, study)

    d = _check_point(study)

    # Warm up imports, sparse solves, eigensolvers, and nonlinear code paths.
    warm_model = TubeModel(loads, study, nonlinear=True)
    warm_model.evaluate(d, gradients=True)
    _central_difference(loads, study, d)

    analytic_times = []
    fd_times = []
    analytic_result = None
    fd_jac = None
    fd_mass_grad = None

    for _ in range(repeats):
        model = TubeModel(loads, study, nonlinear=True)
        t0 = time.perf_counter()
        result = model.evaluate(d, gradients=True)
        analytic_times.append(time.perf_counter() - t0)
        analytic_result = result

        t0 = time.perf_counter()
        jac, mass_grad = _central_difference(loads, study, d)
        fd_times.append(time.perf_counter() - t0)
        fd_jac = jac
        fd_mass_grad = mass_grad

    if analytic_result is None or fd_jac is None or fd_mass_grad is None:
        raise RuntimeError("Benchmark did not produce derivative results")

    constraint_errors = _relative_error(
        analytic_result["jac"], fd_jac, floor=1.0e-7
    )
    mass_errors = _relative_error(
        analytic_result["mass_grad"], fd_mass_grad, floor=1.0e-12
    )

    analytic_ms = 1000.0 * float(np.median(analytic_times))
    fd_ms = 1000.0 * float(np.median(fd_times))
    speedup = fd_ms / analytic_ms

    passed = bool(
        np.max(constraint_errors) < 3.0e-4
        and np.max(mass_errors) < 1.0e-6
    )
    if not passed:
        raise RuntimeError(
            "Production benchmark derivative comparison failed: "
            f"constraint max={np.max(constraint_errors):.3e}, "
            f"mass max={np.max(mass_errors):.3e}"
        )

    return {
        "status": "pass",
        "repeats": int(repeats),
        "point_mm": {
            "outer_diameter": float(d[0]),
            "wall": float(d[1]),
        },
        "linear_elements": int(study["linear_elements"]),
        "nonlinear_elements": int(study["nonlinear_elements"]),
        "constraint_count": int(len(analytic_result["g"])),
        "design_variable_count": 2,
        "analytic_adjoint_median_ms": analytic_ms,
        "central_difference_median_ms": fd_ms,
        "speedup": float(speedup),
        "max_constraint_gradient_relative_error": float(
            np.max(constraint_errors)
        ),
        "max_mass_gradient_relative_error": float(np.max(mass_errors)),
        "analytic_constraint_jacobian": analytic_result["jac"].tolist(),
        "central_difference_constraint_jacobian": fd_jac.tolist(),
        "analytic_mass_gradient": analytic_result["mass_grad"].tolist(),
        "central_difference_mass_gradient": fd_mass_grad.tolist(),
        "constraint_names": list(TubeModel(loads, study, nonlinear=True).names),
        "variable_names": ["OD_mm", "wall_mm"],
    }


def main() -> None:
    loads = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))
    study = json.loads(
        (ROOT / "optimizer_inputs.json").read_text(encoding="utf-8")
    )
    result = run_production_benchmark(loads, study)

    print("Production Jacobian benchmark")
    print(
        "Point: "
        f"D={result['point_mm']['outer_diameter']:.3f} mm, "
        f"t={result['point_mm']['wall']:.3f} mm"
    )
    print(
        "Meshes: "
        f"{result['linear_elements']} linear / "
        f"{result['nonlinear_elements']} nonlinear elements"
    )
    print(f"Enabled constraints: {result['constraint_count']}")
    print(
        "Full analytic/adjoint Jacobian: "
        f"{result['analytic_adjoint_median_ms']:.3f} ms"
    )
    print(
        "Central-difference Jacobian: "
        f"{result['central_difference_median_ms']:.3f} ms"
    )
    print(f"Measured speedup: {result['speedup']:.3f}x")
    print(
        "Max constraint-gradient relative discrepancy: "
        f"{result['max_constraint_gradient_relative_error']:.3e}"
    )
    print(
        "Max mass-gradient relative discrepancy: "
        f"{result['max_mass_gradient_relative_error']:.3e}"
    )
    print("Status: PASS")


if __name__ == "__main__":
    main()
