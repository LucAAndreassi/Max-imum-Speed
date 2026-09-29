"""Diagnostics for DBF Project 2: penalty-induced ill-conditioning.

Produces D1-D4 numerical evidence and saves plots/data under project2_outputs/.

Physical propulsion and Mission-2 battery limits are enforced with SLSQP when
locating design optima.  For the long D3/D4 projected-gradient diagnostics, the
physical feasible set is linearized once at each local optimum and projected
with a cheap Dykstra half-space/box projection.  This preserves the intended
local-conditioning experiment without launching a nonlinear optimizer at every
first-order iteration.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from Project2Models import (
    BOUNDS,
    X0,
    augmented_lagrangian_objective,
    dbf_objective,
    penalized_objective,
    state,
    weight_constraint,
)

OUT = Path(__file__).resolve().parent / "project2_outputs"
OUT.mkdir(exist_ok=True)

LB = BOUNDS[:, 0]
UB = BOUNDS[:, 1]
SCALE = UB - LB

# The original D4 rho sweep adds many extra constrained solves and is not
# required for D1-D4.  Turn this on only when that additional plot is desired.
RUN_EXTRA_D4_SWEEP = True


def to_x(z):
    return LB + SCALE * np.asarray(z, dtype=float)


def to_z(x):
    return (np.asarray(x, dtype=float) - LB) / SCALE


def fd_gradient(fun, z, h=1e-6):
    """Forward finite-difference gradient in normalized design coordinates."""
    z = np.asarray(z, dtype=float)
    g = np.empty_like(z)
    f0 = fun(z)
    for i in range(len(z)):
        step = h * max(1.0, abs(z[i]))
        zp = z.copy()
        zp[i] += step
        g[i] = (fun(zp) - f0) / step
    return g


def fd_hessian(fun, z, h=2e-4):
    z = np.asarray(z, dtype=float)
    n = len(z)
    H = np.zeros((n, n))
    f0 = fun(z)

    for i in range(n):
        hi = h * max(1.0, abs(z[i]))
        zp = z.copy()
        zm = z.copy()
        zp[i] += hi
        zm[i] -= hi
        H[i, i] = (fun(zp) - 2.0 * f0 + fun(zm)) / hi**2

        for j in range(i + 1, n):
            hj = h * max(1.0, abs(z[j]))
            zpp = z.copy()
            zpm = z.copy()
            zmp = z.copy()
            zmm = z.copy()
            zpp[i] += hi
            zpp[j] += hj
            zpm[i] += hi
            zpm[j] -= hj
            zmp[i] -= hi
            zmp[j] += hj
            zmm[i] -= hi
            zmm[j] -= hj
            val = (fun(zpp) - fun(zpm) - fun(zmp) + fun(zmm)) / (4.0 * hi * hj)
            H[i, j] = H[j, i] = val

    return 0.5 * (H + H.T)


def positive_spectrum(H, floor=1e-10):
    ev = np.linalg.eigvalsh(H)
    pos = ev[ev > floor]
    return (ev, math.inf) if len(pos) < 2 else (ev, float(pos[-1] / pos[0]))


def jacobi_condition(H):
    d = np.diag(H).copy()
    if np.any(d <= 0):
        return math.inf, np.full_like(H, np.nan), np.array([])
    invsqrt = 1.0 / np.sqrt(d)
    Hs = (invsqrt[:, None] * H) * invsqrt[None, :]
    ev, kappa = positive_spectrum(Hs)
    return kappa, Hs, ev


# -----------------------------------------------------------------------------
# Physical feasibility
# -----------------------------------------------------------------------------

def physical_constraints_z(z):
    """Return dimensionless physical margins; every entry must be >= 0.

    Ordering:
      0 M2 straight propulsion margin / 50 ft/s
      1 M2 turn propulsion margin / 50 ft/s
      2 M3 straight propulsion margin / 50 ft/s
      3 M3 turn propulsion margin / 50 ft/s
      4 M2 battery margin / usable battery capacity

    This reads one mission state only, which is much faster than separately
    calling propulsion_margins() and battery_margins().
    """
    s = state(to_x(z))
    prop = np.array(
        [
            s["m2_speed_margin_straight"],
            s["m2_speed_margin_turn"],
            s["m3_speed_margin_straight"],
            s["m3_speed_margin_turn"],
        ],
        dtype=float,
    ) / 50.0

    batt = np.array(
        [s["m2_battery_margin"] / max(s["battery_capacity"], 1.0)],
        dtype=float,
    )
    return np.concatenate((prop, batt))


def optimize_penalty(rho, z0):
    """Find a fixed-weight-penalty optimum with physical limits enforced."""
    fun = lambda z: penalized_objective(to_x(z), rho)
    constraints = {"type": "ineq", "fun": physical_constraints_z}

    res = minimize(
        fun,
        np.clip(z0, 0.0, 1.0),
        method="SLSQP",
        bounds=[(0.0, 1.0)] * 6,
        constraints=constraints,
        options={"ftol": 1e-7, "maxiter": 2000, "disp": False},
    )
    return res, to_x(res.x)


def augmented_lagrangian(z0, rho0=5.0, max_outer=12):
    """Weight augmented Lagrangian with physical limits enforced by SLSQP."""
    z = np.clip(np.asarray(z0, float), 0.0, 1.0)
    lam = 0.0
    rho = float(rho0)
    history = []
    prev_violation = math.inf
    constraints = {"type": "ineq", "fun": physical_constraints_z}

    for outer in range(max_outer):
        fun = lambda zz: augmented_lagrangian_objective(to_x(zz), lam, rho)
        res = minimize(
            fun,
            z,
            method="SLSQP",
            bounds=[(0.0, 1.0)] * 6,
            constraints=constraints,
            options={"ftol": 1e-7, "maxiter": 2000, "disp": False},
        )

        z = res.x
        x = to_x(z)
        g = weight_constraint(x)
        violation = max(0.0, g)
        lam = max(0.0, lam + rho * g)
        history.append((outer, dbf_objective(x), g, rho, lam, res.nit))

        if violation < 1e-6:
            break
        if violation > 0.5 * prev_violation:
            rho *= 5.0
        prev_violation = violation

    return z, np.array(history, dtype=float), rho, lam


# -----------------------------------------------------------------------------
# Fast local feasible-set projection for D3/D4
# -----------------------------------------------------------------------------

def fd_constraint_jacobian(fun, z, h=2e-5):
    """Finite-difference Jacobian of vector constraints in z coordinates."""
    z = np.asarray(z, dtype=float)
    c0 = np.asarray(fun(z), dtype=float)
    J = np.zeros((len(c0), len(z)))

    for j in range(len(z)):
        hp = h * max(1.0, abs(z[j]))
        zp = z.copy()
        zm = z.copy()

        # Stay inside the normalized box when constructing the local model.
        zp[j] = min(1.0, z[j] + hp)
        zm[j] = max(0.0, z[j] - hp)
        denom = zp[j] - zm[j]

        if denom <= 1e-14:
            continue

        cp = np.asarray(fun(zp), dtype=float)
        cm = np.asarray(fun(zm), dtype=float)
        J[:, j] = (cp - cm) / denom

    return c0, J


class LocalFeasibleProjector:
    """Projection onto a local linearization of the physical feasible set.

    Around z_ref,
        c(z) ~= c_ref + J (z-z_ref) >= 0,
    which can be written as a set of linear half-spaces A z >= b.

    Projection uses Dykstra's algorithm over the normalized box and those
    half-spaces.  It is inexpensive: no SLSQP solve occurs inside D3/D4.
    """

    def __init__(self, z_ref, constraint_fun=physical_constraints_z, h=2e-5):
        self.z_ref = np.asarray(z_ref, dtype=float).copy()
        c_ref, J = fd_constraint_jacobian(constraint_fun, self.z_ref, h=h)
        self.c_ref = c_ref

        rows = []
        rhs = []
        for i in range(J.shape[0]):
            a = J[i].copy()
            nrm = float(np.linalg.norm(a))
            if nrm <= 1e-12:
                continue

            # c_ref + a^T(z-z_ref) >= 0
            # -> a^T z >= a^T z_ref - c_ref
            b = float(np.dot(a, self.z_ref) - c_ref[i])

            # Normalize each half-space to improve projection numerics.
            rows.append(a / nrm)
            rhs.append(b / nrm)

        self.A = np.asarray(rows, dtype=float) if rows else np.empty((0, len(z_ref)))
        self.b = np.asarray(rhs, dtype=float) if rhs else np.empty(0)

    def min_linear_margin(self, z):
        if len(self.b) == 0:
            return math.inf
        return float(np.min(self.A @ np.asarray(z, dtype=float) - self.b))

    def project(self, z_target, tol=1e-11, max_cycles=60):
        """Approximate Euclidean projection with Dykstra's algorithm."""
        x = np.asarray(z_target, dtype=float).copy()

        # Convex sets: one box plus one half-space per physical constraint.
        corrections = [np.zeros_like(x) for _ in range(1 + len(self.b))]

        for _ in range(max_cycles):
            x_start = x.copy()

            # Set 0: normalized design box [0,1]^n.
            y = x + corrections[0]
            x_new = np.clip(y, 0.0, 1.0)
            corrections[0] = y - x_new
            x = x_new

            # Remaining sets: local linearized physical half-spaces.
            for i, (a, b) in enumerate(zip(self.A, self.b), start=1):
                y = x + corrections[i]
                val = float(np.dot(a, y))
                if val < b:
                    # Rows are normalized, so ||a||^2 = 1.
                    x_new = y + (b - val) * a
                else:
                    x_new = y
                corrections[i] = y - x_new
                x = x_new

            if np.linalg.norm(x - x_start, ord=np.inf) <= tol:
                if (
                    np.all(x >= -1e-10)
                    and np.all(x <= 1.0 + 1e-10)
                    and self.min_linear_margin(x) >= -1e-9
                ):
                    break

        return np.clip(x, 0.0, 1.0)


def projected_gd(
    fun,
    z0,
    step,
    projector,
    tol=1e-6,
    maxit=200000,
    label=None,
    progress_every=5000,
):
    """Projected GD on the local linearized physical feasible set."""
    if step <= 0.0 or not np.isfinite(step):
        raise ValueError(f"Projected-GD step must be positive and finite; got {step}.")

    z = projector.project(z0)
    hist = []

    for k in range(maxit):
        f = fun(z)
        grad = fd_gradient(fun, z)

        z_trial = z - step * grad
        z_next = projector.project(z_trial)

        # Standard projected-gradient mapping for the chosen step length.
        pg = (z - z_next) / step
        pgn = float(np.linalg.norm(pg))
        hist.append((k, f, pgn))

        if pgn <= tol:
            break

        z = z_next

        if label and progress_every and (k + 1) % progress_every == 0:
            print(f"  {label}: iteration {k+1}/{maxit}, projected-grad={pgn:.3e}")

    return z, np.asarray(hist, dtype=float)



def main():
    print("Running smooth DBF Project 2 diagnostics...")

    # ------------------------------------------------------------------
    # Fixed-penalty optima used by D1-D3.
    # ------------------------------------------------------------------
    rhos = [0.1, 1.0, 10.0, 100.0, 500.0, 1000.0, 10000.0]
    zstart = to_z(X0)
    rows = []
    optima = {}

    for rho in rhos:
        res, xstar = optimize_penalty(rho, zstart)
        zstar = res.x
        zstart = zstar

        # Smooth active branch used for the local Hessian conditioning study.
        fun_h = lambda z, rr=rho: (
            dbf_objective(to_x(z))
            + 0.5 * rr * weight_constraint(to_x(z)) ** 2
        )
        H = fd_hessian(fun_h, zstar, h=5e-5)
        ev, kappa = positive_spectrum(H)
        kj, Hs, evj = jacobi_condition(H)

        s = state(xstar)
        g = weight_constraint(xstar)
        min_phys = float(np.min(physical_constraints_z(zstar)))

        rows.append(
            {
                "rho": rho,
                "success": res.success,
                "iterations": res.nit,
                "objective": res.fun,
                "base_objective": dbf_objective(xstar),
                "constraint_g": g,
                "m2_weight": s["m2_takeoff_weight"],
                "kappa": kappa,
                "kappa_jacobi": kj,
                "lambda_min": float(ev[ev > 1e-10][0]) if np.any(ev > 1e-10) else np.nan,
                "lambda_max": float(ev[-1]),
                "V2": xstar[4],
                "V3": xstar[5],
                "min_physical_margin": min_phys,
            }
        )
        optima[rho] = (xstar, zstar, H, ev)

        print(
            f"rho={rho:8g} f={res.fun: .6f} g={g:+.3e} "
            f"kappa={kappa:.3e} jacobi={kj:.3e} nit={res.nit} "
            f"V2={xstar[4]:.2f} V3={xstar[5]:.2f} min_phys={min_phys:+.2e}"
        )

    with (OUT / "conditioning_vs_rho.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ------------------------------------------------------------------
    # D1: Hessian spectrum.
    # ------------------------------------------------------------------
    rho_spec = rhos[-1]
    _, _, _, ev = optima[rho_spec]
    plt.figure(figsize=(7, 4.5))
    plt.semilogy(np.arange(1, len(ev) + 1), np.maximum(np.abs(ev), 1e-14), "o-")
    plt.xlabel("Eigenvalue index")
    plt.ylabel("|Hessian eigenvalue|")
    plt.title(f"D1: Penalized DBF Hessian spectrum ($\\rho={rho_spec:g}$)")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUT / "D1_hessian_spectrum.png", dpi=180)
    plt.close()

    # ------------------------------------------------------------------
    # D2: intrinsic-conditioning test.
    # ------------------------------------------------------------------
    rr = np.array([r["rho"] for r in rows])
    kk = np.array([r["kappa"] for r in rows])
    kj = np.array([r["kappa_jacobi"] for r in rows])

    plt.figure(figsize=(7, 4.5))
    plt.loglog(rr, kk, "o-", label="Original Hessian")
    plt.loglog(rr, kj, "s--", label="After Jacobi rescaling")
    plt.xlabel(r"Penalty weight $\rho$")
    plt.ylabel(r"Condition number $\kappa$")
    plt.title("D2: Intrinsic ill-conditioning test")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "D2_kappa_vs_rho.png", dpi=180)
    plt.close()

    # ------------------------------------------------------------------
    # D3: same first-order diagnostic at several penalty weights.
    # Physical feasibility is represented by a LOCAL linearization around each
    # constrained optimum.  This is intentionally a local conditioning test.
    # ------------------------------------------------------------------
    print("\nRunning D3 local projected-gradient diagnostics...")
    gd_rhos = [500, 1000, 10000]
    plt.figure(figsize=(7, 4.5))
    gd_summary = []

    z_common = to_z(np.array([7.0, 6.0, 10.0, 2.0, 100.0, 90.0]))
    base_center = to_z(np.array([6.0, 6.0, 10.0, 2.0, 100.0, 90.0]))
    direction = z_common - base_center
    d = direction / max(np.linalg.norm(direction), 1e-12)

    for rho in gd_rhos:
        xstar, zstar, Hstar, evstar = optima[rho]
        positive = evstar[evstar > 1e-10]
        L, mu = positive[-1], positive[0]
        step = 2.0 / (L + mu)

        fun = lambda z, rr=rho: penalized_objective(to_x(z), rr)
        projector = LocalFeasibleProjector(zstar)
        z_gd0 = projector.project(zstar + 0.08 * d)

        zend, hist = projected_gd(
            fun,
            z_gd0,
            step,
            projector,
            tol=1e-5,
            maxit=15000,
            label=f"D3 rho={rho:g}",
        )

        f_slsqp = fun(zstar)
        f_ref = min(f_slsqp, float(np.min(hist[:, 1])))
        if f_ref < f_slsqp - 1e-8:
            print(
                f"NOTE D3 rho={rho:g}: local projected GD improved the "
                f"SLSQP reference by {f_slsqp-f_ref:.3e}."
            )

        gap = np.maximum(hist[:, 1] - f_ref, 1e-14)
        actual_margin = float(np.min(physical_constraints_z(zend)))

        plt.semilogy(hist[:, 0], gap, label=rf"$\rho={rho:g}$")
        gd_summary.append(
            (rho, len(hist), hist[-1, 2], fun(zend) - f_ref, actual_margin)
        )

        print(
            f"  D3 rho={rho:g}: iterations={len(hist)}, "
            f"final_pg={hist[-1,2]:.3e}, actual_min_phys={actual_margin:+.3e}"
        )

    plt.xlabel("Iteration")
    plt.ylabel("Objective Gap")
    plt.title("D3: Objective Gap vs Iteration for Various Penalty Weights")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "D3_gradient_descent_convergence.png", dpi=180)
    plt.close()

    with (OUT / "gradient_descent_summary.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "rho",
                "iterations",
                "final_projected_grad",
                "final_objective_gap",
                "final_actual_min_physical_margin",
            ]
        )
        w.writerows(gd_summary)

    # ------------------------------------------------------------------
    # D4: augmented Lagrangian solution and local conditioning.
    # ------------------------------------------------------------------
    print("\nRunning D4 augmented-Lagrangian diagnostics...")
    z_seed = to_z(np.array([7.0, 6.0, 10.0, 2.0, 100.0, 90.0]))
    z_al, hist_al, rho_final, lam_final = augmented_lagrangian(z_seed, rho0=5.0)
    x_al = to_x(z_al)
    g_al = weight_constraint(x_al)
    s_al = state(x_al)

    print("\nFinal augmented-Lagrangian aircraft metrics")
    print(f"Empty aircraft weight: {s_al['empty_weight']:.4f} lb")
    print(f"M2 gross weight:       {s_al['m2_takeoff_weight']:.4f} lb")
    print(f"M2 five-lap time:      {s_al['m2_time']:.4f} s")
    print(f"M3 continuous laps:    {s_al['laps3']:.4f}")

    x_pen, z_pen, H_pen, ev_pen = optima[10000.0]

    fun_pen_local = lambda z: (
        dbf_objective(to_x(z))
        + 0.5 * 10000.0 * weight_constraint(to_x(z)) ** 2
    )
    fun_al_local = lambda z: (
        dbf_objective(to_x(z))
        + 0.5
        * rho_final
        * (weight_constraint(to_x(z)) + lam_final / rho_final) ** 2
        - 0.5 * lam_final**2 / rho_final
    )

    H_al = fd_hessian(fun_al_local, z_al, h=5e-5)
    ev_al, k_al = positive_spectrum(H_al)
    k_pen = positive_spectrum(H_pen)[1]

    with (OUT / "augmented_lagrangian_history.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "outer_iteration",
                "base_objective",
                "constraint_g",
                "rho",
                "lambda",
                "inner_iterations",
            ]
        )
        w.writerows(hist_al)

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(hist_al[:, 0], np.maximum(np.abs(hist_al[:, 2]), 1e-12), "o-")
    plt.xlabel("Augmented-Lagrangian outer iteration")
    plt.ylabel("|weight-constraint residual|")
    plt.title("D4: Augmented Lagrangian enforces the weight constraint")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUT / "D4_augmented_lagrangian_constraint.png", dpi=180)
    plt.close()

    # D4 before/after first-order comparison using separate local linearized
    # physical feasible sets about the fixed-penalty and AL solutions.
    pen_pos = ev_pen[ev_pen > 1e-10]
    al_pos = ev_al[ev_al > 1e-10]
    step_pen = 2.0 / (pen_pos[-1] + pen_pos[0])
    step_al = 2.0 / (al_pos[-1] + al_pos[0])

    projector_pen = LocalFeasibleProjector(z_pen)
    projector_al = LocalFeasibleProjector(z_al)

    # The fixed-penalty Hessian is extremely stiff at rho=10,000.  Keep the
    # D4 perturbation inside the genuinely local region so the first-order
    # comparison measures conditioning rather than large-step nonlinearity.
    # The same perturbation magnitude/direction is used for both formulations.
    d4_perturbation = 0.005
    z_pen0 = projector_pen.project(z_pen + d4_perturbation * d)
    z_al0 = projector_al.project(z_al + d4_perturbation * d)

    maxit_d4 = 20000
    tol_d4 = 1e-5

    # Use exactly the same projected-gradient algorithm for both formulations.
    zend_pen, hist_pen = projected_gd(
        fun_pen_local,
        z_pen0,
        step_pen,
        projector_pen,
        tol=tol_d4,
        maxit=maxit_d4,
        label="D4 fixed penalty",
    )
    zend_al, hist_al_gd = projected_gd(
        fun_al_local,
        z_al0,
        step_al,
        projector_al,
        tol=tol_d4,
        maxit=maxit_d4,
        label="D4 augmented Lagrangian",
    )

    # Plot the best objective attained up to each iteration.  This is a
    # standard convergence-performance view and avoids the misleading V-shape
    # that appears if a fixed step overshoots slightly after already reaching
    # its best local objective.  Nothing is removed from the optimizer run; the
    # curve simply answers: "what is the best gap achieved by iteration k?"
    best_pen_obj = np.minimum.accumulate(hist_pen[:, 1])
    best_al_obj = np.minimum.accumulate(hist_al_gd[:, 1])

    fref_pen = min(float(fun_pen_local(z_pen)), float(best_pen_obj[-1]))
    fref_al = min(float(fun_al_local(z_al)), float(best_al_obj[-1]))

    gap_pen = np.maximum(best_pen_obj - fref_pen, 1e-14)
    gap_al = np.maximum(best_al_obj - fref_al, 1e-14)

    # Normalize each curve by its own initial gap so the figure compares
    # convergence rate instead of the very different absolute curvature scales.
    norm_gap_pen = gap_pen / max(gap_pen[0], 1e-14)
    norm_gap_al = gap_al / max(gap_al[0], 1e-14)

    plt.figure(figsize=(7, 4.5))
    plt.semilogy(
        hist_pen[:, 0],
        norm_gap_pen,
        label="Fixed penalty (rho=10,000)",
    )
    plt.semilogy(
        hist_al_gd[:, 0],
        norm_gap_al,
        label=f"Augmented Lagrangian ($\\rho={rho_final:g}$)",
    )
    plt.xlabel("Projected-gradient iteration")
    plt.ylabel("Best-so-far normalized objective gap")
    plt.title("D4: Before/after convergence with the same first-order method")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "D4_before_after_convergence.png", dpi=180)
    plt.close()

    # ------------------------------------------------------------------
    # D4 local geometry near the weight boundary.
    # ------------------------------------------------------------------
    rho_fix = 10000.0
    rho_al = rho_final
    lam = lam_final
    g_activation = -lam / rho_al
    g = np.linspace(-0.025, 0.010, 1400)

    P_fix = 0.5 * rho_fix * np.maximum(0.0, g) ** 2
    P_al = (
        0.5 * rho_al * np.maximum(0.0, g + lam / rho_al) ** 2
        - 0.5 * lam**2 / rho_al
    )

    plt.figure(figsize=(7, 4.5))
    plt.plot(g, P_fix, linewidth=2.0, label=fr"Fixed penalty ($\rho={rho_fix:g}$)")
    plt.plot(
        g,
        P_al,
        linewidth=2.0,
        label=fr"Augmented Lagrangian ($\rho={rho_al:g},\ \lambda={lam:.3f}$)",
    )
    plt.axvline(0.0, linestyle="--", linewidth=1.4, label=r"Weight boundary $g(x)=0$")
    plt.axvline(
        g_activation,
        linestyle=":",
        linewidth=1.4,
        label=fr"AL activation $g=-\lambda/\rho={g_activation:.4f}$",
    )
    plt.axhline(0.0, linewidth=1.0)
    plt.yscale("symlog", linthresh=1e-4)
    plt.xlabel(r"Weight-constraint residual $g(x)=W_{M2}-20$")
    plt.ylabel("Local constraint contribution to objective")
    plt.title("D4: Local geometry near the weight-constraint boundary")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "D4_local_boundary_geometry.png", dpi=180)
    plt.close()

    # ------------------------------------------------------------------
    # Optional extra D4 rho sweep.  Disabled by default because it requires
    # seven additional fixed-penalty solves and seven additional AL solves.
    # ------------------------------------------------------------------
    if RUN_EXTRA_D4_SWEEP:
        print("\nRunning optional D4 rho sweep...")
        rho_vals = np.array([0.1, 1, 5, 10, 100, 1000, 10000.0])
        z_cmp = to_z(np.array([7.0, 6.0, 10.0, 2.0, 100.0, 90.0]))
        g_fix_vals = []
        g_al_vals = []

        for rho in rho_vals:
            _, x_fix = optimize_penalty(rho, z_cmp)
            z_a, _, _, _ = augmented_lagrangian(z_cmp, rho0=rho)
            g_fix_vals.append(weight_constraint(x_fix))
            g_al_vals.append(weight_constraint(to_x(z_a)))

        v_fix = np.maximum(g_fix_vals, 1e-8)
        v_al = np.maximum(g_al_vals, 1e-8)

        plt.figure(figsize=(7, 4.5))
        plt.loglog(rho_vals, v_fix, "o-", label="Fixed quadratic penalty")
        plt.loglog(rho_vals, v_al, "s--", label="Augmented Lagrangian")
        plt.axhline(1e-6, linestyle=":", label=r"$10^{-6}$ feasibility tolerance")
        plt.xlabel(r"Penalty parameter $\rho$")
        plt.ylabel("Weight-constraint violation [lb]")
        plt.title("D4: Weight-constraint violation comparison")
        plt.grid(True, which="both", alpha=0.25)
        plt.legend()
        plt.tight_layout()
        plt.savefig(OUT / "D4_weight_violation_fixed_vs_AL.png", dpi=180)
        plt.close()

    d4_rows = [
        [
            "fixed_penalty",
            10000.0,
            k_pen,
            weight_constraint(x_pen),
            len(hist_pen),
            hist_pen[-1, 2],
            float(gap_pen[-1]),
            step_pen,
            float(np.min(physical_constraints_z(zend_pen))),
        ],
        [
            "augmented_lagrangian",
            rho_final,
            k_al,
            g_al,
            len(hist_al_gd),
            hist_al_gd[-1, 2],
            float(gap_al[-1]),
            step_al,
            float(np.min(physical_constraints_z(zend_al))),
        ],
    ]

    with (OUT / "d4_before_after_summary.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "formulation",
                "rho",
                "kappa",
                "constraint_g",
                "iterations",
                "final_projected_grad",
                "final_objective_gap",
                "step",
                "final_actual_min_physical_margin",
            ]
        )
        w.writerows(d4_rows)

    print("\nD4 comparison")
    print(
        f"Fixed penalty rho=10000: g={weight_constraint(x_pen):+.3e}, "
        f"kappa={k_pen:.3e}, iters={len(hist_pen)}, "
        f"pg={hist_pen[-1,2]:.3e}, gap={gap_pen[-1]:.3e}"
    )
    print(
        f"Augmented Lagrangian:      g={g_al:+.3e}, rho_final={rho_final:g}, "
        f"kappa={k_al:.3e}, iters={len(hist_al_gd)}, "
        f"pg={hist_al_gd[-1,2]:.3e}, gap={gap_al[-1]:.3e}"
    )
    print("AL design:", np.array2string(x_al, precision=5))
    print("Outputs written to", OUT)


if __name__ == "__main__":
    main()
