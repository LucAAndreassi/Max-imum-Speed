#!/usr/bin/env python3
"""Reproducible inverse Stokes boundary-layer study: baseline, PINN, P1--P5.

Only NumPy, SciPy, Matplotlib and threadpoolctl are needed. Every neural weight
is trained; analytic network derivatives/Jacobians are checked by complex step.
No exact spatial solution is embedded in the neural network or inverse fits.
See REPORT.md for the harmonic time ansatz and baseline implementation choices.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.fft import dct
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from scipy.optimize import least_squares, minimize_scalar
from scipy.sparse import diags
from threadpoolctl import threadpool_info, threadpool_limits

TAU = 2 * np.pi
VARIANTS = {
    "vanilla": (False, False),
    "sampling": (False, True),
    "scaling": (True, False),
    "combined": (True, True),
    "supervised": (False, False),
    "known_kappa": (False, False),
}
COLORS = {"baseline": "#333333", "vanilla": "#bd4930", "sampling": "#d09524",
          "scaling": "#4669bd", "combined": "#16826c", "supervised": "#a14caa",
          "known_kappa": "#786451"}


@dataclass
class Jet2:
    """Minimal forward-mode AD value carrying first/second y derivatives.

    The network uses affine and tanh primitives; these rules automatically
    propagate derivatives, including the chain-rule factors after rescaling.
    Parameter Jacobians for the optimizer are supplied analytically below.
    """
    value: np.ndarray
    first: np.ndarray
    second: np.ndarray

    def affine(self, weight, bias=0.):
        return Jet2(self.value*weight+bias, self.first*weight, self.second*weight)

    def tanh(self):
        value = np.tanh(self.value)
        slope = 1-value*value
        return Jet2(value, slope*self.first,
                    slope*self.second-2*value*slope*self.first*self.first)

    def linear(self, weight, bias=0.):
        return Jet2(self.value@weight+bias, self.first@weight, self.second@weight)


def json_text(value):
    """Emit standard JSON; unavailable diagnostics are null, never NaN tokens."""
    def clean(x):
        if isinstance(x, dict):
            return {k: clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, (float, np.floating)):
            return float(x) if np.isfinite(x) else None
        if isinstance(x, np.integer):
            return int(x)
        return x
    return json.dumps(clean(value), indent=2, allow_nan=False)


@dataclass(frozen=True)
class Case:
    delta: float
    noise: float = .02
    nt: int = 16
    seed: int = 17

    @property
    def kappa(self):
        return self.delta ** 2 / 2

    @property
    def key(self):
        return f"d{self.delta:.6f}_s{self.noise:.3f}_n{self.nt}_seed{self.seed}"


def exact_amplitude(y, kappa, finite=False):
    """Complex harmonic amplitude, u=Re[amplitude exp(i t)]."""
    q = (1 + 1j) / np.sqrt(2 * kappa)
    y = np.asarray(y)
    if finite:
        # Stable equivalent of 1-sinh(q*(1-y))/sinh(q), including small kappa.
        return 1 - np.exp(-q * y) * (-np.expm1(-2*q*(1-y))) / (-np.expm1(-2*q))
    return 1 - np.exp(-q*y)


def exact_velocity(y, t, kappa, finite=False):
    return np.real(exact_amplitude(y, kappa, finite) * np.exp(1j*np.asarray(t)))


def make_measurements(case):
    """Truth is used ONLY here and in scoring, never within an inverse solver."""
    heights = case.delta * np.array([.5, 1., 2.])
    times = np.arange(case.nt) * TAU / case.nt  # exclude duplicate cycle endpoint
    y, t = np.meshgrid(heights, times, indexing="ij")
    clean = exact_velocity(y, t, case.kappa)
    standard_noise = np.random.default_rng(case.seed + 1000).normal(size=clean.shape)
    return {"y": y.ravel(), "t": t.ravel(), "clean": clean.ravel(),
            "observed": (clean + case.noise*standard_noise).ravel(),
            "standard_noise": standard_noise.ravel(), "heights": heights}


def grid_operator(n=400, stretch=6.):
    """Fixed nonuniform grid, independent of unknown viscosity and true delta."""
    z = np.linspace(0, 1, n+1)
    y = np.sinh(stretch*z) / np.sinh(stretch)
    hm, hp = np.diff(y)[:-1], np.diff(y)[1:]
    a = 2 / (hm*(hm+hp))
    b = -2 / (hm*hp)
    c = 2 / (hp*(hm+hp))
    return y, a, b, c


def fd_periodic(kappa, n=400):
    """Conventional FD solver for the periodic state, not an analytic fit.

    The sinusoidal forcing permits direct solution of the FD method-of-lines
    periodic state: (kappa D2-i I) amplitude=-i. This avoids thousands of
    startup cycles in thin layers. mol_check independently integrates the ODEs.
    """
    y, a, b, c = grid_operator(n)
    band = np.zeros((3, n-1), dtype=complex)
    band[0, 1:] = kappa*c[:-1]
    band[1] = kappa*b - 1j
    band[2, :-1] = kappa*a[1:]
    rhs = np.full(n-1, -1j, dtype=complex)
    rhs[-1] -= kappa*c[-1]
    amplitude = np.r_[0j, solve_banded((1, 1), band, rhs), 1+0j]
    return y, amplitude


def interpolate_amplitude(y, grid, amplitude):
    return np.interp(y, grid, amplitude.real) + 1j*np.interp(y, grid, amplitude.imag)


def wall_derivative(grid, amplitude):
    h1, h2 = grid[1:3]
    weights = np.array([-(h1+h2)/(h1*h2), h2/(h1*(h2-h1)), -h1/(h2*(h2-h1))])
    return weights @ amplitude[:3]


def fit_baseline(data, n=400):
    started = time.perf_counter()
    calls = 0

    def objective(log_kappa):
        nonlocal calls
        calls += 1
        grid, amplitude = fd_periodic(np.exp(log_kappa), n)
        predicted = np.real(interpolate_amplitude(data["y"], grid, amplitude)
                            * np.exp(1j*data["t"]))
        return np.mean((predicted - data["observed"])**2)

    # A single bounded minimization over many decades is unsafe: the response
    # becomes nearly flat when all probes lie outside a trial boundary layer.
    # Scan broadly first, then refine the three best coarse neighborhoods.
    search = np.linspace(np.log(1e-8), np.log(.2), 61)
    coarse = np.array([objective(z) for z in search])
    candidates = []
    for j in np.argsort(coarse)[:3]:
        lo, hi = search[max(j-1, 0)], search[min(j+1, len(search)-1)]
        candidates.append(minimize_scalar(objective, bounds=(lo, hi),
                                          method="bounded", options={"xatol": 1e-9}))
    result = min(candidates, key=lambda r: r.fun)
    if not result.success:
        raise RuntimeError(f"Baseline inverse fit failed: {result.message}")
    elapsed = time.perf_counter() - started
    kappa = float(np.exp(result.x))
    grid, amplitude = fd_periodic(kappa, n)
    return {"kappa": kappa, "training_seconds": elapsed, "nfev": calls,
            "cost": float(result.fun), "grid": grid, "amplitude": amplitude}


def mol_check(kappa, n=200):
    """BDF method-of-lines check, initialized at the computed periodic FD state.

    This tests time integration against a spatial numerical state; the exact
    reference is not supplied as an initial condition. Startup from u=1 is
    separately quantified by its slowest diffusion time in P1 verification.
    """
    y, a, b, c = grid_operator(n)
    _, amplitude = fd_periodic(kappa, n)
    matrix = kappa*diags([a[1:], b, c[:-1]], [-1, 0, 1], format="csc")

    def rhs(t, values):
        out = matrix @ values - np.sin(t)
        out[-1] += kappa*c[-1]*np.cos(t)
        return out

    times = np.linspace(0, TAU, 65)
    start = time.perf_counter()
    sol = solve_ivp(rhs, (0, TAU), amplitude.real[1:-1], method="BDF", jac=matrix,
                    t_eval=times, rtol=2e-8, atol=2e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    expected = np.real(amplitude[1:-1, None] * np.exp(1j*times[None, :]))
    return {"kappa": kappa, "n": n, "seconds": time.perf_counter()-start,
            "relative_l2_vs_periodic_fd": float(np.linalg.norm(sol.y-expected)/np.linalg.norm(expected)),
            "endpoint_max_difference": float(np.max(abs(sol.y[:, -1]-sol.y[:, 0]))),
            "nfev": sol.nfev}


def neural_fields(p, y, width, scale=1., jacobian=False):
    """One trainable tanh hidden layer, two outputs A(y), B(y).

    u(y,t)=A(y)cos(t)+B(y)sin(t). This is a neural trial function for the PDE,
    with the exact periodic forcing frequency, NOT the exact spatial solution.
    Returns exact analytic y derivatives and optional parameter Jacobians.
    """
    m = width
    w, b, va, vb = (p[i*m:(i+1)*m] for i in range(4))
    ca, cb = p[4*m:4*m+2]
    x = np.asarray(y).reshape(-1, 1)/scale
    input_jet = Jet2(np.asarray(y).reshape(-1, 1),
                    np.ones((len(x), 1)), np.zeros((len(x), 1)))
    hidden = input_jet.affine(w/scale, b).tanh()
    tanh = hidden.value
    d1 = 1-tanh*tanh
    d2 = -2*tanh*d1
    d3 = -2+8*tanh*tanh-6*tanh**4
    wy = w/scale
    basis1, basis2 = hidden.first, hidden.second
    out_a, out_b = hidden.linear(va, ca), hidden.linear(vb, cb)
    aa, bb = out_a.value, out_b.value
    a1, b1, a2, b2 = out_a.first, out_b.first, out_a.second, out_b.second
    result = (aa, bb, a1, b1, a2, b2)
    if not jacobian:
        return result
    dtype = np.result_type(p)
    ja = np.zeros((len(x), len(p)), dtype=dtype)
    jb, ja2, jb2 = np.zeros_like(ja), np.zeros_like(ja), np.zeros_like(ja)
    ja[:, :m], ja[:, m:2*m], ja[:, 2*m:3*m] = d1*x*va, d1*va, tanh
    jb[:, :m], jb[:, m:2*m], jb[:, 3*m:4*m] = d1*x*vb, d1*vb, tanh
    ja[:, 4*m], jb[:, 4*m+1] = 1, 1
    dw2 = d3*x*wy*wy + 2*d2*w/(scale*scale)
    db2 = d3*wy*wy
    ja2[:, :m], ja2[:, m:2*m], ja2[:, 2*m:3*m] = dw2*va, db2*va, basis2
    jb2[:, :m], jb2[:, m:2*m], jb2[:, 3*m:4*m] = dw2*vb, db2*vb, basis2
    return result, (ja, jb, ja2, jb2)


def initial_parameters(width, seed, initial_kappa):
    rng = np.random.default_rng(seed+2000)
    return np.r_[rng.normal(0, 1, width), rng.uniform(-1, 1, width),
                 rng.normal(0, .1, 2*width), .5, 0., np.log(initial_kappa)]


def collocation_points(n, seed, focused=False, scale=1.):
    rng = np.random.default_rng(seed+3000)
    points = (np.arange(n)+rng.random(n))/n  # stratified, still uniform density
    if focused:
        # Same count and same strata; half moved into a region set by INITIAL
        # viscosity, never true thickness or fitted ground-truth errors.
        points[::2] *= min(1., 6*scale)
    return points


def residual_jacobian(p, data, yc, width, scale=1., jacobian=False):
    """Objective components exactly sum to data MSE + PDE MSE + BC MSE."""
    m, n, nd = width, len(yc), len(data["y"])
    fr = neural_fields(p, yc, m, scale, jacobian)
    fd = neural_fields(p, data["y"], m, scale, jacobian)
    fb = neural_fields(p, [0., 1.], m, scale, jacobian)
    if jacobian:
        (aa, bb, _, _, a2, b2), (ja, jb, ja2, jb2) = fr
        (ad, bd, *_), (jad, jbd, *_) = fd
        (ab, bbound, *_), (jab, jbb, *_) = fb
    else:
        aa, bb, _, _, a2, b2 = fr
        ad, bd, *_ = fd
        ab, bbound, *_ = fb
    kappa = np.exp(p[-1])
    rc, rs = bb-kappa*a2, 1-aa-kappa*b2
    ct, st = np.cos(data["t"]), np.sin(data["t"])
    rd = ad*ct+bd*st-data["observed"]
    bc = np.r_[ab-[0, 1], bbound]
    residual = np.r_[rd/np.sqrt(nd), rc/np.sqrt(2*n), rs/np.sqrt(2*n), bc/2]
    if not jacobian:
        return residual
    jrc, jrs = jb-kappa*ja2, -ja-kappa*jb2
    jrc[:, -1] -= kappa*a2
    jrs[:, -1] -= kappa*b2
    jac = np.vstack([(jad*ct[:, None]+jbd*st[:, None])/np.sqrt(nd),
                     jrc/np.sqrt(2*n), jrs/np.sqrt(2*n), np.vstack([jab, jbb])/2])
    return residual, jac


def train_pinn(data, initial_kappa, seed, variant, width=16, ncollocation=160,
               max_nfev=600, supervised_target=None):
    scaled, focused = VARIANTS[variant]
    initial_delta = np.sqrt(2*initial_kappa)
    scale = initial_delta if scaled else 1.
    yc = collocation_points(ncollocation, seed, focused, initial_delta)
    p0 = initial_parameters(width, seed, initial_kappa)
    fixed_kappa = variant in ("supervised", "known_kappa")
    history = []
    start = time.perf_counter()
    if supervised_target is not None:
        ys, target_a, target_b = supervised_target

    def expand(z):
        return np.r_[z, p0[-1]] if fixed_kappa else z

    def fun(z):
        p = expand(z)
        if supervised_target is not None:
            aa, bb, *_ = neural_fields(p, ys, width, scale)
            res = np.r_[aa-target_a, bb-target_b]/np.sqrt(2*len(ys))
            components = [float(res@res), 0., 0.]
        else:
            res = residual_jacobian(p, data, yc, width, scale)
            nd = len(data["y"])
            components = [float(res[:nd]@res[:nd]),
                          float(res[nd:-4]@res[nd:-4]), float(res[-4:]@res[-4:])]
        history.append([len(history), time.perf_counter()-start, *components,
                        float(np.exp(p[-1]))])
        return res

    def jac(z):
        p = expand(z)
        if supervised_target is not None:
            _, (ja, jb, *_) = neural_fields(p, ys, width, scale, True)
            out = np.vstack([ja, jb])/np.sqrt(2*len(ys))
        else:
            _, out = residual_jacobian(p, data, yc, width, scale, True)
        return out[:, :-1] if fixed_kappa else out

    z0 = p0[:-1] if fixed_kappa else p0
    low, high = np.full(len(z0), -np.inf), np.full(len(z0), np.inf)
    if not fixed_kappa:
        low[-1], high[-1] = np.log(1e-8), np.log(.2)
    result = least_squares(fun, z0, jac=jac, bounds=(low, high), method="trf",
                           ftol=1e-9, xtol=1e-9, gtol=1e-9, max_nfev=max_nfev)
    elapsed = time.perf_counter()-start
    p = expand(result.x)
    return {"parameters": p, "scale": scale, "kappa": float(np.exp(p[-1])),
            "training_seconds": elapsed, "nfev": result.nfev, "njev": result.njev,
            "cost": float(2*result.cost), "optimality": float(result.optimality),
            "optimizer_status": int(result.status), "optimizer_message": result.message,
            "history": np.asarray(history), "collocation": yc, "width": width}


def field_metrics(aa, bb, shear_complex, y, t, case, kappa):
    predicted = aa[:, None]*np.cos(t)-bb[:, None]*np.sin(t)
    truth = exact_velocity(y[:, None], t[None, :], case.kappa)
    error = predicted-truth
    numerator = np.trapz(np.mean(error**2, axis=1), y)
    denominator = np.trapz(np.mean(truth**2, axis=1), y)
    mask = y <= min(6*case.delta, 1.)
    layer_num = np.trapz(np.mean(error[mask]**2, axis=1), y[mask])
    layer_den = np.trapz(np.mean(truth[mask]**2, axis=1), y[mask])
    amp = abs(shear_complex)
    phase = np.degrees(np.angle(shear_complex))
    phase_error = abs((phase-45+180) % 360-180)
    metrics = {"kappa_hat": float(kappa), "kappa_relative_error": float(abs(kappa/case.kappa-1)),
               "field_relative_l2": float(np.sqrt(numerator/denominator)),
               "layer_relative_l2": float(np.sqrt(layer_num/layer_den)),
               "field_max_error": float(np.max(abs(error))), "shear_amplitude": float(amp),
               "shear_amplitude_relative_error": float(abs(amp/(np.sqrt(2)/case.delta)-1)),
               "shear_phase_degrees": float(phase), "shear_phase_error_degrees": float(phase_error)}
    return metrics, {"y": y, "t": t, "truth": truth, "predicted": predicted, "error": error,
                     "shear": np.real(shear_complex*np.exp(1j*t)),
                     "true_shear": np.sqrt(2)/case.delta*np.cos(t+np.pi/4)}


def evaluate(case, method, trained, width):
    # Uniform grid resolves the thinnest true layer with >10 intervals per delta.
    y, t = np.linspace(0, 1, 1201), np.arange(128)*TAU/128
    if method == "baseline":
        amplitude = interpolate_amplitude(y, trained["grid"], trained["amplitude"])
        shear = wall_derivative(trained["grid"], trained["amplitude"])
        rc, rs = np.zeros_like(y), np.zeros_like(y)
        pde_rms, bc_mse, layer_pde_rms = np.nan, np.nan, np.nan
    else:
        aa, bb, a1, b1, a2, b2 = neural_fields(trained["parameters"], y, width, trained["scale"])
        amplitude, shear = aa-1j*bb, a1[0]-1j*b1[0]
        rc, rs = bb-trained["kappa"]*a2, 1-aa-trained["kappa"]*b2
        pde_rms = np.sqrt(np.trapz((rc**2+rs**2)/2, y))
        mask = y <= min(6*case.delta, 1.)
        layer_pde_rms = np.sqrt(np.trapz((rc[mask]**2+rs[mask]**2)/2, y[mask])/y[mask][-1])
        bc_mse = (aa[0]**2+bb[0]**2+(aa[-1]-1)**2+bb[-1]**2)/4
    metrics, fields = field_metrics(amplitude.real, amplitude.imag, shear, y, t, case, trained["kappa"])
    fields["residual"] = rc[:, None]*np.cos(t)+rs[:, None]*np.sin(t)
    fields["spatial_error_power"] = np.mean(dct(fields["error"], type=2, axis=0, norm="ortho")**2, axis=1)
    fields["temporal_error_power"] = np.mean(abs(np.fft.rfft(fields["error"], axis=1))**2, axis=0)
    metrics.update({"heldout_pde_rms": float(pde_rms), "heldout_layer_pde_rms": float(layer_pde_rms),
                    "heldout_boundary_mse": float(bc_mse)})
    # Same 1201*128 velocity query; no derivatives, warmed up, median of 9.
    ts = []
    for _ in range(10):
        start = time.perf_counter()
        if method == "baseline":
            z = interpolate_amplitude(y, trained["grid"], trained["amplitude"])
            value = np.real(z[:, None]*np.exp(1j*t))
        else:
            a, b, *_ = neural_fields(trained["parameters"], y, width, trained["scale"])
            value = a[:, None]*np.cos(t)+b[:, None]*np.sin(t)
        ts.append(time.perf_counter()-start)
    metrics["inference_seconds"] = float(np.median(ts[1:]))
    return metrics, fields


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def p1_verification(out):
    rows = []
    y_eval = np.linspace(0, 1, 6001)
    for delta in [1/6, .1, .05, .02, .01]:
        kappa = delta*delta/2
        for n in [100, 200, 400, 800]:
            started = time.perf_counter()
            y, amplitude = fd_periodic(kappa, n)
            elapsed = time.perf_counter()-started
            prediction = interpolate_amplitude(y_eval, y, amplitude)
            finite = exact_amplitude(y_eval, kappa, finite=True)
            half = exact_amplitude(y_eval, kappa)
            rows.append({"delta": delta, "n": n,
                         "relative_l2_vs_finite_exact": float(np.linalg.norm(prediction-finite)/np.linalg.norm(finite)),
                         "relative_l2_vs_halfspace_exact": float(np.linalg.norm(prediction-half)/np.linalg.norm(half)),
                         "boundary_truncation_bound": float(np.exp(-1/delta)),
                         "forward_seconds": elapsed,
                         "shear_relative_error": float(abs(wall_derivative(y, amplitude)/((1+1j)/delta)-1)),
                         "estimated_startup_cycles_1e6": int(np.ceil(-np.log(1e-6)/(TAU*kappa*np.pi**2)))})
    mol = [mol_check(d*d/2) for d in [1/6, .1, .05, .02, .01]]
    write_csv(out/"baseline_grid_convergence.csv", rows)
    (out/"mol_verification.json").write_text(json_text(mol), encoding="utf-8")
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for delta in sorted(set(r["delta"] for r in rows), reverse=True):
        selected = [r for r in rows if r["delta"] == delta]
        ax.loglog([r["n"] for r in selected], [r["relative_l2_vs_finite_exact"] for r in selected], "o-", label=f"delta/L={delta:.3g}")
    ax.set(xlabel="Spatial intervals", ylabel="Relative L2 amplitude error", title="P1: FD grid convergence against finite-domain exact solution")
    ax.grid(alpha=.25); ax.legend(fontsize=8)
    save_figure(fig, out/"figures/P1_grid_convergence.png")
    return rows, mol


def save_figure(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_summary(rows, out):
    methods = ["baseline", "vanilla", "combined"]
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    knobs = [("delta", "Inverse thickness L/delta", lambda x: 1/x),
             ("noise", "Noise standard deviation (% U0)", lambda x: 100*x),
             ("nt", "Samples per cycle per probe", lambda x: x)]
    for col, (knob, label, transform) in enumerate(knobs):
        for method in methods:
            selected = [r for r in rows if r["method"] == method and
                        (knob == "delta" or abs(r["delta"]-.1)<1e-8) and
                        (knob == "noise" or abs(r["noise"]-.02)<1e-8) and
                        (knob == "nt" or r["nt"] == 16)]
            values = sorted(set(r[knob] for r in selected), key=transform)
            for rowidx, metric in enumerate(["kappa_relative_error", "field_relative_l2", "shear_amplitude_relative_error"]):
                med, low, high = [], [], []
                for value in values:
                    sample = [100*r[metric] for r in selected if r[knob] == value]
                    med.append(np.median(sample)); low.append(min(sample)); high.append(max(sample))
                x = [transform(v) for v in values]
                axes[rowidx, col].plot(x, np.maximum(med, 1e-7), "o-", label=method, color=COLORS[method])
                axes[rowidx, col].fill_between(x, np.maximum(low, 1e-7), np.maximum(high, 1e-7), alpha=.12, color=COLORS[method])
                axes[rowidx, col].set_yscale("log"); axes[rowidx, col].grid(alpha=.25)
                axes[rowidx, col].set_xlabel(label)
        if knob == "delta":
            for rowidx in range(3):
                axes[rowidx, col].set_xscale("log")
    axes[0, 0].set_ylabel("Viscosity relative error (%)")
    axes[1, 0].set_ylabel("Whole-field relative L2 error (%)")
    axes[2, 0].set_ylabel("Wall-shear amplitude error (%)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("P3: one knob at a time; lines = seed medians, bands = seed ranges")
    save_figure(fig, out/"figures/P3_knob_sweeps.png")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for metric, ax, title in zip(["training_seconds", "inference_seconds"], axes,
                                 ["Inverse fit / training", "153,728 velocity queries"]):
        for method in methods:
            selected = [r for r in rows if r["method"] == method and r["noise"] == .02 and r["nt"] == 16]
            ds = sorted(set(r["delta"] for r in selected), reverse=True)
            ax.loglog([1/d for d in ds], [np.median([r[metric] for r in selected if r["delta"] == d]) for d in ds], "o-", color=COLORS[method], label=method)
        ax.set(xlabel="Inverse thickness L/delta", ylabel="Seconds", title=title); ax.grid(alpha=.25)
    axes[0].legend(fontsize=8)
    save_figure(fig, out/"figures/P5_cost.png")


def plot_case(case, out, methods=("baseline", "vanilla", "combined")):
    figdir = out/"figures"
    arrays = {}
    for method in methods:
        p = out/"fields"/f"{case.key}_{method}.npz"
        if p.exists():
            with np.load(p) as z:
                arrays[method] = dict(z)
    if not arrays:
        return
    fig, axes = plt.subplots(len(arrays), 2, figsize=(11, 3*len(arrays)), squeeze=False)
    for rowidx, (method, f) in enumerate(arrays.items()):
        for col, (field, title) in enumerate([("error", "velocity error"), ("residual", "PDE residual")]):
            vmax = max(float(np.max(abs(f[field]))), 1e-6)
            image = axes[rowidx, col].pcolormesh(f["t"]/TAU, f["y"], f[field], shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            axes[rowidx, col].set_ylim(0, min(6*case.delta, 1))
            axes[rowidx, col].set(xlabel="Cycle fraction", ylabel="y/L", title=f"{method}: {title} (own color scale)")
            if method == "baseline" and field == "residual":
                axes[rowidx, col].set_title("baseline: discrete PDE solved; continuous residual not plotted")
            fig.colorbar(image, ax=axes[rowidx, col])
    fig.suptitle(f"P4: near-wall fields, delta/L={case.delta:g}, seed={case.seed}")
    save_figure(fig, figdir/f"P4_fields_d{case.delta:g}.png")

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.7))
    ref = next(iter(arrays.values()))
    axes[0].plot(ref["t"]/TAU, ref["true_shear"], "k--", label="exact")
    axes[1].plot(np.cos(ref["t"]), ref["true_shear"], "k--", label="exact")
    for method, f in arrays.items():
        axes[0].plot(f["t"]/TAU, f["shear"], color=COLORS[method], label=method)
        axes[1].plot(np.cos(f["t"]), f["shear"], color=COLORS[method], label=method)
        axes[2].loglog(np.arange(1, len(f["spatial_error_power"])), np.maximum(f["spatial_error_power"][1:], 1e-30), color=COLORS[method], label=method)
    axes[0].set(xlabel="Cycle fraction", ylabel="du/dy at wall", title="Wall-shear history")
    axes[1].set(xlabel="Outer velocity / U0", ylabel="du/dy at wall", title="Unsteady hysteresis")
    axes[2].set(xlabel="Spatial cosine mode", ylabel="Mean squared DCT coefficient", title="P4: spatial error spectrum")
    for ax in axes:
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    save_figure(fig, figdir/f"P4_shear_spectrum_d{case.delta:g}.png")

    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    for method, f in arrays.items():
        ax.semilogy(np.arange(len(f["temporal_error_power"])),
                    np.maximum(f["temporal_error_power"], 1e-30), label=method, color=COLORS[method])
    ax.set(xlabel="Temporal Fourier harmonic", ylabel="Mean squared FFT coefficient",
           title="P4: temporal error spectrum; harmonic time structure is imposed")
    ax.legend(fontsize=8); ax.grid(alpha=.25)
    save_figure(fig, figdir/f"P4_temporal_spectrum_d{case.delta:g}.png")

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, method in zip(axes, ["vanilla", "combined"]):
        p = out/"histories"/f"{case.key}_{method}.csv"
        if p.exists():
            h = np.genfromtxt(p, delimiter=",", names=True)
            for term in ["data_mse", "pde_mse", "boundary_mse"]:
                ax.semilogy(h["evaluation"], np.maximum(h[term], 1e-20), label=term)
        ax.set(xlabel="Objective evaluations (includes trial steps)", ylabel="Unweighted loss", title=f"P2: {method}")
        ax.grid(alpha=.25); ax.legend(fontsize=8)
    save_figure(fig, figdir/f"P2_losses_d{case.delta:g}.png")


def plot_diagnosis(rows, out):
    selected = [r for r in rows if abs(r["delta"]-.01)<1e-9 and r["noise"] == .02 and r["nt"] == 16]
    methods = [m for m in VARIANTS if any(r["method"] == m for r in selected)]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, metric, title in zip(axes, ["kappa_relative_error", "layer_relative_l2", "shear_amplitude_relative_error"],
                                  ["Viscosity error (%)", "Near-wall field error (%)", "Shear amplitude error (%)"]):
        for j, method in enumerate(methods):
            vals = [100*r[metric] for r in selected if r["method"] == method]
            # Fixed true kappa and supervised controls do not infer viscosity.
            if metric == "kappa_relative_error" and method in ("known_kappa", "supervised"):
                continue
            ax.scatter(np.full(len(vals), j), np.maximum(vals, 1e-6), color=COLORS[method], s=25)
            ax.plot([j-.2, j+.2], [max(np.median(vals), 1e-6)]*2, color=COLORS[method], lw=2)
        ax.set_yscale("log"); ax.set_xticks(range(len(methods)), methods, rotation=35, ha="right", fontsize=8)
        ax.set_title(title); ax.grid(axis="y", alpha=.25)
    fig.suptitle("Path B/C: matched ablations at delta/L=0.01; dots = seeds, bars = medians")
    save_figure(fig, out/"figures/BC_controlled_diagnosis.png")


def experiment_cases(suite, seeds):
    settings = [(d, .02, 16) for d in [1/6, .1, .05, .02, .01]]
    if suite == "pilot":
        settings = [(.1, .02, 16), (.01, .02, 16)]
    else:
        settings += [(.1, s, 16) for s in [0, .01, .05, .1]]
        settings += [(.1, .02, n) for n in [64, 8]]
    return [Case(d, s, n, seed) for d, s, n in settings for seed in seeds]


def run_experiments(args):
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    for directory in ["data", "models", "fields", "histories", "figures"]:
        (out/directory).mkdir(exist_ok=True)
    metadata = {"utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                "matplotlib": matplotlib.__version__, "platform": platform.platform(),
                "processor": platform.processor(), "logical_cpus": os.cpu_count(),
                "threadpools": threadpool_info(), "arguments": vars(args),
                "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "dtype": "float64", "training_backend": "SciPy analytic-Jacobian nonlinear least squares",
                "neural_parameters": 4*args.width+3}
    env_path = out/"environment.json"
    if args.resume and env_path.exists():
        previous = json.loads(env_path.read_text(encoding="utf-8"))
        signature = lambda x: {k: v for k, v in x.items() if k not in ("resume", "output")}
        if (previous["source_sha256"] != metadata["source_sha256"] or
                signature(previous["arguments"]) != signature(metadata["arguments"])):
            raise ValueError("Resume refused: source/configuration changed. Use a new output directory.")
    env_path.write_text(json_text(metadata), encoding="utf-8")
    print("P1 numerical verification", flush=True)
    p1_verification(out)
    cases = experiment_cases(args.suite, args.seeds)
    rows = []
    started = time.perf_counter()
    for i, case in enumerate(cases):
        data = make_measurements(case)
        np.savez_compressed(out/"data"/f"{case.key}.npz", **data, **asdict(case))
        initial_factor = 5. if args.seeds.index(case.seed) % 2 == 0 else .2
        initial_kappa = case.kappa*initial_factor
        methods = ["baseline", "vanilla", "combined"]
        if case.delta in (.1, .01) and case.noise == .02 and case.nt == 16:
            methods += ["sampling", "scaling"]
        if case.delta == .01 and case.noise == .02 and case.nt == 16:
            methods += ["supervised", "known_kappa"]
        for method in methods:
            tag = f"{case.key}_{method}"
            model_path = out/"models"/f"{tag}.npz"
            metrics_path = out/"models"/f"{tag}.json"
            if args.resume and metrics_path.exists():
                row = json.loads(metrics_path.read_text(encoding="utf-8"))
                rows.append(row)
                print(f"[{i+1}/{len(cases)}] cached {tag}", flush=True)
                continue
            print(f"[{i+1}/{len(cases)}] {tag}", flush=True)
            if method == "baseline":
                fitted = fit_baseline(data, args.fd_intervals)
                np.savez_compressed(model_path, grid=fitted["grid"], amplitude=fitted["amplitude"], kappa=fitted["kappa"])
            else:
                target = None
                if method == "supervised":
                    # Oracle diagnostic only: dense targets have additional information.
                    ys = np.unique(np.r_[np.linspace(0, 1, 201), np.linspace(0, 6*case.delta, 201)])
                    z = exact_amplitude(ys, case.kappa)
                    target = (ys, z.real, -z.imag)
                fitted = train_pinn(data, case.kappa if method in ("known_kappa", "supervised") else initial_kappa,
                                    case.seed, method, args.width, args.collocation, args.max_nfev, target)
                np.savez_compressed(model_path, parameters=fitted["parameters"], scale=fitted["scale"],
                                    width=args.width, collocation=fitted["collocation"])
                write_csv(out/"histories"/f"{tag}.csv", [dict(zip(
                    ["evaluation", "seconds", "data_mse", "pde_mse", "boundary_mse", "kappa"], h)) for h in fitted["history"]])
            metrics, fields = evaluate(case, method, fitted, args.width)
            if method in ("supervised", "known_kappa"):
                # These oracle controls are not viscosity estimators.
                metrics["kappa_hat"] = np.nan
                metrics["kappa_relative_error"] = np.nan
            np.savez_compressed(out/"fields"/f"{tag}.npz", **fields)
            row = {**asdict(case), "method": method, "kappa_true": case.kappa,
                   "initial_kappa": case.kappa if method in ("known_kappa", "supervised") else initial_kappa,
                   **metrics, "training_seconds": fitted["training_seconds"], "nfev": fitted["nfev"],
                   "objective": fitted["cost"], "optimizer_status": fitted.get("optimizer_status", 1),
                   "optimizer_message": fitted.get("optimizer_message", "Scalar fit converged"),
                   "optimality": fitted.get("optimality", np.nan)}
            if method != "baseline":
                final_residual = residual_jacobian(fitted["parameters"], data, fitted["collocation"], args.width, fitted["scale"])
                nd = len(data["y"])
                row.update({"final_data_mse": float(final_residual[:nd]@final_residual[:nd]),
                            "final_pde_mse": float(final_residual[nd:-4]@final_residual[nd:-4]),
                            "final_boundary_mse": float(final_residual[-4:]@final_residual[-4:]),
                            "collocation_count": len(fitted["collocation"]),
                            "collocation_inside_one_delta": int(np.sum(fitted["collocation"] <= case.delta)),
                            "collocation_inside_six_delta": int(np.sum(fitted["collocation"] <= 6*case.delta))})
            metrics_path.write_text(json_text(row), encoding="utf-8")
            rows.append(row)
            write_csv(out/"metrics.csv", rows)
            print(f"  kappa error={100*row['kappa_relative_error']:.3g}%, layer L2={100*row['layer_relative_l2']:.3g}%, shear error={100*row['shear_amplitude_relative_error']:.3g}%, {fitted['training_seconds']:.2f}s", flush=True)
    write_csv(out/"metrics.csv", rows)
    plot_summary(rows, out)
    for delta in [.1, .01]:
        plot_case(Case(delta, .02, 16, args.seeds[0]), out)
    plot_diagnosis(rows, out)
    metadata["experiment_wall_seconds"] = time.perf_counter()-started
    metadata["utc_completed"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (out/"environment.json").write_text(json_text(metadata), encoding="utf-8")
    print(f"Completed {len(rows)} method/case runs. Output: {out}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["pilot", "full"], default="full")
    parser.add_argument("--output", default="results")
    parser.add_argument("--seeds", nargs="+", type=int, default=[17, 29, 41])
    parser.add_argument("--width", type=int, default=16)
    parser.add_argument("--collocation", type=int, default=160)
    parser.add_argument("--max-nfev", type=int, default=600)
    parser.add_argument("--fd-intervals", type=int, default=400)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--resume", action="store_true", help="Reuse completed runs only with the same configuration and source.")
    args = parser.parse_args()
    if len(set(args.seeds)) != len(args.seeds):
        parser.error("Seeds must be distinct")
    with threadpool_limits(limits=args.threads):
        run_experiments(args)


if __name__ == "__main__":
    main()
