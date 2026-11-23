"""Nonrotating uniform cantilever with point tip mass: NACA TN 3459, p. 17.

Yntema (1955), equations (13)-(14):
https://ntrs.nasa.gov/api/citations/19930084187/downloads/19930084187.pdf
SI units. The beam is stationary; motor RPM is not a beam rotation rate.
No tip rotary inertia, axial preload, damping or gyroscopic terms are included.
"""
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq


METHOD = "naca_tn_3459_eq13_14"
SOURCE = "https://ntrs.nasa.gov/api/citations/19930084187/downloads/19930084187.pdf"


def characteristic(beta, ratio):
    """Equation (13) divided by cosh(beta), avoiding hyperbolic overflow."""
    e = np.exp(-beta)
    sech = 2*e/(1+e*e)
    return sech + np.cos(beta) + ratio*beta*(np.cos(beta)*np.tanh(beta)-np.sin(beta))


@lru_cache(maxsize=128)
def tip_mass_roots(ratio, count=4):
    """One positive root per interval ((n-1)*pi,n*pi); no zero rigid mode."""
    if not np.isfinite(ratio) or ratio < 0:
        raise ValueError("Tip/beam mass ratio must be finite and nonnegative")
    if not isinstance(count, (int, np.integer)) or count < 1:
        raise ValueError("Mode count must be a positive integer")
    return tuple(brentq(characteristic, n*np.pi, (n+1)*np.pi,
                        args=(ratio,), xtol=1e-14, rtol=1e-14) for n in range(count))


def _parameters(length, EI, mass_per_length, tip_mass):
    if any(not np.isfinite(v) or v <= 0 for v in (length, EI, mass_per_length)):
        raise ValueError("Length, EI and line mass must be positive and finite")
    if not np.isfinite(tip_mass) or tip_mass < 0:
        raise ValueError("Tip mass must be finite and nonnegative")
    return tip_mass/(mass_per_length*length)


def tip_mass_frequencies(length, EI, mass_per_length, tip_mass, count=4):
    """Return f_n in Hz from the exact positive roots of equation (13)."""
    ratio = _parameters(length, EI, mass_per_length, tip_mass)
    beta = np.array(tip_mass_roots(ratio, count))
    return beta**2/(2*np.pi*length**2)*np.sqrt(EI/mass_per_length)


def first_frequency_partials(length, EI, mass_per_length, tip_mass):
    """Implicit-root sensitivities wrt (length, EI, line mass, tip mass).

    Differentiate H(beta,r)=0, r=tip_mass/(line_mass*length). These are
    continuum frequency derivatives, not FE eigenvector sensitivities.
    """
    ratio = _parameters(length, EI, mass_per_length, tip_mass)
    beta = tip_mass_roots(ratio, 1)[0]
    f = beta**2/(2*np.pi*length**2)*np.sqrt(EI/mass_per_length)
    e = np.exp(-beta)
    sech, tanh = 2*e/(1+e*e), np.tanh(beta)
    c, s = np.cos(beta), np.sin(beta)
    a = c*tanh-s
    hb = -sech*tanh-s+ratio*(a-beta*tanh*(s+c*tanh))
    beta_r = -beta*a/hb
    sensitivity = 2*ratio*beta_r/beta
    return f, dict(length=f/length*(-2-sensitivity), EI=f/(2*EI),
                   mass_per_length=f/mass_per_length*(-.5-sensitivity),
                   tip_mass=f*2*beta_r/(beta*mass_per_length*length))


def _shape(beta, fraction):
    """Equation (14) and its derivative wrt z=beta*x/L, stably evaluated."""
    z = beta*np.asarray(fraction)
    e = np.exp(-beta)
    c, s = np.cos(beta), np.sin(beta)
    den = 1+e*e+2*e*c
    a = (1-e*e+2*e*s)/den
    growing = np.exp(z-beta)*(e+c-s)/den
    decaying = np.exp(-z)*(1+e*(c+s))/den
    y = growing-decaying-np.sin(z)+a*np.cos(z)
    dz = growing+decaying-np.cos(z)-a*np.sin(z)
    return y, dz


@lru_cache(maxsize=8)
def _quadrature(count):
    x, w = np.polynomial.legendre.leggauss(max(64, 4*count))
    return (x+1)/2, w/2


def tip_mass_modes(length, EI, mass_per_length, tip_mass, positions, count=4):
    """Exact frequencies and mass-normalized (displacement, slope) samples.

    Normalization is integral(mu*y^2 dx)+tip_mass*y(L)^2=1, including the tip
    mass once. Returned shape columns follow [w0,slope0,w1,slope1,...].
    """
    ratio = _parameters(length, EI, mass_per_length, tip_mass)
    x = np.asarray(positions, float)
    if x.ndim != 1 or not np.isfinite(x).all() or np.any(x < 0) or np.any(x > length):
        raise ValueError("Mode sample positions must lie in [0,length]")
    roots = tip_mass_roots(ratio, count)
    gx, gw = _quadrature(count)
    shapes = np.empty((2*len(x), count))
    for j, beta in enumerate(roots):
        y, dz = _shape(beta, x/length)
        gy, _ = _shape(beta, gx)
        tip, _ = _shape(beta, 1.)
        norm = np.sqrt(mass_per_length*length*(gw@gy**2)+tip_mass*tip**2)
        sign = 1 if tip >= 0 else -1
        shapes[::2, j] = sign*y/norm
        shapes[1::2, j] = sign*dz*beta/(length*norm)
    shapes[np.repeat(x == 0, 2)] = 0.
    return tip_mass_frequencies(length, EI, mass_per_length, tip_mass, count), shapes


def modal_summary(length, EI, mass_per_length, tip_mass, count=4):
    ratio = _parameters(length, EI, mass_per_length, tip_mass)
    roots = tip_mass_roots(ratio, count)
    return dict(method=METHOD, source=SOURCE, printed_page=17, equations=[13, 14],
                beam_rotation_rad_s=0., tip_to_beam_mass_ratio=float(ratio),
                beta=list(roots), beta_squared=[b*b for b in roots],
                frequencies_Hz=tip_mass_frequencies(length, EI, mass_per_length, tip_mass, count).tolist(),
                assumptions="Uniform EI and line mass; fixed root; point tip mass without rotary inertia; unloaded stationary beam")
