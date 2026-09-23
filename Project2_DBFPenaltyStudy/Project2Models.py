"""Smooth DBF aircraft model for Project 2 ill-conditioning study.

This module preserves the main sizing, drag, mission timing, propulsion, and
scoring relationships from Project 1, while replacing discrete operations
(floor, throttle-grid lookup, stepwise turn load-factor reduction) with smooth
or continuous counterparts so finite-difference gradients/Hessians are
meaningful.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq


# Original Project 1 bounds, ordered as:
# [wing_area, fuselage_length, SW2, SW3_reduction, V2, V3]
PROJECT1_BOUNDS = np.array(
    [
        [2.0, 10.0],
        [5.0, 7.0],
        [4.0, 12.0],
        [0.0, 6.0],
        [40.0, 140.0],
        [40.0, 140.0],
    ],
    dtype=float,
)

# Project 2 conditioning-study bounds.  The SW2 upper search bound is widened
# from 12 to 16 lb so the 20-lb gross-weight constraint, rather than the
# original payload box bound, controls the high-payload optimum.  This is a
# deliberate study modification and should be disclosed in the report.
BOUNDS = PROJECT1_BOUNDS.copy()
BOUNDS[2, 1] = 16.0

X0 = np.array([6.0, 6.0, 5.0, 5.0, 80.0, 80.0], dtype=float)


def softmin(a: float, b: float, tau: float = 0.10) -> float:
    """Differentiable approximation to min(a, b), stable for positive a,b."""
    m = min(a, b)
    return m - tau * math.log(math.exp(-(a - m) / tau) + math.exp(-(b - m) / tau))


class SmoothPlaneModel:
    def __init__(self, wing_area: float, fuselage_length: float, sw2: float):
        self.wing_area = float(wing_area)
        self.fuselage_length = float(fuselage_length)
        self.sw2 = float(sw2)
        self.sensor_size()
        self.geo_sizing()

    def sensor_size(self) -> None:
        aspect_ratio = 6.0
        leadshot_density = 430.0  # lb/ft^3
        sensor_volume = self.sw2 / leadshot_density
        sensor_length = (sensor_volume * (aspect_ratio * 2.0) ** 2 / math.pi) ** (1.0 / 3.0)

        cd0_nose = 0.2
        cf_body = 0.0042
        self.sensor_diameter = sensor_length / aspect_ratio
        self.sensor_drag_area = (
            cd0_nose * self.sensor_diameter**2 * math.pi / 4.0
            + cf_body * math.pi * self.sensor_diameter * sensor_length
        )

    def geo_sizing(self) -> None:
        fuselage_diameter = self.sensor_diameter * 1.5
        self.aspect_ratio = 6.0**2 / self.wing_area
        self.e = 0.75
        self.k_wing = 1.0 / (math.pi * self.e * self.aspect_ratio)
        self.chord = 6.0 / self.aspect_ratio

        fuselage_weight = (
            3.5 * math.pi * (fuselage_diameter / 2.0) ** 2 * self.fuselage_length
        )
        electronics_weight = 2.5
        wing_weight = (2.0 / 7.0) * self.wing_area
        empennage_weight = 1.0

        self.empty_weight = (
            fuselage_weight + electronics_weight + wing_weight + empennage_weight
        )

        # Preserve the original Project 1 CG expression.
        self.cgx = (
            electronics_weight * 0.5
            + empennage_weight * self.fuselage_length
            + fuselage_weight * self.fuselage_length / 2.0
            + wing_weight * self.chord / 10.0
        ) / (self.empty_weight - wing_weight)

        lt = self.fuselage_length - (self.cgx + self.chord / 10.0)
        self.sv = 0.07 * self.wing_area * 6.0 / lt
        self.sh = 0.7 * self.wing_area * self.chord / lt

        cd0_fuselage = 0.3
        cf_skin = 0.0042
        self.drag_area = (
            cd0_fuselage * fuselage_diameter**2 * math.pi / 4.0
            + cf_skin
            * (
                math.pi * fuselage_diameter * self.fuselage_length
                + 2.0 * (self.sh + self.sv)
            )
        )

    def lift_coefficient(self, velocity: float, lift: float) -> float:
        rho = 0.0023769
        q = 0.5 * rho * velocity**2
        return lift / (q * self.wing_area)

    def drag(self, velocity: float, lift: float, sensor_deployed: bool) -> float:
        rho = 0.0023769
        q = 0.5 * rho * velocity**2
        cl = self.lift_coefficient(velocity, lift)
        cd_induced = self.k_wing * cl**2
        cd_wing = 0.02 + cd_induced
        wing_drag = cd_wing * q * self.wing_area
        body_drag = self.drag_area * q
        sensor_drag = self.sensor_drag_area * q if sensor_deployed else 0.0
        return wing_drag + body_drag + sensor_drag

    @staticmethod
    def _prop_speed(throttle: float, drag_over_v2: float) -> float:
        """Original Project 1 prop-speed relation evaluated continuously."""
        thrust_static = 21.175
        vmax = 170.0
        a = 1.4 * throttle - 0.4
        term = thrust_static * a / (drag_over_v2 * vmax * throttle)
        return (-term + math.sqrt(term * term + 4.0 * thrust_static * a / drag_over_v2)) / 2.0

    def calculate_current(self, velocity: float, lift: float, sensor_deployed: bool) -> float:
        total_drag = self.drag(velocity, lift, sensor_deployed)
        drag_over_v2 = total_drag / velocity**2

        def residual(throttle: float) -> float:
            return self._prop_speed(throttle, drag_over_v2) - velocity

        lo, hi = 0.300001, 0.999999
        r_lo, r_hi = residual(lo), residual(hi)

        # A design demanding more than full throttle is deliberately made costly,
        # but the function remains continuous enough for the local analysis region.
        if r_hi < 0.0:
            throttle = hi + min(0.25, -r_hi / 200.0)
        elif r_lo > 0.0:
            throttle = lo
        else:
            throttle = brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=100)

        return 85.0 * (1.8541 * throttle**2 - 1.0736 * throttle + 0.2136)


class SmoothCourse:
    def __init__(
        self,
        wing_area: float,
        fuselage_length: float,
        sw2: float,
        sw3_reduction: float,
        v2: float,
        v3: float,
    ):
        self.wing_area = float(wing_area)
        self.fuselage_length = float(fuselage_length)
        self.sw2 = float(sw2)
        self.sw3 = float(sw2 - sw3_reduction)
        self.v2 = float(v2)
        self.v3 = float(v3)
        self.aircraft = SmoothPlaneModel(wing_area, fuselage_length, sw2)

    @property
    def m2_takeoff_weight(self) -> float:
        # Preserve the Project 1 10% shipping-container weight allowance.
        return self.aircraft.empty_weight + 1.1 * self.sw2

    def _wing_lift(self, gross_weight: float) -> float:
        a = self.aircraft.cgx + self.aircraft.chord / 10.0
        return gross_weight - a * gross_weight / (a - self.fuselage_length)

    def _turn_load_factor(self, velocity: float, level_lift: float) -> float:
        # Original target was 2.5 g. Replace the 0.1-g decrement loop with a
        # continuous CLmax-limited load factor.
        clmax = 1.2
        rho = 0.0023769
        q = 0.5 * rho * velocity**2
        n_clmax = clmax * q * self.wing_area / level_lift
        # Smooth min with a small tau; keep safely above 1 for a finite turn radius.
        n = softmin(2.5, n_clmax, tau=0.02)
        return max(1.05, n)

    def mission_state(self) -> dict[str, float]:
        m2_weight = self.m2_takeoff_weight
        m3_weight = self.aircraft.empty_weight + self.sw3

        lift_m2 = self._wing_lift(m2_weight)
        lift_m3 = self._wing_lift(m3_weight)

        n2 = self._turn_load_factor(self.v2, lift_m2)
        n3 = self._turn_load_factor(self.v3, lift_m3)

        amps2 = self.aircraft.calculate_current(self.v2, lift_m2, False)
        amps3 = self.aircraft.calculate_current(self.v3, lift_m3, True)
        amps2_turn = self.aircraft.calculate_current(self.v2, n2 * lift_m2, False)
        amps3_turn = self.aircraft.calculate_current(self.v3, n3 * lift_m3, True)

        straight_time2 = 2000.0 / self.v2
        straight_time3 = 2000.0 / self.v3
        r2 = self.v2**2 / (32.2 * (n2**2 - 1.0))
        r3 = self.v3**2 / (32.2 * (n3**2 - 1.0))
        turn_time2 = 4.0 * math.pi * r2 / self.v2
        turn_time3 = 4.0 * math.pi * r3 / self.v3

        lap_time2 = straight_time2 + turn_time2
        lap_time3 = straight_time3 + turn_time3
        cap_draw2 = amps2 * straight_time2 + amps2_turn * turn_time2
        cap_draw3 = amps3 * straight_time3 + amps3_turn * turn_time3

        battery_capacity = 3300.0 * 0.75 * 3600.0 / 1000.0
        laps_time = 300.0 / lap_time3
        laps_energy = battery_capacity / cap_draw3
        laps3_continuous = softmin(laps_time, laps_energy, tau=0.05)

        m2_time = 5.0 * lap_time2
        m2_score = self.sw2 / m2_time
        m3_score = self.sw3 * laps3_continuous
        gm_score = self.sw2

        return {
            "empty_weight": self.aircraft.empty_weight,
            "m2_takeoff_weight": m2_weight,
            "m2_time": m2_time,
            "lap_time_m3": lap_time3,
            "laps3": laps3_continuous,
            "m2_score": m2_score,
            "m3_score": m3_score,
            "gm_score": gm_score,
            "n2": n2,
            "n3": n3,
        }


def dbf_objective(x: np.ndarray) -> float:
    """Negative normalized DBF competition score (smooth Project 2 version)."""
    c = SmoothCourse(*np.asarray(x, dtype=float))
    s = c.mission_state()
    return -(s["m2_score"] / 0.15 + s["m3_score"] / 50.0 + s["gm_score"] / 12.0)


def weight_constraint(x: np.ndarray, limit: float = 20.0) -> float:
    """g(x) <= 0 form of the Project 1 team target for M2 gross weight."""
    c = SmoothCourse(*np.asarray(x, dtype=float))
    return c.m2_takeoff_weight - limit


def penalized_objective(x: np.ndarray, rho: float, limit: float = 20.0) -> float:
    """Quadratic exterior penalty for the M2 gross-weight constraint."""
    base = dbf_objective(x)
    g = weight_constraint(x, limit=limit)
    violation = max(0.0, g)
    return base + 0.5 * float(rho) * violation**2


def augmented_lagrangian_objective(
    x: np.ndarray,
    lam: float,
    rho: float,
    limit: float = 20.0,
) -> float:
    """Powell-Hestenes-Rockafellar AL for inequality g(x) <= 0.

    This form is differentiable away from the switching point and avoids driving
    rho to extremely large values solely to enforce feasibility.
    """
    g = weight_constraint(x, limit=limit)
    shifted = max(0.0, g + lam / rho)
    return dbf_objective(x) + 0.5 * rho * shifted**2 - 0.5 * lam**2 / rho


def state(x: np.ndarray) -> dict[str, float]:
    return SmoothCourse(*np.asarray(x, dtype=float)).mission_state()
