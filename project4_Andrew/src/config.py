from dataclasses import dataclass

@dataclass(frozen=True)
class Config:
    # Geometry / discretization
    arm_length: float = 0.6096        # m = 2 ft
    n_elem: int = 20
    wall_thickness: float = 0.0015    # m
    radius_initial: float = 0.035     # m, feasible uniform starting design
    radius_min: float = 0.008         # m
    radius_max: float = 0.040         # m
    shape_bound: float = 0.020        # m, +/- z offset of four control points

    # 6061-T6 aluminum baseline
    E: float = 68.9e9                 # Pa
    nu: float = 0.33
    rho: float = 2700.0               # kg/m^3
    yield_strength: float = 276e6     # Pa
    yield_fs: float = 2.0

    # Industrial coaxial propulsion-station model
    # Two counter-rotating motors share one arm. Each motor produces 15 lbf,
    # so the combined maximum thrust is 30 lbf = 133.45 N.
    motor_mass: float = 0.680388555    # kg = 1.5 lb combined motors/props/mount
    rotor_thrust: float = 133.4466     # N = 30 lbf total at one coaxial station
    rotor_rpm: float = 3500.0          # rpm, nominal maximum operating speed used for separation basis
    maneuver_accel_g: float = 2.0      # g, secondary 2g horizontal-maneuver case aligned with the arm axis
    motor_torque_z: float = 0.0        # N*m; coaxial counter-rotation cancels net reaction torque
    gravity: float = 9.80665

    # Constraints
    tip_deflection_limit: float = 0.006096  # m = 1% of the 2-ft span
    frequency_limit: float = 65.0           # Hz, about 11.4% above 3500/60 = 58.33 Hz shaft frequency
    buckling_factor_limit: float = 3.0
    ks_rho: float = 50.0

    # Optimization
    maxiter: int = 160
    ftol: float = 1e-8

    @property
    def G(self):
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def stress_limit(self):
        return self.yield_strength / self.yield_fs

    @property
    def maneuver_axial_force(self):
        return self.motor_mass * self.maneuver_accel_g * self.gravity

    @property
    def rotor_shaft_hz(self):
        return self.rotor_rpm / 60.0

CFG = Config()
