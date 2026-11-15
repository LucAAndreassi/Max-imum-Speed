from dataclasses import dataclass

@dataclass(frozen=True)
class Config:
    # Geometry / discretization
    arm_length: float = 0.30          # m
    n_elem: int = 20
    wall_thickness: float = 0.0015    # m
    radius_initial: float = 0.016     # m
    radius_min: float = 0.006         # m
    radius_max: float = 0.018         # m
    shape_bound: float = 0.020        # m, +/- z offset of four control points

    # 6061-T6 aluminum baseline (MatWeb typical values)
    E: float = 68.9e9                 # Pa
    nu: float = 0.33
    rho: float = 2700.0               # kg/m^3
    yield_strength: float = 276e6     # Pa
    yield_fs: float = 2.0

    # Motor / rotor load model (baseline assumptions; replace with team hardware)
    motor_mass: float = 0.099          # kg
    rotor_thrust: float = 12.65        # N
    thrust_tilt_deg: float = 25.0     # deg inward from vertical, creates compression
    motor_torque_z: float = 0.0      # N*m about global z
    gravity: float = 9.80665

    # Constraints
    tip_deflection_limit: float = 0.0030   # m
    frequency_limit: float = 140.0          # Hz
    buckling_factor_limit: float = 3.0
    ks_rho: float = 50.0

    # Optimization
    maxiter: int = 120
    ftol: float = 1e-8

    @property
    def G(self):
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def stress_limit(self):
        return self.yield_strength / self.yield_fs

CFG = Config()
