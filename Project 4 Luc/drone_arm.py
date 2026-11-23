"""Linear static and free-vibration FEM for a clamped, straight circular tube.

The long axis is +x and +z is up. Each node stores [ux, uy, uz, rx, ry, rz].
All solver inputs and outputs use SI units. The CLI accepts tube sizes in inches.
"""

import argparse
from dataclasses import dataclass
from math import isfinite, pi

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh, spsolve

import arm_inputs as inputs


INCH = 0.0254
LBM = 0.45359237
STANDARD_GRAVITY = 9.80665
LBF = LBM * STANDARD_GRAVITY
DOFS_PER_NODE = 6


class NondifferentiableError(ValueError):
    """The selected structural response has no unique sizing gradient."""


@dataclass(frozen=True)
class Arm:
    """Physical parameters in SI. Thrust is the force of one of two motors."""

    length: float = inputs.LENGTH_IN * INCH
    outer_diameter: float = inputs.INITIAL_OD_IN * INCH
    wall: float = inputs.INITIAL_WALL_IN * INCH
    young: float = inputs.YOUNG_MODULUS_GPA * 1e9
    poisson: float = inputs.POISSON_RATIO
    density: float = inputs.DENSITY_KG_M3
    tip_mass: float = inputs.MOTOR_MASS_COMBINED_LBM * LBM
    thrust_per_motor: float = inputs.THRUST_PER_MOTOR_LBF * LBF
    rpm: float = inputs.ROTOR_RPM
    elements: int = inputs.BEAM_ELEMENTS
    gravity: float = inputs.GRAVITY_M_S2

    def __post_init__(self):
        positive = (self.length, self.outer_diameter, self.wall, self.young, self.density)
        if any(not isfinite(value) or value <= 0 for value in positive):
            raise ValueError("Length, diameter, wall, Young's modulus and density must be positive and finite.")
        if 2 * self.wall >= self.outer_diameter:
            raise ValueError("Wall thickness must be smaller than the outer radius.")
        if not -1.0 < self.poisson < 0.5:
            raise ValueError("Poisson's ratio must lie between -1 and 0.5.")
        nonnegative = (self.tip_mass, self.thrust_per_motor, self.rpm, self.gravity)
        if any(not isfinite(value) or value < 0 for value in nonnegative):
            raise ValueError("Tip mass, thrust, rpm and gravity must be nonnegative and finite.")
        if isinstance(self.elements, bool) or not isinstance(self.elements, (int, np.integer)) or self.elements < 1:
            raise ValueError("The number of beam elements must be a positive integer.")

    @property
    def inner_diameter(self):
        return self.outer_diameter - 2 * self.wall

    @property
    def area(self):
        return pi * self.wall * (self.outer_diameter - self.wall)

    @property
    def inertia(self):
        return self.area * (self.outer_diameter**2 + self.inner_diameter**2) / 16

    def section_derivatives(self):
        """Return area and bending inertia derivatives with respect to OD and wall."""
        d, inner = self.outer_diameter, self.inner_diameter
        d_area = pi * np.array([self.wall, inner])
        # Differentiate the factored inertia expression to retain thin-wall accuracy.
        d_inertia = (d_area * (d**2 + inner**2)
                     + self.area * np.array([2 * (d + inner), -4 * inner])) / 16
        return d_area, d_inertia

    @property
    def polar_inertia(self):
        return 2 * self.inertia

    @property
    def tube_mass(self):
        return self.density * self.area * self.length


@dataclass
class Model:
    arm: Arm
    x: np.ndarray
    stiffness: csr_matrix
    mass: csr_matrix
    load: np.ndarray
    element_stiffness: np.ndarray
    element_mass: np.ndarray
    element_load: np.ndarray


def element_matrices(arm):
    """Two-node frame with axial extension, twist and bending in two planes."""
    h = arm.length / arm.elements
    stiffness = np.zeros((12, 12))
    mass = np.zeros((12, 12))
    shear_modulus = arm.young / (2 * (1 + arm.poisson))

    # Linear interpolation of axial displacement and twist.
    linear_stiffness = np.array([[1, -1], [-1, 1]]) / h
    linear_mass = h / 6 * np.array([[2, 1], [1, 2]])
    axial = np.ix_([0, 6], [0, 6])
    twist = np.ix_([3, 9], [3, 9])
    stiffness[axial] = arm.young * arm.area * linear_stiffness
    mass[axial] = arm.density * arm.area * linear_mass
    stiffness[twist] = shear_modulus * arm.polar_inertia * linear_stiffness
    mass[twist] = arm.density * arm.polar_inertia * linear_mass

    # Cubic Hermite interpolation in displacement/slope coordinates.
    bending_stiffness = arm.young * arm.inertia / h**3 * np.array([
        [12, 6*h, -12, 6*h],
        [6*h, 4*h**2, -6*h, 2*h**2],
        [-12, -6*h, 12, -6*h],
        [6*h, 2*h**2, -6*h, 4*h**2],
    ])
    bending_mass = arm.density * arm.area * h / 420 * np.array([
        [156, 22*h, 54, -13*h],
        [22*h, 4*h**2, 13*h, -3*h**2],
        [54, 13*h, 156, -22*h],
        [-13*h, -3*h**2, -22*h, 4*h**2],
    ])

    # Right-hand-rule rotations give uy' = rz and uz' = -ry.
    lateral = np.ix_([1, 5, 7, 11], [1, 5, 7, 11])
    vertical = np.ix_([2, 4, 8, 10], [2, 4, 8, 10])
    slope_sign = np.array([1, -1, 1, -1])
    sign_product = np.outer(slope_sign, slope_sign)
    stiffness[lateral] = bending_stiffness
    mass[lateral] = bending_mass
    stiffness[vertical] = bending_stiffness * sign_product
    mass[vertical] = bending_mass * sign_product

    q = -arm.density * arm.area * arm.gravity
    element_load = np.zeros(12)
    element_load[[2, 4, 8, 10]] = q * np.array([h/2, -h**2/12, h/2, h**2/12])
    return stiffness, mass, element_load


def element_partials(model):
    """Exact partials of the element operators for the two section dimensions."""
    arm = model.arm
    d_area, d_inertia = arm.section_derivatives()
    area_rate = d_area / arm.area
    inertia_rate = d_inertia / arm.inertia

    # The axial stiffness block is proportional to A. Both bending stiffness
    # blocks and the twist stiffness block are proportional to I (J = 2 I).
    d_stiffness = inertia_rate[:, None, None] * model.element_stiffness
    axial_rows = [0, 6]
    d_stiffness[:, axial_rows] = area_rate[:, None, None] * model.element_stiffness[axial_rows]
    # All mass blocks scale with A except the twist block, proportional to J.
    d_mass = area_rate[:, None, None] * model.element_mass
    twist_rows = [3, 9]
    d_mass[:, twist_rows] = inertia_rate[:, None, None] * model.element_mass[twist_rows]
    d_load = area_rate[:, None] * model.element_load
    return d_area, inertia_rate, d_stiffness, d_mass, d_load


def build_model(arm=Arm()):
    """Assemble sparse global matrices, consistent tube weight and tip loads."""
    element_stiffness, element_mass, element_load = element_matrices(arm)
    x = np.linspace(0.0, arm.length, arm.elements + 1)
    size = DOFS_PER_NODE * len(x)
    starts = DOFS_PER_NODE * np.arange(arm.elements)
    dofs = starts[:, None] + np.arange(12)

    # Shared nodes have duplicate COO indices, summed during conversion to CSR.
    rows = np.repeat(dofs, 12, axis=1).ravel()
    cols = np.tile(dofs, (1, 12)).ravel()

    def assemble(element_matrix):
        values = np.tile(element_matrix.ravel(), arm.elements)
        matrix = coo_matrix((values, (rows, cols)), shape=(size, size)).tocsr()
        matrix.eliminate_zeros()
        return matrix

    stiffness = assemble(element_stiffness)
    mass = assemble(element_mass)
    tip_translations = size - DOFS_PER_NODE + np.arange(3)
    tip_values = np.full(3, arm.tip_mass)
    mass += coo_matrix((tip_values, (tip_translations, tip_translations)), shape=mass.shape).tocsr()

    load = np.zeros(size)
    for element_dofs in dofs:
        load[element_dofs] += element_load
    load[tip_translations[2]] += 2 * arm.thrust_per_motor - arm.tip_mass * arm.gravity
    return Model(arm, x, stiffness, mass, load, element_stiffness, element_mass, element_load)


def statics(model, load=None):
    """Solve with the first six DOFs clamped. Return full u and residual K u - F.

    A supplied load replaces the entire model load. The residual at the first six
    entries gives support reactions. The remaining entries are equilibrium errors.
    """
    load = model.load if load is None else np.asarray(load, dtype=float)
    if load.shape != model.load.shape:
        raise ValueError("The applied load must contain one value per global DOF.")
    free = slice(DOFS_PER_NODE, None)
    displacement = np.zeros_like(load)
    displacement[free] = spsolve(model.stiffness[free, free], load[free])
    return displacement, model.stiffness @ displacement - load


def end_forces(model, displacement, element_load=None):
    """Element resisting end actions in [Fx,Fy,Fz,Mx,My,Mz] order at each end.

    Defaults to the model's tube gravity. Pass zeros(12) for point-load-only
    tests or an array of shape (elements, 12) for other distributed loads.
    """
    element_load = model.element_load if element_load is None else np.asarray(element_load)
    dofs = DOFS_PER_NODE * np.arange(model.arm.elements)[:, None] + np.arange(12)
    return displacement[dofs] @ model.element_stiffness.T - element_load


def bending_stress(model, displacement, element_load=None):
    """Absolute outer-fiber bending stress at each element's two end sections."""
    actions = end_forces(model, displacement, element_load).reshape(-1, 2, 6)
    moments = np.hypot(actions[:, :, 4], actions[:, :, 5])
    return moments * model.arm.outer_diameter / (2 * model.arm.inertia)


def modes(model, count=6):
    """Lowest undamped frequencies [Hz] and full, mass-normalized mode columns."""
    if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 1:
        raise ValueError("Mode count must be a positive integer.")
    free = slice(DOFS_PER_NODE, None)
    stiffness = model.stiffness[free, free]
    mass = model.mass[free, free]
    count = min(count, mass.shape[0])
    if count == mass.shape[0]:
        # ARPACK requires k < dimension; e.g. all six modes of a one-element mesh.
        values, vectors = eigh(stiffness.toarray(), mass.toarray())
    else:
        start = np.random.default_rng(inputs.EIGENSOLVER_SEED).normal(size=mass.shape[0])
        values, vectors = eigsh(stiffness, k=count, M=mass, sigma=0.0, which="LM", v0=start)
    order = np.argsort(values)
    values = values[order]
    if np.any(values <= 0):
        raise ValueError("Nonpositive eigenvalue found in the clamped frame.")
    full_modes = np.zeros((model.mass.shape[0], count))
    full_modes[free] = vectors[:, order]
    return np.sqrt(values) / (2 * pi), full_modes


def tip_stiffness(model):
    """Vertical incremental point stiffness [N/m] from a 1 N tip test."""
    load = np.zeros_like(model.load)
    tip_z = -4
    load[tip_z] = 1.0
    displacement, _ = statics(model, load)
    return 1.0 / displacement[tip_z]


def main():
    defaults = Arm()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--od-in", type=float, default=defaults.outer_diameter / INCH,
                        help="tube outer diameter in inches")
    parser.add_argument("--wall-in", type=float, default=defaults.wall / INCH,
                        help="tube wall thickness in inches")
    parser.add_argument("--elements", type=int, default=defaults.elements,
                        help="number of equal-length frame elements")
    args = parser.parse_args()
    try:
        arm = Arm(outer_diameter=args.od_in * INCH, wall=args.wall_in * INCH, elements=args.elements)
    except ValueError as exc:
        parser.error(str(exc))
    model = build_model(arm)
    displacement, reaction = statics(model)
    stress = bending_stress(model, displacement)
    frequencies, _ = modes(model)
    gross_thrust = 2 * arm.thrust_per_motor
    print(f"Tube OD x wall = {args.od_in:g} x {args.wall_in:g} in, L = {arm.length:.4f} m")
    print(f"Mesh = {arm.elements} elements, {len(model.x)} nodes")
    print(f"Tube mass = {arm.tube_mass:.6f} kg, combined motor mass = {arm.tip_mass:.6f} kg")
    print(f"Gross thrust = {gross_thrust:.6f} N, net tip force (+z) = {gross_thrust - arm.tip_mass * arm.gravity:.6f} N")
    print(f"Distributed tube load (+z) = {-arm.density * arm.area * arm.gravity:.6f} N/m")
    print(f"Root reaction Fz = {reaction[2]:.6f} N, My = {reaction[4]:.6f} N m")
    print(f"Tip displacement uz = {1e3 * displacement[-4]:.6f} mm")
    print(f"Tip slope duz/dx = {-displacement[-2]:.8f} rad")
    print(f"Peak nominal bending stress = {stress.max() / 1e6:.6f} MPa")
    print(f"Vertical tip stiffness = {tip_stiffness(model):.6f} N/m")
    print("Lowest natural frequencies (Hz) = " + ", ".join(f"{value:.4f}" for value in frequencies))
    rotor_hz = arm.rpm / 60
    print(f"Rotor 1x frequency = {rotor_hz:.4f} Hz")
    if rotor_hz > 0:
        print(f"f1 / rotor 1x = {frequencies[0] / rotor_hz:.4f}")


if __name__ == "__main__":
    main()
