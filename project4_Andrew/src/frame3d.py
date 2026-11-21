import numpy as np


def section_props(r, t):
    ri = r - t
    if np.any(np.asarray(ri) <= 0):
        raise ValueError('Wall thickness must be smaller than outer radius.')
    A = np.pi * (r**2 - ri**2)
    I = 0.25 * np.pi * (r**4 - ri**4)
    J = 2.0 * I
    return A, I, I, J


def local_axes(x1, x2):
    dx = np.asarray(x2, dtype=float) - np.asarray(x1, dtype=float)
    L = np.linalg.norm(dx)
    ex = dx / L
    # Choose a reference direction that is not nearly parallel to ex.
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ex, ref)) > 0.90:
        ref = np.array([0.0, 1.0, 0.0])
    ey = np.cross(ref, ex)
    ey /= np.linalg.norm(ey)
    ez = np.cross(ex, ey)
    R = np.vstack([ex, ey, ez])  # global vector -> local vector
    return L, R


def transform12(R):
    T = np.zeros((12, 12))
    for i in range(4):
        T[3*i:3*i+3, 3*i:3*i+3] = R
    return T


def local_stiffness(E, G, A, Iy, Iz, J, L):
    k = np.zeros((12, 12))
    EA = E*A/L
    GJ = G*J/L
    EIy = E*Iy
    EIz = E*Iz

    # Axial
    k[0,0] = k[6,6] = EA
    k[0,6] = k[6,0] = -EA
    # Torsion
    k[3,3] = k[9,9] = GJ
    k[3,9] = k[9,3] = -GJ

    # Bending: local v with rz, about z
    ids = [1,5,7,11]
    kb = EIz / L**3 * np.array([
        [12, 6*L, -12, 6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2],
    ])
    k[np.ix_(ids, ids)] += kb

    # Bending: local w with ry, about y. Sign follows ry convention.
    ids = [2,4,8,10]
    kb = EIy / L**3 * np.array([
        [12, -6*L, -12, -6*L],
        [-6*L, 4*L**2, 6*L, 2*L**2],
        [-12, 6*L, 12, 6*L],
        [-6*L, 2*L**2, 6*L, 4*L**2],
    ])
    k[np.ix_(ids, ids)] += kb
    return k


def local_mass(rho, A, J, L):
    m = np.zeros((12, 12))
    mt = rho*A*L

    # Axial consistent mass
    ma = mt/6.0 * np.array([[2.0,1.0],[1.0,2.0]])
    m[np.ix_([0,6],[0,6])] += ma

    # Torsional rotary inertia
    mr = rho*J*L/6.0 * np.array([[2.0,1.0],[1.0,2.0]])
    m[np.ix_([3,9],[3,9])] += mr

    # Euler-Bernoulli consistent bending mass, each plane
    mb = mt/420.0 * np.array([
        [156, 22*L, 54, -13*L],
        [22*L, 4*L**2, 13*L, -3*L**2],
        [54, 13*L, 156, -22*L],
        [-13*L, -3*L**2, -22*L, 4*L**2],
    ])
    m[np.ix_([1,5,7,11],[1,5,7,11])] += mb

    # w/ry has opposite slope sign; apply a sign transformation.
    S = np.diag([1.0,-1.0,1.0,-1.0])
    mbw = S @ mb @ S
    m[np.ix_([2,4,8,10],[2,4,8,10])] += mbw
    return m


def local_geometric_stiffness(N, L):
    # N > 0 tension; N < 0 compression. This matrix is added to K.
    kg = np.zeros((12,12))
    base = N/(30.0*L) * np.array([
        [36, 3*L, -36, 3*L],
        [3*L, 4*L**2, -3*L, -L**2],
        [-36, -3*L, 36, -3*L],
        [3*L, -L**2, -3*L, 4*L**2],
    ])
    kg[np.ix_([1,5,7,11],[1,5,7,11])] += base
    S = np.diag([1.0,-1.0,1.0,-1.0])
    kg[np.ix_([2,4,8,10],[2,4,8,10])] += S @ base @ S
    return kg


def element_matrices(x1, x2, r, t, E, G, rho):
    A, Iy, Iz, J = section_props(r, t)
    L, R = local_axes(x1, x2)
    T = transform12(R)
    kl = local_stiffness(E,G,A,Iy,Iz,J,L)
    ml = local_mass(rho,A,J,L)
    return T.T @ kl @ T, T.T @ ml @ T, T, kl, L, (A,Iy,Iz,J)
