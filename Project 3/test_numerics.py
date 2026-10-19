"""Scientific checks independent of neural optimizer convergence.

Run python -m unittest -v test_numerics.py from this directory.
"""
import unittest
import numpy as np
import oscillating_bl as bl


class ScientificChecks(unittest.TestCase):
    def test_reference_boundaries_and_pde(self):
        k = .005
        y = np.array([.005, .04, .1, .3])
        t = np.array([.2, 1.3, 2.4, 5.])
        q = (1+1j)/np.sqrt(2*k)
        amp = bl.exact_amplitude(y, k)
        ut = np.real(1j*amp*np.exp(1j*t))
        uyy = np.real(-q*q*np.exp(-q*y)*np.exp(1j*t))
        np.testing.assert_allclose(ut+np.sin(t)-k*uyy, 0, atol=2e-15)
        np.testing.assert_allclose(bl.exact_velocity(0, t, k), 0, atol=1e-15)
        finite = bl.exact_amplitude([0., 1.], k, finite=True)
        np.testing.assert_allclose(finite, [0, 1], atol=1e-15)

    def test_full_parameter_jacobian_by_complex_step(self):
        m = 5
        p = bl.initial_parameters(m, 17, .003)
        data = bl.make_measurements(bl.Case(.1))
        yc = np.linspace(.001, .999, 17)
        for scale in [1., .07]:
            r, jac = bl.residual_jacobian(p, data, yc, m, scale, True)
            numerical = np.zeros_like(jac)
            for j in range(len(p)):
                pc = p.astype(complex)
                pc[j] += 1e-28j
                numerical[:, j] = bl.residual_jacobian(pc, data, yc, m, scale).imag/1e-28
            np.testing.assert_allclose(jac, numerical, rtol=5e-11, atol=1e-11)

    def test_second_spatial_derivative(self):
        p = bl.initial_parameters(6, 29, .001)
        y = np.array([.003, .02, .1, .7])
        for scale in [1., .03]:
            _, _, _, _, a2, b2 = bl.neural_fields(p, y, 6, scale)
            fc = bl.neural_fields(p, y.astype(complex)+1e-25j, 6, scale)
            np.testing.assert_allclose(fc[2].imag/1e-25, a2, rtol=1e-10, atol=1e-10)
            np.testing.assert_allclose(fc[3].imag/1e-25, b2, rtol=1e-10, atol=1e-10)

    def test_fd_second_order_convergence(self):
        k = .005
        errors = []
        for n in [100, 200, 400]:
            y, a = bl.fd_periodic(k, n)
            errors.append(np.linalg.norm(a-bl.exact_amplitude(y, k, finite=True))/np.linalg.norm(a))
        self.assertGreater(errors[0]/errors[1], 3.5)
        self.assertGreater(errors[1]/errors[2], 3.5)

    def test_noiseless_inverse_recovers_viscosity(self):
        case = bl.Case(.1, 0., 16)
        fit = bl.fit_baseline(bl.make_measurements(case), 800)
        self.assertLess(abs(fit['kappa']/case.kappa-1), .001)

    def test_noise_and_probe_coordinates_reproducible(self):
        c = bl.Case(.01, .02, 16, 17)
        d1, d2 = bl.make_measurements(c), bl.make_measurements(c)
        np.testing.assert_array_equal(d1['observed'], d2['observed'])
        np.testing.assert_allclose(np.unique(d1['y']), [.005, .01, .02])
        clean = bl.make_measurements(bl.Case(.01, 0., 16, 17))
        np.testing.assert_array_equal(clean['observed'], clean['clean'])


if __name__ == '__main__':
    unittest.main()
