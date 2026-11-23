"""Independent limits, FE comparisons and derivative checks for NACA modes."""
from dataclasses import replace
import unittest

import numpy as np

from arm import Beam
from vibration import (characteristic, first_frequency_partials, tip_mass_frequencies,
                       tip_mass_modes, tip_mass_roots)


class TipMassVibrationTests(unittest.TestCase):
    def setUp(self):
        self.beam=Beam(1.,70e9,.4,.8,100.,elements=24)
        self.D,self.t=.04,.002

    def test_no_tip_mass_recovers_classical_cantilever_roots(self):
        known=[1.875104068711961,4.694091132974175,7.854757438237612,10.995540734875467]
        np.testing.assert_allclose(tip_mass_roots(0.,4),known,rtol=1e-13)

    def test_four_modes_match_independent_fe_across_mass_ratios(self):
        for ratio in [0.,.01,.5,2.,10.]:
            beam=replace(self.beam,tip_mass=ratio*self.beam.mass_per_length*self.beam.length)
            analytical=beam.modes(self.D,self.t,4)[0]
            fe=beam.fe_modes(self.D,self.t,4)[0]
            np.testing.assert_allclose(analytical,fe,rtol=4e-5)

    def test_large_tip_mass_recovers_massless_beam_spring_limit(self):
        L,EI,mu,mt=1.3,2500.,.5,1e8
        expected=np.sqrt(3*EI/(mt*L**3))/(2*np.pi)
        actual=tip_mass_frequencies(L,EI,mu,mt,1)[0]
        self.assertAlmostEqual(actual/expected,1.,delta=2e-8)

    def test_exact_shapes_are_orthonormal_with_tip_mass_and_clamped(self):
        L,EI,mu,mt=1.3,2500.,.5,.8
        x,w=np.polynomial.legendre.leggauss(96)
        x=np.r_[0.,(x+1)*L/2,L]
        f,modes=tip_mass_modes(L,EI,mu,mt,x,4)
        y=modes[::2]
        gram=mu*L/2*(y[1:-1].T@(w[:,None]*y[1:-1]))+mt*np.outer(y[-1],y[-1])
        np.testing.assert_allclose(gram,np.eye(4),atol=3e-12)
        np.testing.assert_array_equal(modes[:2],0.)
        self.assertTrue(np.all(y[-1]>0))
        # Independent complex-step derivative of the printed equation (14).
        ratio=mt/(mu*L)
        for j,beta in enumerate(tip_mass_roots(ratio,4)):
            a=(np.sinh(beta)+np.sin(beta))/(np.cosh(beta)+np.cos(beta))
            def raw(s):
                z=beta*s/L
                return np.sinh(z)-np.sin(z)+a*(np.cos(z)-np.cosh(z))
            scaling=y[-1,j]/raw(L)
            slope=raw(x.astype(complex)+1e-25j).imag/1e-25*scaling
            np.testing.assert_allclose(modes[1::2,j],slope,atol=3e-10)
            self.assertLess(abs(characteristic(beta,ratio)),1e-11)

    def test_implicit_root_partials_match_independent_resolves(self):
        for mt in [.02,.8,30.]:
            args=np.array([1.3,2500.,.5,mt])
            f,partials=first_frequency_partials(*args)
            self.assertAlmostEqual(f,tip_mass_frequencies(*args,1)[0])
            for j,name in enumerate(["length","EI","mass_per_length","tip_mass"]):
                h=args[j]*1e-4
                delta=np.eye(4)[j]*h
                fd=(tip_mass_frequencies(*(args+delta),1)[0]-tip_mass_frequencies(*(args-delta),1)[0])/(2*h)
                self.assertAlmostEqual(partials[name]/fd,1.,delta=2e-6)

    def test_frequency_is_mesh_independent_and_drops_with_tip_mass(self):
        coarse=replace(self.beam,elements=2).modes(self.D,self.t,4)[0]
        fine=replace(self.beam,elements=80).modes(self.D,self.t,4)[0]
        np.testing.assert_array_equal(coarse,fine)
        light=replace(self.beam,tip_mass=0.).modes(self.D,self.t,4)[0]
        self.assertTrue(np.all(coarse<light))

    def test_invalid_inputs_and_custom_supports_are_explicit(self):
        for args in [(1.,1.,0.,1.),(0.,1.,1.,1.),(1.,1.,1.,-1.)]:
            with self.assertRaises(ValueError):tip_mass_frequencies(*args)
        with self.assertRaises(ValueError):tip_mass_roots(0.,0)
        with self.assertRaises(ValueError):self.beam.modes(self.D,self.t,free=[2,3])
        np.testing.assert_array_equal(self.beam.modes(self.D,self.t,method="fem")[0],
                                      self.beam.fe_modes(self.D,self.t)[0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
