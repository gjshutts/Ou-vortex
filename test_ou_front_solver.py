"""Analytical limits and numerical checks, independent of MATLAB execution."""

import math
import unittest

from ou_front_solver import (
    _rk45, b_Q, branch_derivative, front_rhs, integrate_front,
    jump_slope, launch_front, lid_radius, psi_of_Q, solve_Qt,
    solve_zero_pv_vortex_front, tip_residual,
)


class OuFrontTests(unittest.TestCase):
    def test_finite_jump_and_infinitesimal_limit(self):
        qm, qp, alpha, b0 = 1., 2., 1., 7.6
        S = jump_slope(qm, qp, alpha, b0)
        self.assertAlmostEqual(S * (b0*math.exp(-qp)-b0*math.exp(-qm)), -(qp-qm))
        self.assertAlmostEqual(jump_slope(qm, qm+1e-8, alpha, b0),
                               -1/b_Q(qm, alpha, b0), delta=2e-9)

    def test_dry_weak_buoyancy_lid_limit(self):
        # A surface in a fluid initially at rest: r=(4 Q)**(1/4).
        self.assertAlmostEqual(lid_radius(3., 1., 1e-14), (4*3.)**.25, places=12)

    def test_tip_equation(self):
        for b0 in (2., 7.6, 16.):
            Qt = solve_Qt(1., b0)
            self.assertLess(abs(tip_residual(Qt, 1., b0)), 2e-10)
            self.assertTrue(-1 < psi_of_Q(Qt, 1., b0) < 0)

    def test_launch_branches_and_reference_geometry(self):
        launch = launch_front(1., 7.6, .6965)
        self.assertLess(launch['state'][0], launch['Qt'])
        self.assertGreater(launch['state'][1], launch['Qt'])
        S = -1 / b_Q(launch['Qt'], 1., 7.6)
        Rt = -1/(2*launch['rf_t']**2)
        for z, qm, qp, r in launch['warmup']:
            self.assertAlmostEqual(-1/(2*r*r)-Rt, (z-launch['Zt'])/S, places=14)
            self.assertTrue(qm < qp)

    def test_rk45_against_exact_exponential(self):
        rows, _ = _rk45(lambda z, y: tuple(-v for v in y),
                        0., (1., 2., 3.), 1e-10, 1e-12)
        self.assertEqual(rows[-1][0], 1.)
        for actual, initial in zip(rows[-1][1:], (1., 2., 3.)):
            self.assertAlmostEqual(actual, initial/math.e, delta=3e-11)

    def test_rhs_has_margules_radius_slope(self):
        launch = launch_front(1., 7.6, .6965)
        qm, qp, r = launch['state']
        dz = front_rhs(launch['z_start'], launch['state'], 1., 7.6)
        # R_r=r^-3, hence dR/dz=(dr/dz)/r^3=1/S.
        self.assertAlmostEqual(dz[2]/r**3, 1/jump_slope(qm, qp, 1., 7.6))

    def test_invalid_or_incomplete_shots_are_rejected(self):
        with self.assertRaises(ValueError):
            launch_front(1., 7.6, .99999)
        with self.assertRaises(ArithmeticError):
            integrate_front(1., 7.6, .6)
        with self.assertRaises(ValueError):
            solve_zero_pv_vortex_front(1., 2., .85, .9, verbose=False)

    def test_shooting_and_tolerance_refinement(self):
        results = []
        for rtol, atol in ((1e-8, 1e-10), (1e-10, 1e-12)):
            shot, sol = solve_zero_pv_vortex_front(
                1., 7.6, .68, .70, verbose=False, rtol=rtol, atol=atol,
                tip_tolerance=1e-11)
            self.assertEqual(sol['z'][-1], 1.)
            self.assertLess(abs(shot['residual']), 1e-8)
            self.assertLess(abs(shot['front_lid_radius']-shot['rlid_p']), 3e-7)
            self.assertTrue(all(a < b for a, b in zip(sol['Qm'], sol['Qp'])))
            results.append(shot)
        self.assertLess(abs(results[0]['Zt']-results[1]['Zt']), 2e-8)
        self.assertLess(abs(results[0]['front_lid_radius']-results[1]['front_lid_radius']), 2e-7)


if __name__ == '__main__':
    unittest.main()
