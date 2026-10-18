"""Independent checks of inverses, Jacobians and conserved column mass."""

import math
import unittest
from coordinate_maps import (
    R_from_r, r_from_R, amplification, Zhat_from_z, z_from_Zhat,
    density, vertical_mass_weight, mass_metric, DRY_COEFFICIENTS,
)


def simpson(f, lower, upper, intervals=2000):
    """Composite Simpson quadrature; intervals must be even."""
    step = (upper - lower) / intervals
    total = f(lower) + f(upper)
    for i in range(1, intervals):
        total += (4 if i % 2 else 2) * f(lower + i * step)
    return total * step / 3


class CoordinateMapTests(unittest.TestCase):
    def test_round_trips(self):
        for r in [0.095, 0.1, 0.5, 1.0, 5.0, 7.0]:
            self.assertTrue(math.isclose(r_from_R(R_from_r(r)), r,
                                        rel_tol=3e-12, abs_tol=1e-14))
        for i in range(101):
            z = i / 100
            self.assertAlmostEqual(z_from_Zhat(Zhat_from_z(z)), z, places=14)

    def test_height_derivative_and_normalisation(self):
        # A numerical derivative independently checks dZhat/dz = A.
        dz = 1e-6
        for z in [0.1, 0.5, 0.9]:
            derivative = (Zhat_from_z(z + dz) - Zhat_from_z(z - dz)) / (2 * dz)
            self.assertAlmostEqual(derivative, amplification(z), places=9)
        self.assertAlmostEqual(simpson(amplification, 0, 1), 1, places=12)

    def test_dry_and_constant_density_limits(self):
        for z in [0, 0.2, 0.7, 1]:
            self.assertEqual(amplification(z, DRY_COEFFICIENTS), 1)
            self.assertEqual(Zhat_from_z(z, DRY_COEFFICIENTS), z)
            self.assertEqual(density(z, H0=math.inf, rho_s=2), 2)
            self.assertEqual(vertical_mass_weight(z, DRY_COEFFICIENTS,
                                                  H0=math.inf, rho_s=2), 2)

    def test_transformed_column_mass(self):
        # Integral W dZhat must equal integral rho0 dz, regardless of A.
        for H0 in [0.75, math.inf]:
            exact = 1 if H0 == math.inf else H0 * (-math.expm1(-1 / H0))
            computed = simpson(lambda u: vertical_mass_weight(u, H0=H0), 0, 1)
            self.assertAlmostEqual(computed, exact, places=10)

    def test_radial_metric_against_physical_annular_mass(self):
        # In the dry constant-density limit, integrate dm/dR across an annulus.
        a, b = 0.3, 0.8
        computed = simpson(lambda R: mass_metric(R, 0.5, H0=math.inf,
                                                coefficients=DRY_COEFFICIENTS),
                           R_from_r(a), R_from_r(b))
        self.assertAlmostEqual(computed, math.pi * (b**2 - a**2), places=8)

    def test_invalid_coordinates_and_nonmonotone_profile(self):
        for call in [lambda: R_from_r(0), lambda: R_from_r(0.05),
                     lambda: r_from_R(1 / (2 * 0.095**2)),
                     lambda: z_from_Zhat(1.1), lambda: density(0.5, H0=0),
                     lambda: amplification(0.5, (1, -8, 8))]:
            with self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main(verbosity=2)
