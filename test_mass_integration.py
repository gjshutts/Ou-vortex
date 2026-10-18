"""Exact mass benchmarks and conservation across polygon partitions."""
import math
import unittest
from coordinate_maps import R_from_r, Zhat_from_z, DRY_COEFFICIENTS
from cell_geometry import rectangle, cell_polygons
from mass_integration import integrate_1d, polygon_mass, cell_masses, annular_column_mass


class MassIntegrationTests(unittest.TestCase):
    def test_quadrature_and_failure_reporting(self):
        self.assertAlmostEqual(integrate_1d(lambda x: x**4, 0, 1), 0.2, places=11)
        with self.assertRaises(ArithmeticError):
            integrate_1d(math.exp, 0, 1, absolute_tolerance=1e-15, max_depth=0)

    def test_physical_rectangles_in_all_three_limits(self):
        for coefficients in [DRY_COEFFICIENTS, (1, 0.548584, 4.273361)]:
            for H0 in [math.inf, 0.75]:
                for a, b, zb, zt in [(0.095, 7, 0, 1), (0.5, 2, 0.2, 0.8)]:
                    poly = rectangle(R_from_r(a), R_from_r(b),
                                     Zhat_from_z(zb, coefficients), Zhat_from_z(zt, coefficients))
                    actual = polygon_mass(poly, coefficients=coefficients, H0=H0, rho_s=1.7)
                    exact = annular_column_mass(a, b, zb, zt, H0=H0, rho_s=1.7)
                    self.assertAlmostEqual(actual, exact, places=8)

    def test_sloping_triangle_has_exact_dry_mass(self):
        Rinf = 1 / (2 * 0.095**2)
        right = 20.0
        triangle = [(0, 0), (right, 0), (0, 1)]
        exact = math.pi / (2 * right) * (-math.log1p(-right / Rinf)) - math.pi / (2 * Rinf)
        computed = polygon_mass(triangle, coefficients=DRY_COEFFICIENTS, H0=math.inf)
        self.assertAlmostEqual(computed, exact, places=11)

    def test_partition_conserves_mass(self):
        domain = rectangle(0, 30)
        planes = [(0, 0, 0), (0.05, 0, -0.5), (0, 1, -0.5)]
        cells = cell_polygons(planes, domain)
        for coefficients in [DRY_COEFFICIENTS, (1, 0.548584, 4.273361)]:
            computed = math.fsum(cell_masses(cells, coefficients=coefficients))
            exact = annular_column_mass(0.095, 1 / math.sqrt(1 / 0.095**2 - 60))
            self.assertAlmostEqual(computed, exact, places=9)

    def test_orientation_and_empty_cells(self):
        poly = rectangle(0, 10)
        self.assertAlmostEqual(polygon_mass(poly), polygon_mass(list(reversed(poly))), places=12)
        self.assertEqual(polygon_mass([]), 0)
        self.assertEqual(polygon_mass([(0, 0), (1, 0), (2, 0)]), 0)

    def test_invalid_parameters(self):
        for call in [lambda: polygon_mass(rectangle(0, 60)),
                     lambda: polygon_mass(rectangle(0, 1, 0, 2)),
                     lambda: polygon_mass([], H0=-1),
                     lambda: polygon_mass([], absolute_tolerance=0)]:
            with self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main(verbosity=2)
