import math
import unittest
from cell_geometry import rectangle, cell_polygons
from coordinate_maps import DRY_COEFFICIENTS
from mass_jacobian import mass_jacobian, edge_mass_integral, finite_difference_jacobian


class MassJacobianTests(unittest.TestCase):
    def setUp(self):
        self.domain = rectangle(0, 30)
        self.planes = [(0, 0, 0), (0.05, 0, -0.5), (0, 1, -0.5)]
        self.cells = cell_polygons(self.planes, self.domain)

    def test_exact_two_cell_interface(self):
        planes = [(0, 0, 0), (1, 0, -10)]
        domain = rectangle(0, 20)
        J = mass_jacobian(planes, cell_polygons(planes, domain),
                          coefficients=DRY_COEFFICIENTS, H0=math.inf)
        Rinf = 1 / (2 * 0.095**2)
        exact = math.pi / (2 * (Rinf - 10)**2)
        self.assertAlmostEqual(J[0][0], exact, places=12)
        self.assertAlmostEqual(J[0][1], -exact, places=12)

    def test_symmetry_signs_null_direction_and_positive_quadratic_form(self):
        J = mass_jacobian(self.planes, self.cells)
        for i in range(3):
            self.assertAlmostEqual(math.fsum(J[i]), 0, places=12)
            for j in range(3):
                self.assertEqual(J[i][j], J[j][i])
                if i != j:
                    self.assertLessEqual(J[i][j], 0)
        for vector in [(1, -2, 3), (2, 1, -1)]:
            form = math.fsum(vector[i]*J[i][j]*vector[j] for i in range(3) for j in range(3))
            self.assertGreaterEqual(form, 0)

    def test_finite_differences_in_dry_density_and_moist_cases(self):
        for coefficients, H0 in [(DRY_COEFFICIENTS, math.inf),
                                 (DRY_COEFFICIENTS, 0.75), ((1, .548584, 4.273361), .75)]:
            parameters = dict(coefficients=coefficients, H0=H0)
            J = mass_jacobian(self.planes, self.cells, **parameters)
            fd = finite_difference_jacobian(self.planes, self.domain, step=1e-5,
                                            absolute_tolerance=1e-12, **parameters)
            error = max(abs(J[i][j]-fd[i][j]) for i in range(3) for j in range(3))
            scale = max(abs(v) for row in J for v in row)
            self.assertLess(error / scale, 1e-6)

    def test_common_offset_shift(self):
        planes = [(0, 0, 0), (1, 0, -0.5), (0, 1, -0.5)]
        shifted = [(Q, bm, S+1000) for Q,bm,S in planes]
        domain = rectangle(0, 1)
        self.assertEqual(mass_jacobian(planes, cell_polygons(planes, domain)),
                         mass_jacobian(shifted, cell_polygons(shifted, domain)))

    def test_horizontal_edge_and_reversal(self):
        start, end = (0, .5), (10, .5)
        value = edge_mass_integral(start, end, coefficients=DRY_COEFFICIENTS, H0=math.inf)
        Rinf = 1 / (2 * .095**2)
        exact = math.pi / 2 * (1/(Rinf-10)-1/Rinf)
        self.assertGreater(value, 0)
        self.assertAlmostEqual(value, exact, places=11)
        self.assertAlmostEqual(value, edge_mass_integral(end, start,
                              coefficients=DRY_COEFFICIENTS, H0=math.inf), places=12)

    def test_relative_quadrature_for_large_integral(self):
        from mass_integration import integrate_1d
        computed = integrate_1d(lambda x: 1e12*math.exp(x), 0, 1,
                                absolute_tolerance=1e-10, relative_tolerance=1e-9)
        exact = 1e12*math.expm1(1)
        self.assertLess(abs(computed-exact)/exact, 1e-9)

    def test_empty_cells_and_parallel_planes(self):
        planes = [(1, 2, 0), (1, 2, 1)]
        J = mass_jacobian(planes, cell_polygons(planes, self.domain))
        self.assertEqual(J, [[0, 0], [0, 0]])

    def test_duplicate_planes_and_invalid_inputs(self):
        with self.assertRaises(ValueError):
            mass_jacobian([(1,2,0), (1,2,0)], [[], []])
        with self.assertRaises(ValueError):
            mass_jacobian(self.planes, self.cells, edge_tolerance=0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
