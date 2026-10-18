"""Small exact equilibria, step control, status reporting and initial data."""
import math
import unittest
from coordinate_maps import R_from_r, r_from_R, DRY_COEFFICIENTS
from cell_geometry import rectangle
from mass_integration import annular_column_mass
from equilibrium_solver import solve_offsets, solve_linear_system, solve_zero_pv_equilibrium
from initial_state import generate_initial_state


class EquilibriumSolverTests(unittest.TestCase):
    def exact_two_cell_problem(self, boundary, start, **options):
        r0 = .095
        domain = rectangle(0,20)
        targets = [annular_column_mass(r0,r_from_R(boundary,r0),H0=math.inf),
                   annular_column_mass(r_from_R(boundary,r0),r_from_R(20,r0),H0=math.inf)]
        return solve_offsets([0,1],[0,0],targets,domain,[0,-start],r0=r0,
                             coefficients=DRY_COEFFICIENTS,H0=math.inf,verbose=False, **options)

    def test_linear_solver_pivoting_and_singular_detection(self):
        solution = solve_linear_system([[0,2],[1,3]],[4,7])
        self.assertEqual(solution,[1,2])
        with self.assertRaises(ArithmeticError):
            solve_linear_system([[1,2],[2,4]],[1,2])

    def test_known_two_cell_equilibrium(self):
        result = self.exact_two_cell_problem(10,9,damping=1,mass_tolerance=1e-10)
        self.assertTrue(result["converged"])
        self.assertAlmostEqual(result["offsets"][1],-10,places=8)
        self.assertEqual(result["offsets"][0],0)
        self.assertLess(abs(result["total_mass"]-result["total_target_mass"]),1e-9)

    def test_backtracking_prevents_disappearance(self):
        result = self.exact_two_cell_problem(19,.1,damping=1,mass_tolerance=1e-9)
        self.assertTrue(result["converged"])
        self.assertFalse(result["empty_cells"])
        self.assertTrue(any(h.get("accepted_step",1)<1 for h in result["history"]))
        errors = [h["max_mass_error"] for h in result["history"]]
        self.assertTrue(all(b<a for a,b in zip(errors,errors[1:])))

    def test_iteration_limit_is_not_success(self):
        result = self.exact_two_cell_problem(10,9,max_iterations=0)
        self.assertFalse(result["converged"])
        self.assertEqual(result["status"],"iteration_limit")

    def test_initial_labels_and_total_mass(self):
        data = generate_initial_state(4,.19,3,2,1)
        for i,q in enumerate(data["Q"]):
            center = (data["r_edges"][i]+data["r_edges"][i+1])/2
            self.assertEqual(q,center**4/4)
            self.assertEqual(data["bm"][i],2*math.exp(-q))
        self.assertAlmostEqual(math.fsum(data["target_masses"]),annular_column_mass(.19,3),places=12)

    def test_small_moist_case(self):
        result = solve_zero_pv_equilibrium(4,.19,3,2,1,verbose=False)
        self.assertTrue(result["converged"])
        self.assertLess(result["max_mass_error"],1e-6)
        self.assertEqual(result["empty_cells"],[])

    def test_inconsistent_targets_are_rejected(self):
        with self.assertRaises(ValueError):
            solve_offsets([0,1],[0,0],[1,1],rectangle(0,1),[0,-.5],verbose=False)


if __name__ == "__main__":
    unittest.main(verbosity=2)
