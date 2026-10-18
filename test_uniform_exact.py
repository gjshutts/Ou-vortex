"""Validate the uniform-density shortcut against general quadrature and saved results."""
import json,math,unittest
from pathlib import Path
from coordinate_maps import DRY_COEFFICIENTS
from mass_integration import polygon_mass
from mass_jacobian import edge_mass_integral


class ExactUniformTests(unittest.TestCase):
    def test_saved_twenty_cell_masses(self):
        d=json.loads(Path(__file__).with_name('equilibrium_dry_20.json').read_text())
        for p,m in zip(d['polygons'],d['masses']):
            exact=polygon_mass(p,r0=.19,coefficients=DRY_COEFFICIENTS,H0=math.inf)
            general=polygon_mass(p,r0=.19,coefficients=DRY_COEFFICIENTS,H0=math.inf,
                                 use_exact_uniform=False)
            self.assertAlmostEqual(exact,m,places=9)
            self.assertAlmostEqual(exact,general,places=9)

    def test_shared_edges_in_both_directions(self):
        for start,end in [((0,0),(10,1)),((40,.2),(40+1e-10,.8)),((52,.8),(10,.1))]:
            for a,b in [(start,end),(end,start)]:
                exact=edge_mass_integral(a,b,coefficients=DRY_COEFFICIENTS,H0=math.inf)
                general=edge_mass_integral(a,b,coefficients=DRY_COEFFICIENTS,H0=math.inf,
                                          use_exact_uniform=False)
                self.assertTrue(math.isclose(exact,general,rel_tol=1e-9,abs_tol=1e-11))

    def test_scaled_constant_amplification(self):
        poly=[(0,0),(20,0),(20,1),(0,1)]
        a=polygon_mass(poly,coefficients=(3,0,0),H0=math.inf,rho_s=2)
        b=polygon_mass(poly,coefficients=DRY_COEFFICIENTS,H0=math.inf,rho_s=2)
        self.assertEqual(a,b)


if __name__=='__main__':
    unittest.main(verbosity=2)
