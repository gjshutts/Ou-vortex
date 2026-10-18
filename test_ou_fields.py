"""Checks for analytical Q-surface and plotting-grid reconstruction."""
import math
import unittest
from ou_fields import Pchip,OuFields
from ou_front_solver import solve_zero_pv_vortex_front,b_Q


class InterpolationTests(unittest.TestCase):
    def test_linear_data_are_reproduced(self):
        p=Pchip([0.,1.,3.,5.],[2.,5.,11.,17.])
        for x in (.2,.8,1.5,3.8):
            self.assertAlmostEqual(p(x),2+3*x,places=13)

    def test_no_overshoot_on_monotone_data(self):
        p=Pchip([0.,.2,1.,3.],[0.,.01,.8,1.])
        values=[p(3*i/1000) for i in range(1001)]
        self.assertTrue(all(a<=b for a,b in zip(values,values[1:])))
        with self.assertRaises(ValueError):p(4.)


class FieldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _,sol=solve_zero_pv_vortex_front(1.,7.,.68,.72,verbose=False,
                                       rtol=1e-10,atol=1e-12,tip_tolerance=1e-11)
        cls.f=OuFields(sol)

    def test_true_tip_and_lid_are_used(self):
        z,r=self.f.top(self.f.Qt)
        self.assertEqual(z,self.f.sol['Zt'])
        self.assertEqual(r,self.f.sol['rf_t'])
        self.assertEqual(self.f.top(self.f.qm)[0],1.)
        self.assertEqual(self.f.top(self.f.qp)[0],1.)

    def test_surface_balance_and_termination(self):
        for Q in (self.f.qm/2,.7,self.f.Qt,self.f.qp+1):
            zt,rt=self.f.top(Q)
            z=zt/3;r=self.f.radius(Q,z)
            self.assertAlmostEqual(1/r**2-1/rt**2,
                                   2*b_Q(Q,1.,7.)*(z-zt),places=13)
            if zt<1:
                with self.assertRaises(ValueError):self.f.radius(Q,zt+.01)

    def test_strip_volumes_match_initial_annular_mass(self):
        intervals=[(self.f.qm/3,2*self.f.qm/3),
                   ((self.f.qm+self.f.Qt)/2,(self.f.qm+3*self.f.Qt)/4),
                   ((3*self.f.Qt+self.f.qp)/4,(self.f.Qt+self.f.qp)/2),
                   (self.f.qp+1,self.f.qp+2)]
        for qa,qb in intervals:
            poly=self.f.strip(qa,qb,samples=2001)
            # Exact integral over each straight segment in the sampled (r,z) polygon.
            mass=abs(math.pi*math.fsum((a[0]**2+a[0]*b[0]+b[0]**2)/3*(b[1]-a[1])
                     for a,b in zip(poly,poly[1:]+poly[:1])))
            target=2*math.pi*(math.sqrt(qb)-math.sqrt(qa))
            self.assertLess(abs(mass-target)/target,5e-6)

    def test_grid_has_no_gaps_and_does_not_smooth_the_front(self):
        import numpy as np
        r,z,Q=self.f.grid(nr=301,nz=101,nq=1000)
        self.assertTrue(np.isfinite(Q).all())
        self.assertGreaterEqual(float(Q.min()),0.)
        j=95;rf=self.f.front_r(float(z[j]));i=int(np.searchsorted(r,rf))
        self.assertLess(Q[j,i-1],self.f.q_m(float(z[j]))+.001)
        self.assertGreater(Q[j,i],self.f.q_p(float(z[j]))-.001)


if __name__=='__main__':unittest.main()
