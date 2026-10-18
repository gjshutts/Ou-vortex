"""Initial rest-state labels, exact annular target masses and an offset seed."""
import math
from coordinate_maps import R_from_r
from mass_integration import annular_column_mass


def generate_initial_state(n, r_inner, r_outer, bm0, alpha, H0=0.75, rho_s=1.0):
    """Uniform-in-r initial annuli; Q=r_center**4/4, bm=bm0*exp(-alpha*Q).

    The offset seed is the dry-derived construction used in our Mathematica
    model. It is only a numerical starting guess, not a moist equilibrium.
    Each annulus initially spans the entire physical unit-depth column.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer.")
    if not all(math.isfinite(v) for v in (r_inner, r_outer, bm0, alpha)):
        raise ValueError("Initial-state parameters must be finite.")
    if not 0 < r_inner < r_outer or bm0 < 0 or alpha < 0:
        raise ValueError("Require 0 < r_inner < r_outer and nonnegative bm0, alpha.")
    edges = [r_inner + (r_outer-r_inner)*i/n for i in range(n+1)]
    centers = [(a+b)/2 for a,b in zip(edges[:-1], edges[1:])]
    Q_edges = [r**4/4 for r in edges]
    Q = [r**4/4 for r in centers]
    bm = [bm0*math.exp(-alpha*q) for q in Q]
    targets = [annular_column_mass(a, b, H0=H0, rho_s=rho_s)
               for a,b in zip(edges[:-1], edges[1:])]

    def x_coth_x(x):
        # The limit at zero is 1. The series avoids dividing by zero.
        if abs(x) < 1e-4:
            return 1+x**2/3-x**4/45+2*x**6/945
        return x/math.tanh(x)

    Rinf = 1/(2*r_inner**2)
    Rm = []
    for r, q in zip(edges[:-1], Q_edges[:-1]):
        derivative = -alpha*bm0*math.exp(-alpha*q)
        Rm.append(Rinf-x_coth_x(derivative*r**2)/(2*r**2))
    Rm.append(R_from_r(r_outer, r_inner))
    P_edges = [0.0]
    for i in range(n):
        P_edges.append(P_edges[-1]+(Q_edges[i]+Q_edges[i+1])/2*(Rm[i+1]-Rm[i]))
    offsets = [(P_edges[i]+P_edges[i+1])/2-Q[i]*(Rm[i]+Rm[i+1])/2-bm[i]/2
               for i in range(n)]
    return {"r_edges": edges, "Q": Q, "bm": bm, "target_masses": targets,
            "initial_offsets": offsets}
