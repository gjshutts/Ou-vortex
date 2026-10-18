"""Physical mass of polygon cells in (R, Zhat), using Green's theorem.

This is building block 3, not an equilibrium solver. It uses only the Python
standard library and the two preceding modules. Density is prescribed;
physical height and transformed height both span the unit column.
"""

import math
from coordinate_maps import DEFAULT_COEFFICIENTS, vertical_mass_weight, r_from_R, density
from cell_geometry import signed_area


def integrate_1d(function, lower, upper, absolute_tolerance=1e-10, max_depth=24,
                 relative_tolerance=0.0):
    """Adaptive Simpson quadrature, with an estimated absolute error target.

    Compare one Simpson panel with two smaller panels. If they disagree,
    subdivide again and split the error allowance between the two halves.
    The tolerance is an estimate, not a rigorous mathematical error bound.
    An optional relative tolerance supports integrals with large magnitudes.
    Failure to reach the error criterion raises an exception.
    """
    if not math.isfinite(absolute_tolerance) or absolute_tolerance <= 0:
        raise ValueError("absolute_tolerance must be finite and positive.")
    if not math.isfinite(relative_tolerance) or relative_tolerance < 0:
        raise ValueError("relative_tolerance must be finite and nonnegative.")
    if not isinstance(max_depth, int) or max_depth < 0:
        raise ValueError("max_depth must be a nonnegative integer.")
    if not all(math.isfinite(v) for v in (lower, upper)):
        raise ValueError("Integration bounds must be finite.")
    if lower == upper:
        return 0.0

    def evaluate(x):
        result = function(x)
        if not math.isfinite(result):
            raise ValueError("The integrand returned a non-finite value.")
        return result

    def panel(a, b, fa, fm, fb):
        return (b - a) * (fa + 4 * fm + fb) / 6

    def refine(a, b, fa, fm, fb, coarse, tolerance, depth):
        midpoint = (a + b) / 2
        f_left = evaluate((a + midpoint) / 2)
        f_right = evaluate((midpoint + b) / 2)
        left = panel(a, midpoint, fa, f_left, fm)
        right = panel(midpoint, b, fm, f_right, fb)
        difference = left + right - coarse
        if abs(difference) <= 15 * (tolerance + relative_tolerance * abs(left + right)):
            # The correction removes the leading Simpson error term.
            return left + right + difference / 15
        if depth == 0:
            raise ArithmeticError("Quadrature did not reach its tolerance; increase max_depth.")
        return math.fsum((
            refine(a, midpoint, fa, f_left, fm, left, tolerance / 2, depth - 1),
            refine(midpoint, b, fm, f_right, fb, right, tolerance / 2, depth - 1),
        ))

    fa = evaluate(lower)
    fm = evaluate((lower + upper) / 2)
    fb = evaluate(upper)
    return refine(lower, upper, fa, fm, fb, panel(lower, upper, fa, fm, fb),
                  absolute_tolerance, max_depth)


def polygon_mass(polygon, r0=0.095, coefficients=DEFAULT_COEFFICIENTS,
                 H0=0.75, rho_s=1.0, absolute_tolerance=1e-10,
                 use_exact_uniform=True):
    """Integrate physical mass over one transformed polygon.

    J(R,Zhat)=pi*W(Zhat)/[2*(Rinf-R)**2]. A radial primitive is
    F(R)=pi/2 * R/[Rinf*(Rinf-R)], chosen to vanish at R=0.
    Thus mass = integral around the CCW boundary of F(R)*W(Zhat) dZhat.
    Horizontal edges contribute zero. Clockwise input is reversed locally.
    An empty or zero-area cell has mass zero.
    """
    if not math.isfinite(absolute_tolerance) or absolute_tolerance <= 0:
        raise ValueError("absolute_tolerance must be finite and positive.")
    # Validate the physical parameters even for an empty cell.
    r_from_R(0, r0)
    density(0, H0, rho_s)
    vertical_mass_weight(0, coefficients, H0, rho_s)
    for point in polygon:
        if len(point) != 2:
            raise ValueError("Each vertex must be a pair (R, Zhat).")
        R, Zhat = point
        r_from_R(R, r0)
        if not math.isfinite(Zhat) or not 0 <= Zhat <= 1:
            raise ValueError("Vertex Zhat must lie in [0, 1].")

    area = signed_area(polygon)
    if area == 0:
        return 0.0
    vertices = list(polygon) if area > 0 else list(reversed(polygon))
    Rinf = 1 / (2 * r0**2)
    contributions = []
    if use_exact_uniform and H0 == math.inf and coefficients[1] == 0 and coefficients[2] == 0:
        # A constant positive polynomial normalises to A=1. Integrate exactly.
        for i,start in enumerate(vertices):
            end = vertices[(i+1)%len(vertices)]
            dR,dZ = end[0]-start[0],end[1]-start[1]
            if dZ == 0:
                continue
            u = Rinf-start[0]
            t = dR/u
            if abs(t) < 1e-3:
                # [-log(1-t)/t] - 1 = t/2 + t**2/3 + ...
                q_minus_one = math.fsum(t**k/(k+1) for k in range(1,8))
            else:
                q_minus_one = -math.log1p(-t)/t-1
            mean_primitive = math.pi/2*(start[0]/(Rinf*u)+q_minus_one/u)
            contributions.append(rho_s*dZ*mean_primitive)
        mass = math.fsum(contributions)
        if mass < 0:
            raise ArithmeticError("Negative exact mass: check polygon geometry.")
        return mass
    edge_tolerance = absolute_tolerance / len(vertices)

    # Reuse W at repeated height samples, especially shared vertices.
    weights = {}

    def weight(Zhat):
        if Zhat not in weights:
            weights[Zhat] = vertical_mass_weight(Zhat, coefficients, H0, rho_s)
        return weights[Zhat]

    for i, start in enumerate(vertices):
        end = vertices[(i + 1) % len(vertices)]
        dR = end[0] - start[0]
        dZ = end[1] - start[1]
        if dZ == 0:
            continue

        def integrand(t):
            # Interpolate between endpoints; use them exactly at t=0 and 1.
            R = start[0] if t == 0 else end[0] if t == 1 else start[0] + t * dR
            Zhat = start[1] if t == 0 else end[1] if t == 1 else start[1] + t * dZ
            primitive = (math.pi / 2) * R / (Rinf * (Rinf - R))
            return primitive * weight(Zhat) * dZ

        contributions.append(integrate_1d(integrand, 0, 1, edge_tolerance))

    mass = math.fsum(contributions)
    if mass < 0:
        raise ArithmeticError("Negative mass: check polygon simplicity or quadrature tolerance.")
    return mass


def cell_masses(polygons, **parameters):
    """One physical mass per cell; **parameters forwards named arguments."""
    return [polygon_mass(poly, **parameters) for poly in polygons]


def annular_column_mass(r_inner, r_outer, z_bottom=0.0, z_top=1.0,
                        H0=0.75, rho_s=1.0):
    """Exact benchmark: a physical annulus spanning a specified height range."""
    if not all(math.isfinite(r) for r in (r_inner, r_outer)) or not 0 <= r_inner < r_outer:
        raise ValueError("Require finite radii with 0 <= r_inner < r_outer.")
    if not 0 <= z_bottom < z_top <= 1:
        raise ValueError("Require 0 <= z_bottom < z_top <= 1.")
    density(z_bottom, H0, rho_s)
    depth = z_top - z_bottom
    column = rho_s * depth if H0 == math.inf else (
        rho_s * H0 * math.exp(-z_bottom / H0) * (-math.expm1(-depth / H0)))
    return math.pi * (r_outer**2 - r_inner**2) * column
