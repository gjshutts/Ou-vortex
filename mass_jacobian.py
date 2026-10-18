"""Edge-based derivatives J[i][j] = partial mass_i / partial S_j.

Q and bm are held fixed. Only the plane offsets S vary.
The matrix is returned as a list of rows, without NumPy dependencies.
"""

import math
from coordinate_maps import DEFAULT_COEFFICIENTS, mass_metric
from mass_integration import integrate_1d


def edge_mass_integral(start, end, r0=0.095, coefficients=DEFAULT_COEFFICIENTS,
                       H0=0.75, rho_s=1.0, absolute_tolerance=1e-10,
                       relative_tolerance=1e-9, use_exact_uniform=True):
    """Integrate the physical mass metric along a straight transformed edge.

    This is integral J_mass ds, with ds measured in (R,Zhat). It is not
    the signed Green's-theorem integral used to calculate a cell's mass.
    Reversing this edge does not change the result. Horizontal edges can
    contribute here, although they contribute zero to that cell-mass integral.
    """
    if len(start) != 2 or len(end) != 2:
        raise ValueError("Edge endpoints must be pairs (R, Zhat).")
    mass_metric(*start, r0=r0, coefficients=coefficients, H0=H0, rho_s=rho_s)
    mass_metric(*end, r0=r0, coefficients=coefficients, H0=H0, rho_s=rho_s)
    length = math.dist(start, end)
    if length == 0:
        return 0.0
    if use_exact_uniform and H0 == math.inf and coefficients[1] == 0 and coefficients[2] == 0:
        Rinf = 1/(2*r0**2)
        return rho_s*length*math.pi/(2*(Rinf-start[0])*(Rinf-end[0]))

    def integrand(t):
        if t == 0:
            R, Zhat = start
        elif t == 1:
            R, Zhat = end
        else:
            R = start[0] + t * (end[0] - start[0])
            Zhat = start[1] + t * (end[1] - start[1])
        return length * mass_metric(R, Zhat, r0, coefficients, H0, rho_s)

    return integrate_1d(integrand, 0, 1, absolute_tolerance,
                        relative_tolerance=relative_tolerance)


def mass_jacobian(planes, polygons, r0=0.095, coefficients=DEFAULT_COEFFICIENTS,
                  H0=0.75, rho_s=1.0, edge_tolerance=1e-9,
                  quadrature_tolerance=1e-12, quadrature_relative_tolerance=1e-9):
    """Assemble the symmetric Jacobian from shared material interfaces.

    For a shared edge between cells i and j, let
    D = sqrt((Qi-Qj)**2 + (bmi-bmj)**2), K = integral(J_mass ds)/D.
    Add K to both diagonal entries and -K to both off-diagonal entries.
    Outer walls do not move with S and make no contribution.

    Supply the CURRENT polygons produced by cell_polygons for these planes.
    edge_tolerance is distance to a plane-pair line in transformed coordinates.
    quadrature_tolerance estimates error in each edge's K coefficient.
    """
    if not math.isfinite(quadrature_relative_tolerance) or quadrature_relative_tolerance < 0:
        raise ValueError("Relative quadrature tolerance must be finite and nonnegative.")
    n = len(planes)
    if len(polygons) != n:
        raise ValueError("Supply one polygon per plane, including empty cells.")
    if any(len(p) != 3 or not all(math.isfinite(v) for v in p) for p in planes):
        raise ValueError("Planes must contain finite (Q, bm, S) values.")
    for value in (edge_tolerance, quadrature_tolerance):
        if not math.isfinite(value) or value <= 0:
            raise ValueError("Tolerances must be finite and positive.")
    # Validate the physical inputs even if no shared interfaces are present.
    mass_metric(0, 0, r0, coefficients, H0, rho_s)
    if len(set(tuple(p) for p in planes)) != n:
        raise ValueError("Identical planes have ambiguous offset derivatives; remove duplicates.")

    matrix = [[0.0 for j in range(n)] for i in range(n)]
    for i, polygon in enumerate(polygons):
        for k, start in enumerate(polygon):
            end = polygon[(k + 1) % len(polygon)]
            for j in range(i + 1, n):  # Count each cell pair only once.
                if not polygons[j]:
                    continue
                a = planes[i][0] - planes[j][0]
                b = planes[i][1] - planes[j][1]
                c = planes[i][2] - planes[j][2]
                normal = math.hypot(a, b)
                if normal == 0:
                    continue  # Parallel planes with different offsets cannot share an edge.

                def line_distance(point):
                    return abs(math.fsum((a*point[0], b*point[1], c))) / normal

                # Equality at both endpoints implies equality along the edge.
                if max(line_distance(start), line_distance(end)) > edge_tolerance:
                    continue
                integral = edge_mass_integral(start, end, r0, coefficients, H0, rho_s,
                                              absolute_tolerance=quadrature_tolerance * normal,
                                              relative_tolerance=quadrature_relative_tolerance)
                K = integral / normal
                matrix[i][j] -= K
                matrix[j][i] -= K

    # Set each diagonal to the negative sum of its off-diagonal entries.
    # Each row then sums to zero, up to floating-point rounding.
    for i in range(n):
        matrix[i][i] = -math.fsum(matrix[i][j] for j in range(n) if j != i)
    return matrix


def finite_difference_jacobian(planes, domain, step=1e-5, **mass_parameters):
    """An independent numerical check: perturb one S at a time and rebuild.

    This is for verification, not the planned Newton solver. It is slower
    because it reconstructs cells and recalculates masses for every column.
    Both perturbed configurations must stay in a differentiable topology.
    """
    from cell_geometry import cell_polygons
    from mass_integration import cell_masses
    if not math.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive.")
    n = len(planes)
    matrix = [[0.0 for j in range(n)] for i in range(n)]
    for j in range(n):
        plus = list(planes)
        minus = list(planes)
        Q, bm, S = planes[j]
        plus[j] = (Q, bm, S + step)
        minus[j] = (Q, bm, S - step)
        masses_plus = cell_masses(cell_polygons(plus, domain), **mass_parameters)
        masses_minus = cell_masses(cell_polygons(minus, domain), **mass_parameters)
        for i in range(n):
            matrix[i][j] = (masses_plus[i] - masses_minus[i]) / (2 * step)
    return matrix
