"""Convex material cells from the upper envelope of planes.

A point is a pair (R, Zhat). A polygon is a list of its boundary vertices,
without repeating the first vertex at the end. A plane is a tuple (Q, bm, S):
P(R, Zhat) = Q*R + bm*Zhat + S. Cell i is where P_i >= all other planes.

This module uses only Python's standard library. It calculates geometry,
not mass, and does not solve for the offsets S.
"""

import math


def rectangle(R_left, R_right, Z_bottom=0.0, Z_top=1.0):
    """Four vertices, ordered counterclockwise in transformed coordinates."""
    if not all(math.isfinite(v) for v in (R_left, R_right, Z_bottom, Z_top)):
        raise ValueError("Rectangle bounds must be finite.")
    if R_left >= R_right or Z_bottom >= Z_top:
        raise ValueError("Rectangle bounds must have positive width and height.")
    return [(R_left, Z_bottom), (R_right, Z_bottom),
            (R_right, Z_top), (R_left, Z_top)]


def signed_area(polygon):
    """Positive for counterclockwise vertices; zero for a line or empty cell.

    Triangles relative to the first vertex reduce cancellation for small
    cells at large R. This is area in (R,Zhat), not physical volume or mass.
    """
    if len(polygon) < 3:
        return 0.0
    origin = polygon[0]
    terms = []
    for i in range(1, len(polygon) - 1):
        u = (polygon[i][0] - origin[0], polygon[i][1] - origin[1])
        v = (polygon[i + 1][0] - origin[0], polygon[i + 1][1] - origin[1])
        terms.append(u[0] * v[1] - u[1] * v[0])
    return math.fsum(terms) / 2


def _clean_polygon(polygon, tolerance):
    """Remove repeated adjacent vertices and discard zero-area pieces."""
    out = []
    for point in polygon:
        point = tuple(point)
        if not out or math.dist(point, out[-1]) > tolerance:
            out.append(point)
    if len(out) > 1 and math.dist(out[0], out[-1]) <= tolerance:
        out.pop()
    area = signed_area(out)
    if area == 0:
        return []
    if area < 0:
        out.reverse()
    return out


def clip_half_plane(polygon, a, b, c, tolerance=1e-12):
    """Keep the portion satisfying a*R + b*Zhat + c >= 0.

    Each edge has four possibilities: both endpoints inside, leaving,
    entering, or both outside. A crossing is found by linear interpolation.

    The tolerance is a distance in the nondimensional coordinate plane.
    Normalising (a,b) makes clipping insensitive to a common positive scale
    factor in the inequality. Points within tolerance of the line are treated
    as on it. Features smaller than tolerance cannot be resolved reliably.
    """
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be finite and nonnegative.")
    if not all(math.isfinite(v) for v in (a, b, c)):
        raise ValueError("Half-plane coefficients must be finite.")
    if any(len(p) != 2 or not all(math.isfinite(v) for v in p) for p in polygon):
        raise ValueError("Vertices must be finite pairs (R, Zhat).")
    if not polygon:
        return []

    normal = math.hypot(a, b)
    if normal == 0:
        # A constant inequality either keeps everything or rejects everything.
        return _clean_polygon(polygon, tolerance) if c >= 0 else []
    a, b, c = a / normal, b / normal, c / normal

    def distance(point):
        value = math.fsum((a * point[0], b * point[1], c))
        return 0.0 if abs(value) <= tolerance else value

    out = []
    for i, p1 in enumerate(polygon):
        p2 = polygon[(i + 1) % len(polygon)]  # Wrap the last edge to vertex 0.
        d1, d2 = distance(p1), distance(p2)
        inside1, inside2 = d1 >= 0, d2 >= 0

        if inside1 and inside2:
            out.append(p2)
        elif inside1 != inside2:
            # The endpoints are on different sides, so 0 <= t <= 1.
            t = d1 / (d1 - d2)
            crossing = (p1[0] + t * (p2[0] - p1[0]),
                        p1[1] + t * (p2[1] - p1[1]))
            out.append(crossing)
            if inside2:
                out.append(p2)
        # When both endpoints are outside, append nothing.

    return _clean_polygon(out, tolerance)


def plane_value(plane, point):
    """Evaluate P = Q*R + bm*Zhat + S."""
    Q, bm, S = plane
    R, Zhat = point
    return Q * R + bm * Zhat + S


def cell_polygon(index, planes, domain, tolerance=1e-12):
    """Construct one cell by repeatedly clipping a convex domain.

    Python indices start at 0: the first plane is planes[0].
    Exactly identical planes are assigned to the lowest index so they do not
    duplicate a cell. Different material labels normally give distinct planes.
    """
    if not isinstance(index, int) or not 0 <= index < len(planes):
        raise ValueError("index must identify a plane, starting at zero.")
    if any(len(p) != 3 or not all(math.isfinite(v) for v in p) for p in planes):
        raise ValueError("Each plane must be a finite tuple (Q, bm, S).")

    poly = list(domain)
    Qi, bmi, Si = planes[index]
    for j, (Qj, bmj, Sj) in enumerate(planes):
        if j == index:
            continue
        a, b, c = Qi - Qj, bmi - bmj, Si - Sj
        if a == 0 and b == 0 and c == 0 and j < index:
            return []
        poly = clip_half_plane(poly, a, b, c, tolerance)
        if not poly:
            break
    # Validate and orient even when there is only one plane.
    return clip_half_plane(poly, 0.0, 0.0, 1.0, tolerance)


def cell_polygons(planes, domain, tolerance=1e-12):
    """Return one polygon per plane, retaining [] for any empty cell."""
    return [cell_polygon(i, planes, domain, tolerance) for i in range(len(planes))]
