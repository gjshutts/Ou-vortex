"""Inspect the 3-by-3 Jacobian for our simple three-cell partition."""
import math
from cell_geometry import rectangle, cell_polygons
from mass_jacobian import mass_jacobian, finite_difference_jacobian

# The third offset is -0.6 so no interface passes through a domain corner.
planes = [(0.0, 0.0, 0.0), (1.0, 0.0, -0.5), (0.0, 1.0, -0.6)]
domain = rectangle(0.0, 1.0)
cells = cell_polygons(planes, domain)
J = mass_jacobian(planes, cells)
J_fd = finite_difference_jacobian(planes, domain, absolute_tolerance=1e-14)

print("MASS JACOBIAN: rows = cell masses; columns = plane offsets")
for row in J:
    print(" ".join(f"{value: .10e}" for value in row))
print("\nRow sums:", [math.fsum(row) for row in J])
print("Maximum finite-difference difference:",
      max(abs(J[i][j]-J_fd[i][j]) for i in range(3) for j in range(3)))

# Raising only S0 gives each cell's first-order mass change as column 0.
change_in_S = [0.001, 0.0, 0.0]
change_in_mass = [math.fsum(J[i][j]*change_in_S[j] for j in range(3))
                  for i in range(3)]
print("\nFor offset changes", change_in_S)
print("Predicted first-order mass changes:", change_in_mass)
print("Their sum:", math.fsum(change_in_mass))
print("Increasing S0 expands cell 0 at the expense of its neighbours.")
print("Adding the SAME offset to all planes changes no cells and no masses.")
