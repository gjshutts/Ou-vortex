"""Compare coordinate areas with physical masses, and check a known annulus."""

import math
from coordinate_maps import R_from_r, Zhat_from_z, DRY_COEFFICIENTS
from cell_geometry import rectangle, cell_polygons, signed_area
from mass_integration import cell_masses, polygon_mass, annular_column_mass

r0 = 0.095
# A simple geometry exercise with three planes in a unit (R,Zhat) square.
domain = rectangle(0.0, 1.0)
planes = [(0.0, 0.0, 0.0), (1.0, 0.0, -0.5), (0.0, 1.0, -0.5)]
cells = cell_polygons(planes, domain)

cases = [
    ("Dry, constant density", dict(coefficients=DRY_COEFFICIENTS, H0=math.inf)),
    ("Dry, density decay", dict(coefficients=DRY_COEFFICIENTS, H0=0.75)),
    ("Moist, density decay", dict(H0=0.75)),
]
print("SAME TRANSFORMED POLYGONS, THREE PHYSICAL MASS MEASURES")
print("The moist case maps these polygons to different physical heights.")
for label, parameters in cases:
    masses = cell_masses(cells, r0=r0, **parameters)
    print(f"\n{label}")
    for i, (poly, mass) in enumerate(zip(cells, masses)):
        print(f"Cell {i}: coordinate area={signed_area(poly):.6f}, mass={mass:.10f}")
    exact = annular_column_mass(r0, 1 / math.sqrt(1 / r0**2 - 2),
                               H0=parameters["H0"])
    print(f"Sum of cell masses = {math.fsum(masses):.12f}")
    print(f"Exact domain mass  = {exact:.12f}")

print("\nA PHYSICAL ANNULUS: r=0.5 to 2.0, z=0.2 to 0.8")
# The same physical rectangle has different vertical bounds in the moist map.
for label, parameters in cases:
    coefficients = parameters.get("coefficients", None)
    if coefficients is None:
        lower, upper = Zhat_from_z(0.2), Zhat_from_z(0.8)
    else:
        lower = Zhat_from_z(0.2, coefficients)
        upper = Zhat_from_z(0.8, coefficients)
    poly = rectangle(R_from_r(0.5, r0), R_from_r(2.0, r0), lower, upper)
    computed = polygon_mass(poly, r0=r0, **parameters)
    exact = annular_column_mass(0.5, 2.0, 0.2, 0.8, H0=parameters["H0"])
    print(f"{label}: computed={computed:.10f}, exact={exact:.10f}, error={computed-exact:.2e}")
