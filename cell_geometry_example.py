"""Three deliberately simple planes with an exactly known partition."""

from pathlib import Path
from cell_geometry import rectangle, clip_half_plane, cell_polygons, signed_area

# Use a unit square for this geometry exercise, not an atmospheric equilibrium.
domain = rectangle(0.0, 1.0)
planes = [
    (0.0, 0.0, 0.0),   # P0 = 0
    (1.0, 0.0, -0.5),  # P1 = R - 0.5
    (0.0, 1.0, -0.5),  # P2 = Zhat - 0.5
]

# Watch cell 0 being constructed. P0 >= P1 means R <= 0.5.
first_cut = clip_half_plane(domain, -1.0, 0.0, 0.5)
# P0 >= P2 also requires Zhat <= 0.5.
second_cut = clip_half_plane(first_cut, 0.0, -1.0, 0.5)
print("Starting rectangle:", domain)
print("After keeping R <= 0.5:", first_cut)
print("After also keeping Zhat <= 0.5:", second_cut)

cells = cell_polygons(planes, domain)
print("\nUPPER-ENVELOPE CELLS")
for i, poly in enumerate(cells):
    print(f"Cell {i}: vertices={poly}, area={signed_area(poly):.6f}")
print("Total coordinate area:", sum(signed_area(poly) for poly in cells))
print("Expected areas: 0.25, 0.375, 0.375. These are NOT physical masses.")

# Write a vector drawing using the computed vertices: no plotting packages.
colours = ["#e6edf5", "#d9eaf7", "#f5d7d9"]
parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="700" height="650" viewBox="0 0 700 650">',
         '<rect width="700" height="650" fill="white"/>',
         '<g font-family="Arial,sans-serif" fill="#13283e">',
         '<text x="350" y="35" text-anchor="middle" font-size="23">Cell geometry: three competing planes</text>']
for i, poly in enumerate(cells):
    points = " ".join(f"{90 + 480*R},{560 - 480*Z}" for R, Z in poly)
    parts.append(f'<polygon points="{points}" fill="{colours[i]}" stroke="#13283e" stroke-width="2"/>')
for R, Z, label in [(0.25, 0.25, "P0 = 0"), (0.78, 0.32, "P1 = R - 0.5"),
                     (0.32, 0.78, "P2 = Zhat - 0.5")]:
    parts.append(f'<text x="{90 + 480*R}" y="{560 - 480*Z}" text-anchor="middle" font-size="18">{label}</text>')
parts.extend(['<text x="330" y="605" text-anchor="middle" font-size="20">R</text>',
              '<text x="30" y="325" text-anchor="middle" font-size="20">Zhat</text>',
              '<text x="90" y="585" text-anchor="middle" font-size="16">0</text>',
              '<text x="570" y="585" text-anchor="middle" font-size="16">1</text>',
              '<text x="70" y="90" text-anchor="middle" font-size="16">1</text>',
              '<text x="350" y="637" text-anchor="middle" font-size="15">Geometry example in transformed coordinates; no equilibrium solve.</text>',
              '</g></svg>'])
output = Path(__file__).with_name("cell_geometry_example.svg")
output.write_text("\n".join(parts))
print("Vector drawing saved to:", output)
