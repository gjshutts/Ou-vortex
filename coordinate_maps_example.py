"""Run this file to explore the coordinate maps; no installations are needed."""

import math
from coordinate_maps import (
    R_from_r, r_from_R, amplification, Zhat_from_z, z_from_Zhat,
    density, vertical_mass_weight, mass_metric, DRY_COEFFICIENTS,
)

# Assign the physical parameters to names. Change these to experiment.
r0 = 0.095
rmax = 7.0
H0 = 0.75
rho_s = 1.0

print("RADIAL MAP")
print("       r             R       recovered r")
for r in [r0, 0.1, 1.0, 5.0, rmax]:
    R = R_from_r(r, r0)
    recovered = r_from_R(R, r0)
    # An f-string inserts values into text. .6f means six decimal places.
    print(f"{r:8.3f}    {R:12.6f}    {recovered:12.6f}")

print("\nMOIST HEIGHT MAP AND DENSITY")
print("    z       A(z)       Zhat    recovered z      rho0       W(Zhat)")
for z in [0.0, 0.25, 0.5, 0.75, 1.0]:
    Zhat = Zhat_from_z(z)
    print(f"{z:5.2f}   {amplification(z):9.6f}   {Zhat:9.6f}"
          f"   {z_from_Zhat(Zhat):11.6f}   {density(z, H0, rho_s):9.6f}"
          f"   {vertical_mass_weight(Zhat, H0=H0, rho_s=rho_s):9.6f}")

# The same geometry functions cover all three formulations.
z = 0.5
R = R_from_r(1.0, r0)
print("\nTHREE MASS MEASURES at r=1, physical z=0.5")
print("These are densities per unit transformed-coordinate area, not cell masses.")
print("Dry, constant density:",
      mass_metric(R, z, r0, coefficients=DRY_COEFFICIENTS, H0=math.inf))
print("Dry, density decay:  ",
      mass_metric(R, z, r0, coefficients=DRY_COEFFICIENTS, H0=H0))
print("Moist, density decay:",
      mass_metric(R, Zhat_from_z(z), r0, H0=H0))

bm = 16.0
print("\nPHYSICAL BUOYANCY for a parcel with conserved bm=16")
for z in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    b = amplification(z) * bm
    print(f"z={z:.1f}: b={b:.6f}")
