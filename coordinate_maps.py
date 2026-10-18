"""Coordinate maps for the rotating zero-PV model.

This first building block uses only Python's standard library.
All functions take one number at a time (not NumPy arrays).

Physical height z and moist height Zhat both range from 0 to 1.
The unit physical column represents 16 km in our atmospheric example.
The inner radius r0 is positive; r = 0 is outside this coordinate system.
"""

import math

# A tuple is an ordered collection of values: here c0, c1 and c2.
# These are the unnormalised polynomial coefficients from Mathematica.
DEFAULT_COEFFICIENTS = (1.0, 0.548584, 4.273361)
DRY_COEFFICIENTS = (1.0, 0.0, 0.0)


def _positive_number(value, name):
    """Reject zero, negative or non-finite values with a useful message."""
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number.")


def _unit_height(value, name):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must lie between 0 and 1, inclusive.")


def _profile(coefficients):
    """Validate the quadratic, then return its coefficients and integral.

    Names starting with an underscore are helpers for this module.
    Call the public functions below in your own calculations.
    """
    if len(coefficients) != 3:
        raise ValueError("Supply three coefficients: (c0, c1, c2).")
    c0, c1, c2 = coefficients  # Unpack the tuple into three names.
    if not all(math.isfinite(c) for c in coefficients):
        raise ValueError("The amplification coefficients must be finite.")

    # A positive A is essential: otherwise the height map is not invertible.
    values = [c0, c0 + c1 + c2]
    if c2 > 0:
        vertex = -c1 / (2 * c2)
        if 0 < vertex < 1:
            values.append(c0 + c1 * vertex + c2 * vertex**2)
    if min(values) <= 0:
        raise ValueError("The amplification polynomial must be positive on [0, 1].")

    normalisation = c0 + c1 / 2 + c2 / 3
    return c0, c1, c2, normalisation


# 1. Radial coordinate and its inverse ----------------------------------------

def R_from_r(r, r0=0.095):
    """R = (1/r0**2 - 1/r**2)/2, for r >= r0."""
    _positive_number(r0, "r0")
    _positive_number(r, "r")
    if r < r0:
        raise ValueError("r must be at or outside the inner wall r0.")
    return 0.5 * (1 / r0**2 - 1 / r**2)


def r_from_R(R, r0=0.095):
    """Invert the radial map; R must be below Rinf = 1/(2*r0**2)."""
    _positive_number(r0, "r0")
    Rinf = 1 / (2 * r0**2)
    if not math.isfinite(R) or not 0 <= R < Rinf:
        raise ValueError("R must satisfy 0 <= R < 1/(2*r0**2).")
    return 1 / math.sqrt(1 / r0**2 - 2 * R)


# 2. Amplification, moist coordinate and inverse -------------------------------

def amplification(z, coefficients=DEFAULT_COEFFICIENTS):
    """A(z), normalised so its integral from 0 to 1 is one."""
    _unit_height(z, "z")
    c0, c1, c2, normalisation = _profile(coefficients)
    return (c0 + c1 * z + c2 * z**2) / normalisation


def Zhat_from_z(z, coefficients=DEFAULT_COEFFICIENTS):
    """Integrate A exactly: Zhat = (c0*z + c1*z**2/2 + c2*z**3/3)/norm."""
    _unit_height(z, "z")
    c0, c1, c2, normalisation = _profile(coefficients)
    return (c0 * z + c1 * z**2 / 2 + c2 * z**3 / 3) / normalisation


def z_from_Zhat(Zhat, coefficients=DEFAULT_COEFFICIENTS):
    """Invert the monotone cubic by bisection.

    Bisection repeatedly halves a bracket containing the answer. It is slower
    than Newton's method, but particularly simple and reliable for this first
    implementation. Sixty halvings reach floating-point precision on [0, 1].
    """
    _unit_height(Zhat, "Zhat")
    c0, c1, c2, normalisation = _profile(coefficients)
    if Zhat == 0 or Zhat == 1:
        return float(Zhat)

    lower = 0.0
    upper = 1.0
    for _ in range(60):
        midpoint = (lower + upper) / 2
        mapped = (c0 * midpoint + c1 * midpoint**2 / 2
                  + c2 * midpoint**3 / 3) / normalisation
        if mapped < Zhat:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


# 3. Reference density and mass measure ---------------------------------------

def density(z, H0=0.75, rho_s=1.0):
    """rho0(z) = rho_s*exp(-z/H0); H0=math.inf gives constant density."""
    _unit_height(z, "z")
    _positive_number(rho_s, "rho_s")
    if H0 == math.inf:
        return rho_s
    _positive_number(H0, "H0")
    return rho_s * math.exp(-z / H0)


def vertical_mass_weight(Zhat, coefficients=DEFAULT_COEFFICIENTS,
                         H0=0.75, rho_s=1.0):
    """W(Zhat) = rho0(z(Zhat))/A(z(Zhat)), including dz/dZhat = 1/A."""
    z = z_from_Zhat(Zhat, coefficients)
    return density(z, H0, rho_s) / amplification(z, coefficients)


def mass_metric(R, Zhat, r0=0.095, coefficients=DEFAULT_COEFFICIENTS,
                H0=0.75, rho_s=1.0):
    """Return J such that dm = J(R,Zhat) dR dZhat.

    Axisymmetric volume is 2*pi*r*dr*dz. Since dr/dR=r**3 and
    dz/dZhat=1/A, J = 2*pi*r**4*rho0/A.
    This equals pi*W/[2*(Rinf-R)**2] in the Mathematica formulation.
    """
    r = r_from_R(R, r0)
    W = vertical_mass_weight(Zhat, coefficients, H0, rho_s)
    return 2 * math.pi * r**4 * W
