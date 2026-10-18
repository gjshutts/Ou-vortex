"""Dry rotational Ou front: tip launch, front ODEs, and outer shooting.

Reference: MATLAB/solve_zeroPV_vortex_front version used in paper.m.
Uses physical z, with lids at 0 and 1, and R = -1/(2*r**2).
No finite outer wall is imposed. All public inputs are scalar numbers.

The default integrator is an adaptive Dormand--Prince RK45 method, written
using the standard library. Optional method='BDF' uses scipy.integrate;
it is the stiff alternative to MATLAB's ode15s, not the same algorithm.
"""

from math import exp, expm1, fsum, isfinite, log1p, sqrt


def _parameters(alpha, b0):
    if not all(isfinite(v) and v > 0 for v in (alpha, b0)):
        raise ValueError("alpha and b0 must be finite and positive")


def b_Q(Q, alpha, b0):
    """Derivative of b(Q) = b0 exp(-alpha Q)."""
    return -alpha * b0 * exp(-alpha * Q)


def jump_slope(Qm, Qp, alpha, b0):
    """Margules slope dz/dR = -Delta Q/Delta b, positive here."""
    if not 0 < Qm < Qp:
        raise ValueError("Require 0 < Qm < Qp")
    return -(Qp - Qm) / (b0 * exp(-alpha * Qm) * expm1(-alpha * (Qp - Qm)))


def psi_of_Q(Q, alpha, b0):
    return -sqrt(2 * b0 * (0.5 + alpha * Q) * exp(-alpha * Q)
                 / (alpha * Q**1.5))


def tip_residual(Q, alpha, b0):
    if Q <= 0:
        return float('nan')
    psi = psi_of_Q(Q, alpha, b0)
    if not -1 < psi < 0:
        return float('nan')
    return sqrt(Q) * exp(alpha * Q) * (psi - log1p(psi)) - 2 * b0


def _root(f, a, b, tolerance=1e-12):
    """Bracketed bisection: deliberately simple and transparent."""
    fa, fb = f(a), f(b)
    if not all(isfinite(v) for v in (fa, fb)) or fa * fb > 0:
        raise ValueError("Root endpoints must be finite and bracket a root")
    if fa == 0:
        return a
    if fb == 0:
        return b
    for _ in range(100):
        c = (a + b) / 2
        fc = f(c)
        if not isfinite(fc):
            raise ArithmeticError("Non-finite residual inside root bracket")
        if fc == 0 or b - a <= tolerance:
            return c
        if fa * fc < 0:
            b = c
        else:
            a, fa = c, fc
    raise ArithmeticError("Root iteration limit reached")


def _scan_root(f, a, b, count=1000):
    """Locate the first finite sign change, matching the reference scans."""
    previous = None
    for i in range(count):
        x = a + (b - a) * i / (count - 1)
        try:
            fx = f(x)
        except (ValueError, ZeroDivisionError, OverflowError):
            fx = float('nan')
        if isfinite(fx):
            if fx == 0:
                return x
            if previous is not None and previous[1] * fx < 0:
                return _root(f, previous[0], x)
            previous = (x, fx)
        else:
            previous = None
    raise ValueError(f"No finite sign change in [{a}, {b}]")


def solve_Qt(alpha, b0):
    _parameters(alpha, b0)
    return _scan_root(lambda Q: tip_residual(Q, alpha, b0), .05, 20., 4000)


def psi_from_geometry(Q, r, z, alpha, b0):
    bq = b_Q(Q, alpha, b0)
    return 2 * bq * z / (1 / r**2 - 2 * bq * z)


def branch_derivative(Q, S, r, z, alpha, b0):
    psi = psi_from_geometry(Q, r, z, alpha, b0)
    if not -1 < psi < 0:
        raise ArithmeticError("Geometric psi must lie between -1 and 0")
    bq = b_Q(Q, alpha, b0)
    # bQQ/bQ = -alpha for the exponential buoyancy profile.
    denominator = bq / sqrt(Q) + alpha / 2 * (psi - log1p(psi))
    return (bq + 1 / S) * r**2 * psi / denominator


def front_rhs(z, state, alpha, b0):
    Qm, Qp, r = state
    if not (0 < Qm < Qp and r > 0 and all(isfinite(v) for v in state)):
        raise ArithmeticError("Invalid front state")
    S = jump_slope(Qm, Qp, alpha, b0)
    result = (branch_derivative(Qm, S, r, z, alpha, b0),
              branch_derivative(Qp, S, r, z, alpha, b0), r**3 / S)
    if not all(isfinite(v) for v in result):
        raise ArithmeticError("Non-finite front derivative")
    return result


def lid_radius(Q, alpha, b0):
    """Non-frontal Q-surface intersection with z=1."""
    bq = b_Q(Q, alpha, b0)
    if bq == 0:
        return sqrt(2 * sqrt(Q))
    return sqrt(-expm1(-4 * bq * sqrt(Q)) / (2 * bq))


def launch_front(alpha, b0, Zt, dZ_launch=1e-5):
    """Reference Phase A launch and two fixed-tip-slope warmup steps."""
    _parameters(alpha, b0)
    if not (isfinite(Zt) and isfinite(dZ_launch)
            and 0 < Zt < 1 and 0 < dZ_launch < (1 - Zt) / 3):
        raise ValueError("Require 0 < Zt < Zt + 3*dZ_launch < 1")
    Qt = solve_Qt(alpha, b0)
    psit = psi_of_Q(Qt, alpha, b0)
    bqt = b_Q(Qt, alpha, b0)
    S = -1 / bqt
    rt = sqrt(1 / (2 * bqt * Zt * (1 + 1 / psit)))

    def advance(Qn, side, r, z, span):
        def residual(Q):
            dQdz = branch_derivative(Q, S, r, z, alpha, b0)
            return abs((Q - Qn) / dQdz) - dZ_launch
        if side == 'minus':
            return _scan_root(residual, max(1e-10, Qn - span), Qn - 1e-12)
        return _scan_root(residual, Qn + 1e-12, Qn + span)

    Qm = advance(Qt, 'minus', rt, Zt, .1)
    Qp = advance(Qt, 'plus', rt, Zt, .1)
    R = -1 / (2 * rt**2) + dZ_launch / S
    z = Zt + dZ_launch
    r = sqrt(-1 / (2 * R))
    warmup = [(z, Qm, Qp, r)]
    for _ in range(2):
        R += dZ_launch / S
        z += dZ_launch
        if R >= 0:
            raise ArithmeticError("Warmup front radius became unbounded")
        r = sqrt(-1 / (2 * R))
        Qm = advance(Qm, 'minus', r, z, .2)
        Qp = advance(Qp, 'plus', r, z, .2)
        warmup.append((z, Qm, Qp, r))
    return {'Qt': Qt, 'psit': psit, 'rf_t': rt, 'Zt': Zt,
            'warmup': warmup, 'state': (Qm, Qp, r), 'z_start': z}


# Dormand--Prince 5(4) tableau, the familiar RK45 pair.
_C = (0, 1/5, 3/10, 4/5, 8/9, 1, 1)
_A = ((), (1/5,), (3/40, 9/40), (44/45, -56/15, 32/9),
      (19372/6561, -25360/2187, 64448/6561, -212/729),
      (9017/3168, -355/33, 46732/5247, 49/176, -5103/18656),
      (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84))
_B4 = (5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40)


def _rk45(rhs, z, y, rtol, atol):
    rows = [(z, *y)]
    h = min(1e-5, (1 - z) / 100)
    accepted = rejected = 0
    for _ in range(200000):
        if z == 1:
            return rows, {'accepted_steps': accepted, 'rejected_steps': rejected}
        h = min(h, 1 - z)
        if h < 1e-14 * max(1, abs(z)):
            raise ArithmeticError(f"RK45 step too small at z={z}")
        try:
            stages = []
            for c, weights in zip(_C, _A):
                trial = tuple(y[j] + h * fsum(w * k[j] for w, k in zip(weights, stages))
                              for j in range(3))
                stages.append(rhs(z + c * h, trial))
            y5 = tuple(y[j] + h * fsum(w * k[j] for w, k in zip(_A[-1], stages))
                       for j in range(3))
            err = max(abs(h * fsum((a - b) * k[j]
                                  for a, b, k in zip((*_A[-1], 0), _B4, stages)))
                      / (atol + rtol * max(abs(y[j]), abs(y5[j]))) for j in range(3))
        except (ArithmeticError, ValueError):
            err = float('inf')
        if err <= 1:
            z = min(1., z + h)
            y = y5
            rows.append((z, *y))
            accepted += 1
        else:
            rejected += 1
        h *= 5 if err == 0 else max(.1, min(5., .9 * err**(-.2)))
    raise ArithmeticError("RK45 exceeded its step limit")


def integrate_front(alpha, b0, Zt, dZ_launch=1e-5, *,
                    rtol=1e-8, atol=1e-10, method='RK45'):
    """Integrate to the upper lid; failure never counts as a valid shot."""
    if not all(isfinite(v) and v > 0 for v in (rtol, atol)):
        raise ValueError("Tolerances must be finite and positive")
    launch = launch_front(alpha, b0, Zt, dZ_launch)
    rhs = lambda z, y: front_rhs(z, y, alpha, b0)
    if method == 'RK45':
        rows, diagnostics = _rk45(rhs, launch['z_start'], launch['state'], rtol, atol)
    elif method in ('BDF', 'Radau'):
        from scipy.integrate import solve_ivp
        ode = solve_ivp(rhs, (launch['z_start'], 1.), launch['state'],
                        method=method, rtol=rtol, atol=atol)
        if not ode.success or ode.t[-1] != 1.:
            raise ArithmeticError(f"{method} did not reach the upper lid: {ode.message}")
        rows = [(float(z), *(float(v) for v in state)) for z, state in zip(ode.t, ode.y.T)]
        diagnostics = {'nfev': ode.nfev, 'njev': ode.njev, 'nlu': ode.nlu}
    else:
        raise ValueError("method must be RK45, BDF, or Radau")
    sol = {key: [row[i] for row in rows] for i, key in enumerate(('z', 'Qm', 'Qp', 'rf'))}
    sol.update({key: launch[key] for key in ('Qt', 'psit', 'rf_t', 'Zt', 'warmup')})
    sol.update(alpha=alpha, b0=b0, has_front=True, method=method,
               rtol=rtol, atol=atol, dZ_launch=dZ_launch, diagnostics=diagnostics)
    sol['S'] = [jump_slope(qm, qp, alpha, b0) for qm, qp in zip(sol['Qm'], sol['Qp'])]
    for side in ('m', 'p'):
        sol['psi_' + side] = [psi_from_geometry(q, r, z, alpha, b0)
                              for q, r, z in zip(sol['Q' + side], sol['rf'], sol['z'])]
    return sol


def solve_zero_pv_vortex_front(alpha, b0, Zlo, Zhi, dZ_launch=1e-5, *,
                               tip_tolerance=1e-8, verbose=True, **integration_options):
    """Shoot on tip height until the two non-frontal lid radii agree.

    Returns (shot, sol), mirroring MATLAB's two outputs. The initial bracket
    must straddle a root. Invalid or incomplete integrations raise an error.
    """
    if not 0 < Zlo < Zhi < 1:
        raise ValueError("Require 0 < Zlo < Zhi < 1")
    if not isfinite(tip_tolerance) or tip_tolerance <= 0:
        raise ValueError("tip_tolerance must be finite and positive")
    cache = {}

    def residual(Zt):
        if Zt not in cache:
            sol = integrate_front(alpha, b0, Zt, dZ_launch, **integration_options)
            rm = lid_radius(sol['Qm'][-1], alpha, b0)
            rp = lid_radius(sol['Qp'][-1], alpha, b0)
            cache[Zt] = (rp - rm, sol, rm, rp)
            if verbose:
                print(f"Zt={Zt:.10f}  lid-radius residual={rp-rm:+.6e}", flush=True)
        return cache[Zt][0]

    Zt = _root(residual, Zlo, Zhi, tip_tolerance)
    residual(Zt)
    mismatch, sol, rm, rp = cache[Zt]
    shot = {'Zt': Zt, 'residual': mismatch, 'rlid_m': rm, 'rlid_p': rp,
            'front_lid_radius': sol['rf'][-1], 'shooting_evaluations': len(cache),
            'tip_tolerance': tip_tolerance}
    return shot, sol
