"""Small-case equilibrium solver: pinned offset, analytic Jacobian, backtracking.

Correctness-first implementation using standard Python, without acceleration.
No time stepping is performed. Q and bm remain fixed while S is corrected.
"""
import math
from coordinate_maps import DEFAULT_COEFFICIENTS, R_from_r
from cell_geometry import rectangle, cell_polygons
from mass_integration import polygon_mass, cell_masses
from mass_jacobian import mass_jacobian
from initial_state import generate_initial_state


def solve_linear_system(matrix, rhs, pivot_tolerance=1e-13):
    """Gaussian elimination with partial pivoting, for small dense systems.

    A later performance version can use a numerical-library linear solver.
    Work on a copy, choose the largest pivot in each column, eliminate below
    it, then substitute backwards. Never form a matrix inverse explicitly.
    """
    n = len(rhs)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("The matrix must be square and match rhs.")
    if n == 0:
        return []
    if not all(math.isfinite(v) for row in matrix for v in row) or not all(math.isfinite(v) for v in rhs):
        raise ValueError("Linear-system inputs must be finite.")
    rows = [list(row)+[value] for row,value in zip(matrix, rhs)]
    scale = max(abs(v) for row in matrix for v in row)
    if scale == 0:
        raise ArithmeticError("Singular reduced Jacobian.")
    for k in range(n):
        pivot = max(range(k,n), key=lambda i: abs(rows[i][k]))
        if abs(rows[pivot][k]) <= pivot_tolerance*scale:
            raise ArithmeticError("Reduced Jacobian is singular or poorly conditioned.")
        rows[k], rows[pivot] = rows[pivot], rows[k]
        for i in range(k+1,n):
            factor = rows[i][k]/rows[k][k]
            rows[i][k] = 0.0
            for j in range(k+1,n+1):
                rows[i][j] -= factor*rows[k][j]
    solution = [0.0]*n
    for i in range(n-1,-1,-1):
        solution[i] = (rows[i][n]-math.fsum(rows[i][j]*solution[j] for j in range(i+1,n)))/rows[i][i]
    return solution


def solve_offsets(Q, bm, target_masses, domain, initial_offsets,
                  r0=0.095, coefficients=DEFAULT_COEFFICIENTS, H0=0.75, rho_s=1.0,
                  max_iterations=100, mass_tolerance=1e-6, damping=0.25,
                  max_backtracks=20, quadrature_tolerance=1e-10, verbose=True):
    """Correct S so each current cell mass matches its prescribed target.

    Offset 0 is pinned to zero: subtract its initial value from every S.
    Acceptance requires positive masses and a lower maximum mass residual.
    Always inspect result['converged']; an unfinished iterate is not a solution.
    """
    n = len(Q)
    if n == 0 or any(len(v) != n for v in (bm, target_masses, initial_offsets)):
        raise ValueError("Supply equally sized, nonempty label, target and offset lists.")
    if not all(math.isfinite(v) for values in (Q,bm,target_masses,initial_offsets) for v in values):
        raise ValueError("Labels, targets and offsets must be finite.")
    if min(target_masses) <= 0:
        raise ValueError("Every prescribed material mass must be positive.")
    if not isinstance(max_iterations,int) or max_iterations < 0:
        raise ValueError("max_iterations must be a nonnegative integer.")
    if not isinstance(max_backtracks,int) or max_backtracks < 1:
        raise ValueError("max_backtracks must be a positive integer.")
    if not math.isfinite(mass_tolerance) or mass_tolerance <= 0:
        raise ValueError("mass_tolerance must be finite and positive.")
    if not math.isfinite(damping) or not 0 < damping <= 1:
        raise ValueError("damping must lie in (0,1].")
    parameters = dict(r0=r0, coefficients=coefficients, H0=H0, rho_s=rho_s)
    target_total = math.fsum(target_masses)
    domain_mass = polygon_mass(domain, absolute_tolerance=quadrature_tolerance, **parameters)
    if abs(domain_mass-target_total) > max(10*n*quadrature_tolerance, 1e-9*target_total):
        raise ValueError("The target masses do not sum to the domain's physical mass.")
    offsets = [S-initial_offsets[0] for S in initial_offsets]

    def evaluate(S):
        planes = list(zip(Q,bm,S))
        cells = cell_polygons(planes,domain)
        masses = cell_masses(cells,absolute_tolerance=quadrature_tolerance, **parameters)
        residual = [m-t for m,t in zip(masses,target_masses)]
        return planes,cells,masses,residual,max(abs(v) for v in residual)

    planes,cells,masses,residual,error = evaluate(offsets)
    history = []
    status = "iteration_limit"
    for iteration in range(max_iterations+1):
        record = {"iteration": iteration, "max_mass_error": error,
                  "total_mass_error": math.fsum(residual), "minimum_mass": min(masses)}
        history.append(record)
        if verbose:
            print(f"Newton {iteration:3d}: max mass error={error:.6e}, "
                  f"total error={record['total_mass_error']:.3e}", flush=True)
        if min(masses) <= 0:
            status = "empty_initial_cell"
            break
        if error < mass_tolerance:
            status = "converged"
            break
        if iteration == max_iterations:
            break
        J = mass_jacobian(planes,cells, **parameters)
        try:
            correction = [0.0]+solve_linear_system([row[1:] for row in J[1:]],
                                                   [-v for v in residual[1:]])
        except ArithmeticError:
            status = "singular_jacobian"
            break
        step = damping
        accepted = False
        for _ in range(max_backtracks):
            trial = [S+step*dS for S,dS in zip(offsets,correction)]
            try:
                trial_state = evaluate(trial)
                if min(trial_state[2]) > 0 and trial_state[4] < error:
                    accepted = True
                    break
            except ArithmeticError:
                pass  # A quadrature failure cannot justify accepting a trial.
            step /= 2
        if not accepted:
            status = "no_improving_step"
            break
        record["accepted_step"] = step
        if verbose and step < damping:
            print(f"            reduced step={step:.6g}", flush=True)
        offsets = trial
        planes,cells,masses,residual,error = trial_state

    if verbose:
        print(f"Status: {status}; accepted Newton steps: {len(history)-1}", flush=True)
    return {"converged": status == "converged", "status": status,
            "iterations": len(history)-1, "offsets": offsets, "planes": planes,
            "polygons": cells, "masses": masses, "target_masses": list(target_masses),
            "mass_errors": residual, "max_mass_error": error,
            "max_relative_mass_error": max(abs(v)/t for v,t in zip(residual,target_masses)),
            "total_mass": math.fsum(masses), "total_target_mass": target_total,
            "empty_cells": [i for i,m in enumerate(masses) if m <= 0], "history": history}


def solve_zero_pv_equilibrium(n=8, r_inner=0.19, r_outer=3.0, bm0=2.0, alpha=1.0,
                              coefficients=DEFAULT_COEFFICIENTS, H0=0.75, rho_s=1.0,
                              initial_offsets=None, **solver_options):
    """Generate rest-state labels and targets, then solve a unit-depth vortex."""
    initial = generate_initial_state(n,r_inner,r_outer,bm0,alpha,H0,rho_s)
    start = initial["initial_offsets"] if initial_offsets is None else initial_offsets
    domain = rectangle(0,R_from_r(r_outer,r_inner))
    result = solve_offsets(initial["Q"],initial["bm"],initial["target_masses"],domain,start,
                           r0=r_inner,coefficients=coefficients,H0=H0,rho_s=rho_s, **solver_options)
    result["initial_state"] = initial
    result["parameters"] = dict(n=n,r_inner=r_inner,r_outer=r_outer,bm0=bm0,alpha=alpha,
                                 coefficients=list(coefficients),H0=H0,rho_s=rho_s)
    return result
