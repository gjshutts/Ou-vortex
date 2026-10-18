# Rotational Ou vortex: Python equilibrium solvers

Python code for extending the rectilinear adjustment problem of Ou (1984)
to rotating, axisymmetric flow. The initial state is at rest, with zero
potential vorticity and prescribed buoyancy

```text
b(Q) = b0 exp(-alpha Q),     Q = M^2,
```

where M is absolute angular momentum. All coordinates and parameters used
in the supplied examples are nondimensional.

## Two complementary calculations

1. **Ou-style front calculation:** analytical expressions away from the front,
   coupled to numerical integration and shooting for the frontal ODE.
2. **Geometric element model:** a mass-constrained equilibrium constructed from
   the upper envelope of pressure planes in transformed coordinates.

Both core solvers use only Python's standard library. MATLAB and Mathematica
are not required. The geometric solver also contains optional moist
amplification and static density weighting; the paper comparison cases
supplied here are **dry and incompressible**.

These programs calculate idealized balanced equilibria. They do not simulate
the time-dependent turbulent adjustment or the finite thickness of a real
frontal transition zone.

## Getting started

Use Python 3.10 or later. Download or clone this repository, then open a
terminal in the downloaded **Ou-vortex folder**. The Python files are directly
in that folder; there is no additional `code` subfolder to enter.

Run a dry frontal ODE calculation:

```sh
python3 run_ou_front.py --b0 2 --rtol 1e-10 --atol 1e-12 --tip-tolerance 1e-11
```

Run a 20-element geometric equilibrium:

```sh
python3 run_geometric_paper_comparison.py --elements 20
```

This uses r0 = 0.003, rmax = 4.2, b0 = 2 and alpha = 1. It writes results to
`paper_comparison/`, replacing the corresponding saved result.

## Verification

The following tests require only Python's standard library:

```sh
python3 -m unittest -v test_coordinate_maps test_cell_geometry test_mass_integration test_mass_jacobian test_equilibrium_solver test_uniform_exact test_ou_front_solver
```

For the complete test suite and plotting routines, install NumPy, Pillow
and ReportLab. On macOS or Linux, an isolated Python environment can be set
up from the repository folder with:

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy Pillow reportlab
python3 -m unittest discover -v
```

The Python publication preparation copy passed all 52 tests, and a fresh
20-element calculation reproduced its saved front coordinates. The default
ODE integrator is RK45; optional SciPy BDF/Radau paths are not the basis of
the supplied results.

## Reproducing the figures

With the plotting packages installed, run the ODE paper cases b0 = 2, 4 and 7:

```sh
python3 plot_ou_paper_cases.py --pdf --output output/pdf
```

This produces paired buoyancy and tangential-velocity figures in SVG and PDF,
plus an overview. Shading is embedded raster data; contours, fronts and text
are vector objects. The plotting window extends from r = 0 to 4.2; it does
not impose a radial wall on the ODE solution.

Plot the saved 80-element geometric solution with element edges:

```sh
python3 plot_geometric_velocity.py paper_comparison/dry_b2_n80_r0_0.003_rmax_4.2.json --field v --edge-width 0.25 --pdf --output output/pdf/Geometric_b2_N80_velocity_with_edges.svg
python3 plot_geometric_velocity.py paper_comparison/dry_b2_n80_r0_0.003_rmax_4.2.json --field b --edge-width 0.25 --pdf --output output/pdf/Geometric_b2_N80_buoyancy_with_edges.svg
```

Produce the enlarged 20-element front geometry:

```sh
python3 plot_front_geometry_detail.py --output output/pdf/Geometric_b2_N20_upper_front_geometry.svg
```

Six bold edges join seven marked frontal vertices. T identifies the tip;
A identifies the upper-lid attachment. These plots use the geometric-model
solution, not ODE field values.

## Front-location convergence table

The supplied [comparison table](paper_comparison/geometric_ode_convergence_b2.md)
contains 20-, 40-, 80-, 160-, 320- and 400-element results at fixed parameters.
Rebuild it from the saved Python results with:

```sh
python3 build_convergence_table.py
```

The 320- and 400-element final iterates narrowly miss the absolute element-mass
tolerance of 1e-6; both rows are explicitly flagged. They must not be treated
as converged runs. The geometric model has a finite radial wall whereas the
ODE construction does not, so differences can include boundary effects as
well as discretisation and algebraic errors.

## Main files

| File | Purpose |
|---|---|
| `coordinate_maps.py` | Coordinate maps and mass metric |
| `cell_geometry.py` | Convex cells from half-plane clipping |
| `mass_integration.py` | Physical mass of a transformed polygon |
| `mass_jacobian.py` | Derivatives of masses with respect to plane offsets |
| `initial_state.py` | Rest-state parcel labels and target masses |
| `equilibrium_solver.py` | Damped Newton equilibrium solver |
| `front_geometry.py` | Discrete frontal contacts and endpoints |
| `ou_front_solver.py` | Independent dry frontal ODE solver |
| `ou_fields.py` | Field reconstruction using the computed ODE front |
| `plot_ou_paper_cases.py` | ODE buoyancy and velocity figures |
| `plot_geometric_velocity.py` | Geometric-model buoyancy or velocity figures |
| `plot_front_geometry_detail.py` | Enlarged front with marked vertices |
| `build_convergence_table.py` | Front-location comparison table |

Example scripts introduce the calculation in stages: coordinate maps, cell
geometry, mass integration, the mass Jacobian and the equilibrium solver.
Keep the Python modules together so their local imports work.

## Author, licence and citation

Maintainer: Glenn Shutts.

The code is distributed under the [MIT licence](LICENSE).
The accompanying paper citation and a DOI for a fixed software release will
be added when available. The working repository is
[gjshutts/Ou-vortex](https://github.com/gjshutts/Ou-vortex).
