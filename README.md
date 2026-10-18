[README.md](https://github.com/user-attachments/files/33640325/README.md)
# Rotational Ou vortex: Python equilibrium solvers

**Publication preparation copy. This folder has not yet been uploaded to GitHub
or assigned a release DOI. A software licence and final paper citation remain
to be selected before publication.**

Two complementary calculations for a rotating, axisymmetric fluid initially at
rest, with zero potential vorticity and prescribed buoyancy
`b(Q) = b0 exp(-alpha Q)`, where `Q = M^2` is squared absolute angular momentum:

1. **Ou-style front calculation:** analytical expressions away from the front,
   coupled to numerical integration and shooting for the frontal ODE.
2. **Geometric element model:** a mass-constrained equilibrium constructed from
   the upper envelope of pressure planes in transformed coordinates.

Both core solvers use only Python's standard library. MATLAB and Mathematica
are not required to run these calculations. The geometric solver additionally
contains optional moist amplification and static density weighting; the paper
comparison cases supplied here are **dry and incompressible**.

## Getting started

Use Python 3.10 or later. Open a terminal in this repository folder, then:

```sh
cd code
python3 -m unittest discover -v
python3 run_ou_front.py --b0 2 --rtol 1e-10 --atol 1e-12 --tip-tolerance 1e-11
python3 run_geometric_paper_comparison.py --elements 20
```

The core numerical tests need no additional packages. One plotting-grid test
requires NumPy; install the plotting dependencies before running the complete
test suite:

```sh
python3 -m pip install -r ../requirements-plotting.txt
```

To run only the standard-library tests, use:

```sh
python3 -m unittest -v test_coordinate_maps test_cell_geometry test_mass_integration test_mass_jacobian test_equilibrium_solver test_uniform_exact test_ou_front_solver
```

See [REPRODUCING.md](REPRODUCING.md) for the paper cases, figure commands,
convergence table and numerical limitations. Files in `code/` import one
another directly: keep them together when learning from or reusing an example.

## Main files

| File in `code/` | Purpose |
|---|---|
| `coordinate_maps.py` | Radial/vertical coordinate maps and mass metric |
| `cell_geometry.py` | Convex cells obtained by half-plane clipping |
| `mass_integration.py` | Physical mass of a transformed polygon |
| `mass_jacobian.py` | Derivatives of element masses with respect to plane offsets |
| `initial_state.py` | Rest-state parcel labels and target masses |
| `equilibrium_solver.py` | Damped Newton equilibrium solver |
| `front_geometry.py` | Discrete frontal contacts and endpoints |
| `ou_front_solver.py` | Independent dry frontal ODE solver |
| `ou_fields.py` | Analytical field reconstruction using the computed ODE front |
| `plot_ou_paper_cases.py` | Buoyancy/velocity figures for b0 = 2, 4, 7 |
| `plot_geometric_velocity.py` | Element-model buoyancy or velocity plots |
| `plot_front_geometry_detail.py` | Enlarged 20-element front with marked vertices |
| `build_convergence_table.py` | Rebuild the saved front-location comparison |

The explanatory notes in `docs/` introduce the calculation building block by
building block. Their commands assume that `code/` is the working directory.
`code/paper_comparison/` contains Python-generated reference equilibria and
the convergence table; it includes the two final iterates that miss the
specified mass tolerance, clearly flagged as such.

## Citation and release

`CITATION.cff` supplies draft software citation metadata. The repository URL,
paper reference and archived release DOI will be added when available.
GitHub hosts the working source; a fixed release should be archived on Zenodo
and its version-specific DOI cited for reproducibility.

Maintainer: Glenn Shutts. Reuse terms have not yet been selected; see
[PUBLICATION_STATUS.md](PUBLICATION_STATUS.md).
