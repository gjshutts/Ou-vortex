# Geometric/ODE convergence comparison

Dry, constant-density equilibria at r0=0.003, rmax=4.2, b0=2, alpha=1.
Initial elements are uniformly spaced in physical radius. The same absolute
per-element mass tolerance, 1e-6, is used for every resolution.
Differences are geometric minus ODE; all coordinates are nondimensional.

| Ne | Tip r | Tip z | Attachment r | Delta tip r | Delta tip z | Delta attachment r | Max mass residual |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 2.031524272 | 0.814567720 | 2.428447031 | +2.871e-02 | +1.495e-02 | -5.276e-03 | 9.727e-07 |
| 40 | 2.008762486 | 0.803495815 | 2.433267610 | +5.950e-03 | +3.876e-03 | -4.557e-04 | 8.690e-07 |
| 80 | 2.004642305 | 0.800699749 | 2.433482445 | +1.830e-03 | +1.080e-03 | -2.409e-04 | 7.789e-07 |
| 160 | 2.003172829 | 0.799840615 | 2.433637344 | +3.602e-04 | +2.208e-04 | -8.596e-05 | 9.623e-07 |
| 320* | 2.002967641 | 0.799715732 | 2.433717405 | +1.550e-04 | +9.587e-05 | -5.903e-06 | 1.055e-06 |
| 400* | 2.002875734 | 0.799657971 | 2.433706227 | +6.307e-05 | +3.811e-05 | -1.708e-05 | 1.710e-06 |
| ODE | 2.002812667 | 0.799619862 | 2.433723308 | - | - | - | - |

* Asterisk: the run stopped without satisfying the 1e-6 mass tolerance; the
reported coordinates belong to its final iterate. All rows retain every
element and have a single connected, unbranched frontal chain.

Locations are extracted from polygon contacts, not measured from images.
The tip is the lowest endpoint of the nonconsecutive-label contact chain;
the attachment is its intersection with z=1. A finite-element tip is a
discrete approximation to the continuum zero-jump critical point.

The ODE reference uses launch increment 1e-5, integration relative tolerance
1e-10, absolute tolerance 1e-12, and shooting height tolerance 1e-11.
Earlier launch refinement moved its lid radius by about 2.2e-7. The ODE
extends towards the axis and has no finite outer wall; differences from
the fixed-domain geometric solutions can include boundary effects as well
as discretisation and algebraic errors. They should not be called pure
discretisation errors without a separate boundary-sensitivity argument.

Full precision, convergence flags and source filenames are in the companion
JSON. The companion .tex is a table fragment for inclusion in the paper;
it uses standard LaTeX tabular and needs no additional package.
