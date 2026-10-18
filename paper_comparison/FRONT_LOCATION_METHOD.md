# Discrete front locations and the 80-element comparison

Both geometric runs use b0=2, alpha=1, r0=0.003, rmax=4.2 and constant density.
The fresh 80-element run converged in 43 Newton steps (3.58 seconds), with
maximum mass error 7.7885e-7. The 400-element iterate has maximum mass error
1.7100e-6 and does not satisfy the chosen 1e-6 criterion.

For the ordered zero-PV label family, ordinary material neighbours have
consecutive Q labels. `front_geometry.frontal_edges` finds shared polygon
edges between nonconsecutive labels and verifies their plane-equality residual.
The selected edges form one connected, unbranched chain with two endpoints
in both cases. This is a model-specific front detector, not a generic rule
for arbitrary PV distributions or an ocean lens.

Each contact lies on

    (Qi-Qj) R + (bi-bj) z + (Si-Sj) = 0.

The discrete tip is the chain's lowest endpoint, a polygon junction where the
first such contact begins. The attachment is the other endpoint, lying at
z=1 where the last contact intersects the lid. Coordinates come from the
computed polygons and are mapped back to physical r. No image threshold,
colour contour, fitted curve or manually chosen point enters this extraction.

A finite-element tip is not the continuum zero-jump critical point: the two
labels on its first edge still have a finite separation. Its resolution error
must therefore be assessed by changing the element count. The attachment is
also resolution-dependent even though its intersection with the lid is
well-defined for any particular discrete solution.

| Model | Tip r | Tip z | Lid r |
|---|---:|---:|---:|
| ODE | 2.002812667 | 0.799619862 | 2.433723308 |
| 80 elements | 2.004642305 | 0.800699749 | 2.433482445 |
| 400-element final iterate | 2.002875734 | 0.799657971 | 2.433706227 |

Both new plots omit element outlines and keep the black discrete front.
They use the same aspect ratio and colour limits +/-1.05. Within each cell,
v=sqrt(Qi)/r-r/2; the analytic ODE field is not substituted into the plot.
Subpixel overlap of gradient fills removes white antialias seams during
vector rendering. That rendering aid is independent of the raw polygon
coordinates used to extract and draw the front.

SVG and PDF variants are under output/pdf, with names ending in
`velocity_no_edges`. The 80-element field visibly retains its discrete
velocity bands, while the 400-element field is much smoother.
