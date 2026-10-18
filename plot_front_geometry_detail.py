"""Enlarged physical-space view of the discrete geometric-model front.

Uses saved element vertices, rather than the ODE solution. Ordinary boundaries
are grey, nonconsecutive-label contacts are black, and chain endpoints are
labelled T (tip) and A (upper-lid attachment). Exports vector SVG and PDF.
"""
import argparse
import json
import math
from pathlib import Path

from coordinate_maps import r_from_R
from front_geometry import frontal_edges, check_front_chain, front_summary
from plot_geometric_velocity import sample_edge
from plot_ou_paper_cases import VectorFigure


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result', type=Path, nargs='?', default=Path(__file__).with_name(
        'paper_comparison') / 'dry_b2_n20_r0_0.003_rmax_4.2.json')
    parser.add_argument('--output', type=Path, default=root / 'output/pdf' /
                        'Geometric_b2_N20_upper_front_geometry.svg')
    args = parser.parse_args()
    result = json.loads(args.result.read_text())
    p = result['parameters']
    if not result['converged'] or result['empty_cells']:
        raise ValueError('This figure requires a converged run with no empty cells.')
    if p['coefficients'] != [1., 0., 0.] or p['H0'] != 'Infinity':
        raise ValueError('This physical-height view expects the dry incompressible case.')
    edges = frontal_edges(result)
    chain = check_front_chain(edges)
    summary = front_summary(result)
    r0 = p['r_inner']
    nodes = []
    for edge in edges:
        for point in (edge['start'], edge['end']):
            if not any(math.dist(point, other) < 1e-8 for other in nodes):
                nodes.append(point)
    nodes.sort(key=lambda point: point[1])
    tip, attachment = nodes[0], nodes[-1]
    if abs(attachment[1]-1.) > 1e-9:
        raise ValueError('Expected the front endpoint to lie on the upper lid.')

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig = VectorFigure(args.output, 1100, 690, True)
    fig.canvas.setTitle('Geometric model: enlarged upper-front element geometry')
    left, top, width, depth = 100., 120., 910., 420.
    rmin, rmax, zmin, zmax = 1.7, 2.8, .72, 1.
    xy = lambda r, z: (left+(r-rmin)/(rmax-rmin)*width,
                       top+(zmax-z)/(zmax-zmin)*depth)
    mapped = lambda point: xy(r_from_R(point[0], r0), point[1])
    fig.text(550, 32, f'Upper-front element geometry: {p["n"]} elements', 24)
    fig.text(550, 61, f'b0 = {p["bm0"]:g}, alpha = {p["alpha"]:g}; '
             f'r0 = {r0:g}, rmax = {p["r_outer"]:g}; dry, incompressible', 16)

    fig.clip(left, top, width, depth)
    # Draw each shared edge once. All boundaries retain their computed geometry.
    seen = set()
    for polygon in result['polygons']:
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            key = tuple(sorted((tuple(round(v, 8) for v in a),
                                tuple(round(v, 8) for v in b))))
            if key in seen:
                continue
            seen.add(key)
            fig.line([xy(r, z) for r, z in sample_edge(a, b, r0)],
                     width=.75, colour=(.64, .64, .64))
    for edge in edges:
        fig.line([xy(r, z) for r, z in sample_edge(edge['start'], edge['end'], r0)],
                 width=3.)
    fig.endclip()

    fig.line([xy(rmin, zmin), xy(rmax, zmin), xy(rmax, 1.),
              xy(rmin, 1.), xy(rmin, zmin)], width=1.)
    for r in (1.8, 2., 2.2, 2.4, 2.6, 2.8):
        x, y = xy(r, zmin)
        fig.line([(x, y), (x, y+5)])
        fig.text(x, y+24, f'{r:g}', 16)
    for z in (.75, .8, .85, .9, .95, 1.):
        x, y = xy(rmin, z)
        fig.line([(x-5, y), (x, y)])
        fig.text(x-13, y+5, f'{z:g}', 16, 'end')
    fig.text(left+width/2, top+depth+52, 'r', 20)
    fig.text(33, top+depth/2, 'z', 20)

    def marker(point, radius):
        x, y = mapped(point)
        circle = [(x+radius*math.cos(2*math.pi*i/40),
                   y+radius*math.sin(2*math.pi*i/40)) for i in range(41)]
        fig.polygon(circle, (1., 1., 1.))
        fig.line(circle, width=1.5)

    for node in nodes:
        marker(node, 3.5)
    for node, letter, dx, dy in ((tip, 'T', -23, -9),
                                  (attachment, 'A', 15, -13)):
        marker(node, 6.)
        x, y = mapped(node)
        fig.text(x+dx, y+dy, letter, 20)
    fig.text(*xy(2.72, 1.017), 'Upper lid', 15)
    fig.text(115, 623, f'T: tip ({summary["tip_r"]:.5f}, {summary["tip_z"]:.5f})', 17, 'start')
    fig.text(590, 623, f'A: attachment ({summary["upper_r"]:.5f}, {summary["upper_z"]:.5f})', 17, 'start')
    fig.line([(115, 655), (165, 655)], width=.75, colour=(.64, .64, .64))
    fig.text(180, 660, 'Element edges', 15, 'start')
    fig.line([(440, 655), (490, 655)], width=3.)
    fig.text(505, 660, 'Linked front edges', 15, 'start')
    fig.text(810, 660, 'Circles: front vertices', 15, 'start')
    fig.save()
    metadata = {'source': str(args.result.resolve()), 'parameters': p,
                'max_mass_error': result['max_mass_error'],
                'front_chain': chain, 'front': summary,
                'front_vertices_physical': [[r_from_R(a, r0), z] for a, z in nodes],
                'window': {'r': [rmin, rmax], 'z': [zmin, zmax]},
                'coordinates': 'Physical nondimensional r and z; axes independently scaled',
                'ordinary_edge_width': .75, 'front_edge_width': 3.,
                'front_label_pairs_one_based': [[e['i']+1, e['j']+1] for e in edges]}
    args.output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
