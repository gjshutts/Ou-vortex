"""Run the dry Ou front solver and save a readable JSON result.

Example: python3 run_ou_front.py --b0 7.6
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

from ou_front_solver import solve_zero_pv_vortex_front


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b0', type=float, default=7.6)
    parser.add_argument('--alpha', type=float, default=1.)
    parser.add_argument('--zlo', type=float)
    parser.add_argument('--zhi', type=float)
    parser.add_argument('--launch', type=float, default=1e-5)
    parser.add_argument('--rtol', type=float, default=1e-8)
    parser.add_argument('--atol', type=float, default=1e-10)
    parser.add_argument('--tip-tolerance', type=float, default=1e-8)
    parser.add_argument('--method', choices=('RK45', 'BDF', 'Radau'), default='RK45')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if (args.zlo is None) != (args.zhi is None):
        parser.error('Supply both --zlo and --zhi')
    if args.zlo is None:
        brackets = {2.: (.75, .85), 7.6: (.68, .70), 16.: (.68, .70)}
        if args.alpha != 1. or args.b0 not in brackets:
            parser.error('Supply --zlo and --zhi for this parameter combination')
        zlo, zhi = brackets[args.b0]
    else:
        zlo, zhi = args.zlo, args.zhi
    start = perf_counter()
    shot, sol = solve_zero_pv_vortex_front(
        args.alpha, args.b0, zlo, zhi, args.launch, method=args.method,
        rtol=args.rtol, atol=args.atol, tip_tolerance=args.tip_tolerance)
    elapsed = perf_counter() - start
    output = args.output or Path(__file__).parent / f'ou_front_b{args.b0:g}.json'
    output.write_text(json.dumps({'shot': shot, 'sol': sol, 'elapsed_seconds': elapsed},
                                 indent=2, allow_nan=False) + '\n')
    print(f"\nTip: z={shot['Zt']:.10f}, r={sol['rf_t']:.10f}")
    print(f"Lid: r={shot['front_lid_radius']:.10f}; residual={shot['residual']:+.3e}")
    print(f"Elapsed: {elapsed:.3f} seconds\nSaved: {output.resolve()}")


if __name__ == '__main__':
    main()
