"""Small test: eight cells, three formulations, generated entirely in Python."""
import json
import math
from pathlib import Path
from coordinate_maps import DRY_COEFFICIENTS
from equilibrium_solver import solve_zero_pv_equilibrium

cases = [
    ("dry", dict(coefficients=DRY_COEFFICIENTS,H0=math.inf)),
    ("density", dict(coefficients=DRY_COEFFICIENTS,H0=0.75)),
    ("moist", dict(H0=0.75)),
]

# Strict JSON does not have infinity. Encode that parameter as descriptive text.
def json_safe(value):
    if isinstance(value,dict):
        return {k:json_safe(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [json_safe(v) for v in value]
    if isinstance(value,float) and value == math.inf:
        return "Infinity"
    return value


if __name__ == "__main__":
    for name,parameters in cases:
        print(f"\n=== {name.upper()} CASE: n=8, r0=0.19, rmax=3, bm0=2, alpha=1 ===")
        result = solve_zero_pv_equilibrium(8,0.19,3.0,2.0,1.0, **parameters)
        print(f"Converged: {result['converged']}")
        print(f"Total mass: {result['total_mass']:.12f}")
        print(f"Target mass: {result['total_target_mass']:.12f}")
        print("Empty cells:",result["empty_cells"])
        if not result["converged"]:
            raise RuntimeError(f"{name} case failed: {result['status']}")
        output = Path(__file__).with_name(f"equilibrium_small_{name}.json")
        output.write_text(json.dumps(json_safe(result),indent=2,allow_nan=False))
        print("Saved:",output)
