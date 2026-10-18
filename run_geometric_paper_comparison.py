"""Dry geometric-model run matching the ODE paper plotting window.

Default: N=400, r0=0.003, rmax=4.2, b0=2, alpha=1, density=1.
The mass tolerance and Newton method remain those of the reference solver.
"""
import argparse
import json
import math
from pathlib import Path
from time import perf_counter

from coordinate_maps import DRY_COEFFICIENTS, r_from_R
from equilibrium_solver import solve_zero_pv_equilibrium
from front_geometry import frontal_edges, front_summary, check_front_chain
from ou_front_solver import solve_zero_pv_vortex_front
from ou_fields import OuFields


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--elements',type=int,default=400)
    parser.add_argument('--r0',type=float,default=.003)
    args=parser.parse_args()
    root=Path(__file__).parent/'paper_comparison'
    root.mkdir(exist_ok=True)
    stem=f'dry_b2_n{args.elements}_r0_{args.r0:g}_rmax_4.2'
    start=perf_counter()
    result=solve_zero_pv_equilibrium(args.elements,args.r0,4.2,2.,1.,
                coefficients=DRY_COEFFICIENTS,H0=math.inf,max_iterations=150,
                mass_tolerance=1e-6,verbose=True)
    elapsed=perf_counter()-start
    result['parameters']['H0']='Infinity'
    result['elapsed_seconds']=elapsed
    (root/(stem+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    summary={k:result[k] for k in ('converged','status','iterations','max_mass_error',
                                  'max_relative_mass_error','empty_cells')}
    summary.update(parameters=result['parameters'],elapsed_seconds=elapsed,
                   total_mass_error=result['total_mass']-result['total_target_mass'])
    if result['converged']:
        try:
            edges=frontal_edges(result)
            summary['geometric_front']=front_summary(result)
            summary['front_chain']=check_front_chain(edges)
            shot,ode=solve_zero_pv_vortex_front(1.,2.,.75,.85,verbose=False,
                          rtol=1e-10,atol=1e-12,tip_tolerance=1e-11)
            (root/'ode_b2_reference.json').write_text(json.dumps({'shot':shot,'sol':ode},indent=2)+'\n')
            summary['ode_front']={'tip_r':ode['rf_t'],'tip_z':ode['Zt'],
                                  'upper_r':ode['rf'][-1]}
            summary['geometric_minus_ode']={key:summary['geometric_front'][key]-value
                                    for key,value in summary['ode_front'].items()}
            field=OuFields(ode)
            samples=[]
            for z in (.81,.85,.9,.95,1.):
                for edge in edges:
                    a,b=edge['start'],edge['end']
                    if min(a[1],b[1])-1e-12<=z<=max(a[1],b[1])+1e-12 and a[1]!=b[1]:
                        t=(z-a[1])/(b[1]-a[1]);R=a[0]+t*(b[0]-a[0])
                        r=r_from_R(R,args.r0);ro=field.front_r(z)
                        samples.append({'z':z,'geometric_r':r,'ode_r':ro,'difference':r-ro})
                        break
            summary['front_samples']=samples
        except (ValueError,IndexError) as exc:
            summary['front_comparison_error']=str(exc)
    (root/(stem+'_summary.json')).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('FINAL SUMMARY',json.dumps(summary),flush=True)
    print('Saved to',root.resolve(),flush=True)


if __name__=='__main__':main()
