"""Build the fixed-domain b0=2 geometric/ODE front-location table.

Read saved results; never re-solve or silently discard an unconverged row.
Outputs a Markdown table, a LaTeX table fragment and full-precision JSON.
"""
import json
from pathlib import Path

from front_geometry import front_summary,frontal_edges,check_front_chain


def main():
    folder=Path(__file__).parent/'paper_comparison'
    ode=json.loads((folder/'ode_b2_reference.json').read_text())['sol']
    assert ode['alpha']==1. and ode['b0']==2. and ode['z'][-1]==1.
    reference={'tip_r':ode['rf_t'],'tip_z':ode['Zt'],'attachment_r':ode['rf'][-1]}
    rows=[]
    for n in (20,40,80,160,320,400):
        filename=f'dry_b2_n{n}_r0_0.003_rmax_4.2.json'
        result=json.loads((folder/filename).read_text());p=result['parameters']
        assert (p['n'],p['r_inner'],p['r_outer'],p['bm0'],p['alpha'])==(n,.003,4.2,2.,1.)
        assert p['H0']=='Infinity' and p['coefficients']==[1.,0.,0.] and p['rho_s']==1.
        assert not result['empty_cells']
        assert result['converged']==(result['max_mass_error']<1e-6)
        geometry=front_summary(result)
        assert abs(geometry['upper_z']-1.)<1e-10
        chain=check_front_chain(frontal_edges(result))
        values={'tip_r':geometry['tip_r'],'tip_z':geometry['tip_z'],'attachment_r':geometry['upper_r']}
        row=dict(n=n,**values,differences={k:values[k]-v for k,v in reference.items()},
                 max_mass_error=result['max_mass_error'],max_relative_mass_error=result['max_relative_mass_error'],
                 converged=result['converged'],status=result['status'],iterations=result['iterations'],
                 front_segments=geometry['segments'],front_chain=chain,source=filename)
        rows.append(row)
    data={'parameters':{'r0':.003,'rmax':4.2,'b0':2.,'alpha':1.,'density':1.,'mass_tolerance':1e-6},
          'difference_convention':'geometric minus ODE',
          'ode_reference':reference,'rows':rows}
    (folder/'geometric_ode_convergence_b2.json').write_text(json.dumps(data,indent=2)+'\n')
    head='''# Geometric/ODE convergence comparison

Dry, constant-density equilibria at r0=0.003, rmax=4.2, b0=2, alpha=1.
Initial elements are uniformly spaced in physical radius. The same absolute
per-element mass tolerance, 1e-6, is used for every resolution.
Differences are geometric minus ODE; all coordinates are nondimensional.

| Ne | Tip r | Tip z | Attachment r | Delta tip r | Delta tip z | Delta attachment r | Max mass residual |
|---:|---:|---:|---:|---:|---:|---:|---:|
'''
    lines=[]
    for row in rows:
        tag='' if row['converged'] else '*'
        d=row['differences']
        lines.append(f"| {row['n']}{tag} | {row['tip_r']:.9f} | {row['tip_z']:.9f} | {row['attachment_r']:.9f} | {d['tip_r']:+.3e} | {d['tip_z']:+.3e} | {d['attachment_r']:+.3e} | {row['max_mass_error']:.3e} |")
    lines.append(f"| ODE | {reference['tip_r']:.9f} | {reference['tip_z']:.9f} | {reference['attachment_r']:.9f} | - | - | - | - |")
    notes='''

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
'''
    (folder/'geometric_ode_convergence_b2.md').write_text(head+'\n'.join(lines)+notes)

    def scientific(value,signed=False):
        mantissa,power=f'{value:{"+" if signed else ""}.2e}'.split('e')
        return rf'${mantissa}\times10^{{{int(power)}}}$'
    tex=[r'\begin{table*}',r'\centering',r'\small',r'\setlength{\tabcolsep}{4pt}',
         r'\begin{tabular}{r r r r r r r r}',r'\hline',
         r'$N_e$ & $r_{\rm tip}$ & $z_{\rm tip}$ & $r_{\rm a}$ & $\Delta r_{\rm tip}$ & $\Delta z_{\rm tip}$ & $\Delta r_{\rm a}$ & $\max_i|\delta m_i|$ \\',r'\hline']
    for row in rows:
        tag=r'$^{\ast}$' if not row['converged'] else ''
        vals=[str(row['n'])+tag,f"{row['tip_r']:.7f}",f"{row['tip_z']:.7f}",f"{row['attachment_r']:.7f}"]
        vals.extend(scientific(row['differences'][k],True) for k in ('tip_r','tip_z','attachment_r'))
        vals.append(scientific(row['max_mass_error']))
        tex.append(' & '.join(vals)+r' \\')
    tex.extend([r'\hline',f"ODE & {reference['tip_r']:.7f} & {reference['tip_z']:.7f} & {reference['attachment_r']:.7f} & --- & --- & --- & ---"+r' \\',
                r'\hline',r'\end{tabular}',
                r'\caption{Front locations for dry geometric-model equilibria with $r_0=0.003$, $r_{\max}=4.2$, $b_0=2$ and $\alpha=1$. Initial elements are uniformly spaced in radius. The tip and upper-lid attachment ($r_{\rm a}$) are obtained from the polygon contact chain. Differences are geometric minus ODE. The absolute element-mass tolerance is $10^{-6}$. Rows marked $\ast$ report final iterates that did not meet this tolerance.}',
                r'\label{tab:geometric-front-convergence}',r'\end{table*}'])
    (folder/'geometric_ode_convergence_b2.tex').write_text('\n'.join(tex)+'\n')
    print(head+'\n'.join(lines))
    print('Saved Markdown, LaTeX fragment and full-precision JSON in',folder.resolve())


if __name__=='__main__':main()
