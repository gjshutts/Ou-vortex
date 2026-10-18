"""Vector velocity or buoyancy shading of a saved dry geometric run.

Uses v=sqrt(Q_i)/r-r/2 inside each computed element; no ODE field is inserted.
SVG needs only the standard library; optional PDF export needs ReportLab.
"""
import argparse
import json
import math
from pathlib import Path

from coordinate_maps import r_from_R, R_from_r
from front_geometry import frontal_edges,check_front_chain
from plot_ou_paper_cases import VectorFigure,rgb


def sample_edge(a,b,r0):
    ra,rb=r_from_R(a[0],r0),r_from_R(b[0],r0)
    count=max(3,math.ceil(abs(rb-ra)/.001))
    points=[]
    for i in range(count+1):
        r=ra+(rb-ra)*i/count
        t=i/count if a[0]==b[0] else (R_from_r(r,r0)-a[0])/(b[0]-a[0])
        points.append((r,a[1]+t*(b[1]-a[1])))
    return points


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result',type=Path)
    parser.add_argument('--allow-unconverged',action='store_true')
    parser.add_argument('--pdf',action='store_true')
    parser.add_argument('--field',choices=('v','b'),default='v')
    parser.add_argument('--edge-width',type=float,default=.075)
    parser.add_argument('--vlim',type=float,default=1.05)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'output/pdf/Geometric_b2_N400_r0_0.003_rmax4.2_velocity.svg')
    args=parser.parse_args()
    result=json.loads(args.result.read_text());p=result['parameters']
    if not result['converged'] and not args.allow_unconverged:
        raise ValueError('Result has not met its mass tolerance; pass --allow-unconverged explicitly to plot it.')
    if p['H0']!='Infinity' or p['coefficients'][1:]!=[0,0] or result['empty_cells']:
        raise ValueError('Expected a dry constant-density result with no empty cells.')
    if args.vlim<=0 or args.edge_width<0:raise ValueError('Invalid colour limit or edge width')
    r0,rmax=p['r_inner'],p['r_outer']
    physical=[]
    for poly in result['polygons']:
        points=[]
        for a,b in zip(poly,poly[1:]+poly[:1]):points.extend(sample_edge(a,b,r0)[:-1])
        physical.append(points)
    velocity=args.field=='v'
    extrema=([math.sqrt(plane[0])/r-r/2 for plane,poly in zip(result['planes'],physical) for r,z in poly]
             if velocity else [plane[1] for plane in result['planes']])
    if velocity and max(abs(v) for v in extrema)>args.vlim:raise ValueError('Colour limit would clip computed velocities')
    edges=frontal_edges(result);chain=check_front_chain(edges)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig=VectorFigure(args.output,1000,650,args.pdf)
    left,top,width=75,85,790
    depth=width*310/525  # Same field aspect ratio as the side-by-side ODE plots.
    xy=lambda r,z:(left+r/rmax*width,top+(1-z)*depth)
    field_name='tangential velocity' if velocity else 'buoyancy'
    fig.text(500,29,f'Geometric model: {field_name}, {p["n"]} elements',23)
    fig.text(500,54,f'b0 = {p["bm0"]:g}, alpha = {p["alpha"]:g}, r0 = {r0:g}, rmax = {rmax:g}; dry, constant density',16)
    fig.clip(left,top,width,depth)
    for plane,poly in zip(result['planes'],physical):
        if not velocity:
            fig.polygon([xy(r,z) for r,z in poly],rgb(plane[1]/p['bm0']))
            continue
        Q=plane[0];lo=min(r for r,z in poly);hi=max(r for r,z in poly)
        if hi<=lo:continue
        colours=[rgb((math.sqrt(Q)/(lo+(hi-lo)*i/128)-(lo+(hi-lo)*i/128)/2)/args.vlim,True) for i in range(129)]
        fig.polygon([xy(r,z) for r,z in poly],colours[64],(xy(lo,0)[0],xy(hi,0)[0],colours),bleed=.2)
    # Draw boundaries after all fills; light, very thin lines reveal the cells.
    if args.edge_width>0:  # PDF width zero means a hairline, not an omitted edge.
        for poly in physical:
            fig.line([xy(r,z) for r,z in poly+[poly[0]]],width=args.edge_width,colour=(.68,.68,.68))
    for edge in edges:
        fig.line([xy(r,z) for r,z in sample_edge(edge['start'],edge['end'],r0)],width=1.4)
    fig.endclip()
    fig.line([xy(0,0),xy(rmax,0),xy(rmax,1),xy(0,1),xy(0,0)],width=1.)
    for r in (0,1,2,3,4,4.2):
        x,y=xy(r,0);fig.line([(x,y),(x,y+5)])
        fig.text(x,y+24,f'{r:g}',15)
    for z in (0,.2,.4,.6,.8,1):
        x,y=xy(0,z);fig.line([(x-5,y),(x,y)])
        fig.text(x-12,y+5,f'{z:g}',15,'end')
    fig.text(left+width/2,top+depth+48,'r',19)
    fig.text(27,top+depth/2,'z',19)
    fig.text(903,top-12,args.field,18)
    for i in range(256):
        u=1-(i+.5)/256;y=top+depth*i/256
        fig.polygon([(890,y),(909,y),(909,y+depth/256+.1),(890,y+depth/256+.1)],rgb(2*u-1 if velocity else u,velocity))
    if velocity:
        for v in (-args.vlim,-.5,0.,.5,args.vlim):
            fig.text(920,top+depth*(1-v/args.vlim)/2+5,f'{v:g}',15,'start')
    else:
        for u in (0,.25,.5,.75,1):
            fig.text(920,top+depth*(1-u)+5,f'{p["bm0"]*u:g}',15,'start')
    edge_label='thin element edges' if args.edge_width>0 else 'element outlines omitted'
    fig.text(500,631,f'Physical r and z; {edge_label}; black front. Maximum mass residual {result["max_mass_error"]:.2e}.',14)
    fig.save()
    metadata={'source':str(args.result.resolve()),'converged':result['converged'],
              'max_mass_error':result['max_mass_error'],args.field+'min':min(extrema),args.field+'max':max(extrema),
              'colour_limits':[-args.vlim,args.vlim] if velocity else [0,p['bm0']],'element_edge_width':args.edge_width,
              'front_chain':chain,'shading':'Within each element, v=sqrt(Q_i)/r-r/2; SVG/PDF vector gradients' if velocity else 'Constant computed buoyancy label within each element; white to dark red',
              'axis_patch':False,'field_aspect_ratio':width/depth}
    args.output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2));print('Saved',args.output.resolve())


if __name__=='__main__':main()
