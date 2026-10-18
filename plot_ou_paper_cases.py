"""Reproduce frontal paper cases from the Python ODE solver, not elements.

python3 plot_ou_paper_cases.py             # SVG outputs (NumPy and Pillow)
python3 plot_ou_paper_cases.py --pdf       # also vector PDFs (ReportLab)
"""
import argparse
import json
import base64
from io import BytesIO
import math
from pathlib import Path
from time import perf_counter
from html import escape

from ou_front_solver import solve_zero_pv_vortex_front
from ou_fields import OuFields


def rgb(u,velocity=False):
    u=max(-1.,min(1.,u)) if velocity else max(0.,min(1.,u))
    end=(33/255,102/255,172/255) if u<0 else ((178/255,24/255,43/255) if velocity else (139/255,0.,0.))
    return tuple(1+abs(u)*(c-1) for c in end)


def hexcolour(c):
    return '#'+''.join(f'{round(255*v):02x}' for v in c)


class VectorFigure:
    """A small shared SVG/PDF drawing interface; all plotted geometry is vector."""
    def __init__(self,path,width,height,pdf):
        self.path,self.width,self.height=path,width,height
        self.svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                  '<rect width="100%" height="100%" fill="white"/>']
        self.ids=0
        self.canvas=None
        if pdf:
            from reportlab.pdfgen.canvas import Canvas
            self.canvas=Canvas(str(path.with_suffix('.pdf')),pagesize=(width,height),pageCompression=1)
            self.canvas.setTitle('Python ODE: dry rotational Ou paper cases')
            self.canvas.setAuthor('Glenn Shutts')

    def text(self,x,y,label,size=14,anchor='middle'):
        self.svg.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" text-anchor="{anchor}" fill="#172b40">{escape(label)}</text>')
        if self.canvas:
            c=self.canvas;c.setFillColorRGB(.09,.17,.25);c.setFont('Helvetica',size)
            fn=c.drawCentredString if anchor=='middle' else c.drawRightString if anchor=='end' else c.drawString
            fn(x,self.height-y,label)

    def _path(self,points,closed=False):
        p=self.canvas.beginPath();p.moveTo(points[0][0],self.height-points[0][1])
        for x,y in points[1:]:p.lineTo(x,self.height-y)
        if closed:p.close()
        return p

    def line(self,points,width=1.,colour=(0.,0.,0.),closed=False):
        if len(points)<2:return
        coords=' '.join(f'{x:.4f},{y:.4f}' for x,y in points)
        self.svg.append(f'<polyline points="{coords}" fill="none" stroke="{hexcolour(colour)}" stroke-width="{width}" stroke-linejoin="round"/>')
        if self.canvas:
            self.canvas.setStrokeColorRGB(*colour);self.canvas.setLineWidth(width)
            self.canvas.drawPath(self._path(points,closed),stroke=1,fill=0)

    def polygon(self,points,colour,gradient=None,bleed=0.):
        coords=' '.join(f'{x:.4f},{y:.4f}' for x,y in points)
        fill=hexcolour(colour)
        if gradient:
            x0,x1,colours=gradient
            self.ids+=1;ident=f'g{self.ids}'
            stops=''.join(f'<stop offset="{i/(len(colours)-1):.6f}" stop-color="{hexcolour(c)}"/>' for i,c in enumerate(colours))
            self.svg.append(f'<defs><linearGradient id="{ident}" gradientUnits="userSpaceOnUse" x1="{x0}" x2="{x1}" y1="0" y2="0">{stops}</linearGradient></defs>')
            fill=f'url(#{ident})'
        self.svg.append(f'<polygon points="{coords}" fill="{fill}" stroke="{fill}" stroke-width="{max(.1,2*bleed)}" stroke-linejoin="round"/>')
        if self.canvas:
            c=self.canvas;p=self._path(points,True)
            if gradient:
                if bleed>0:
                    # Subpixel overlap avoids white antialias seams between
                    # independently clipped gradient fills. Geometry/linework
                    # still uses the original computed boundary coordinates.
                    area=sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(points,points[1:]+points[:1]))
                    sign=1 if area>0 else -1
                    normals=[]
                    for a,b in zip(points,points[1:]+points[:1]):
                        dx,dy=b[0]-a[0],b[1]-a[1];length=math.hypot(dx,dy)
                        normals.append((sign*dy/length,-sign*dx/length) if length else (0.,0.))
                    expanded=[]
                    for i,(x,y) in enumerate(points):
                        nx=normals[i-1][0]+normals[i][0];ny=normals[i-1][1]+normals[i][1]
                        dot=nx*normals[i][0]+ny*normals[i][1]
                        scale=bleed/max(dot,.5)
                        expanded.append((x+scale*nx,y+scale*ny))
                    p=self._path(expanded,True)
                c.saveState();c.clipPath(p,stroke=0,fill=0)
                c.linearGradient(x0,0,x1,0,list(colours),extend=True)
                c.restoreState()
            else:
                c.setFillColorRGB(*colour);c.setStrokeColorRGB(*colour);c.setLineWidth(.10)
                c.drawPath(p,stroke=1,fill=1)

    def clip(self,x,y,w,h):
        self.ids+=1;ident=f'clip{self.ids}'
        self.svg.append(f'<defs><clipPath id="{ident}"><rect x="{x}" y="{y}" width="{w}" height="{h}"/></clipPath></defs><g clip-path="url(#{ident})">')
        if self.canvas:
            self.canvas.saveState();p=self.canvas.beginPath();p.rect(x,self.height-y-h,w,h)
            self.canvas.clipPath(p,stroke=0,fill=0)

    def endclip(self):
        self.svg.append('</g>')
        if self.canvas:self.canvas.restoreState()

    def image(self,x,y,w,h,png):
        encoded=base64.b64encode(png).decode('ascii')
        self.svg.append(f'<image x="{x}" y="{y}" width="{w}" height="{h}" preserveAspectRatio="none" href="data:image/png;base64,{encoded}"/>')
        if self.canvas:
            from reportlab.lib.utils import ImageReader
            self.canvas.drawImage(ImageReader(BytesIO(png)),x,self.height-y-h,w,h)

    def save(self):
        self.path.write_text('\n'.join(self.svg+['</svg>']))
        if self.canvas:self.canvas.showPage();self.canvas.save()


def prepare(b0,bracket,nr=1600,nz=1200,nq=8000):
    import numpy as np
    from PIL import Image
    start=perf_counter()
    shot,sol=solve_zero_pv_vortex_front(1.,b0,*bracket,verbose=False,
                                      rtol=1e-10,atol=1e-12,tip_tolerance=1e-11)
    solve_time=perf_counter()-start
    field=OuFields(sol)
    r,z,Q=field.grid(nr,nz,nq)
    B=b0*np.exp(-Q)
    V=np.zeros_like(Q)
    V[:,1:]=np.sqrt(Q[:,1:])/r[None,1:]-r[None,1:]/2
    # Retain only the regular-axis value V[:,0]=0. No finite-radius patch.
    vmax=float(np.max(np.abs(V)))
    def png(array,velocity=False):
        u=np.clip(array,-1,1) if velocity else np.clip(array,0,1)
        end=np.where(u[:,:,None]<0,np.array([33,102,172]),np.array([178,24,43])) if velocity else np.array([139,0,0])
        colours=np.rint(255+np.abs(u[:,:,None])*(end-255)).astype(np.uint8)
        out=BytesIO();Image.fromarray(colours[::-1]).save(out,format='PNG')
        return out.getvalue()
    dv=.1 if vmax<=1 else .2
    levels=[dv*i for i in range(-int(vmax/dv),int(vmax/dv)+1) if i]
    return {'b0':b0,'shot':shot,'sol':sol,'fields':field,'vmax':vmax,
            'nr':nr,'nz':nz,
            'b_image':png(B/b0),'v_image':png(V/vmax,True),
            'field_min_Q':float(Q.min()),'field_max_Q':float(Q.max()),
            'inner_velocity_min':float(V[:,r<=.2].min()),
            'inner_velocity_max':float(V[:,r<=.2].max()),
            'vmax_with_axis_patch':float(np.max(np.abs(V[:,r>.2]))),
            'b_contours':field.buoyancy_contours(),
            'v_contours':field.velocity_contours(levels),'velocity_levels':levels,
            'solve_seconds':solve_time}


def panel(fig,data,left,top,width,depth,velocity=False):
    b0=data['b0'];scale=data['vmax'] if velocity else b0
    xy=lambda r,z:(left+r/4.2*width,top+(1-z)*depth)
    fig.text(left+width/2,top-20,f"{'Tangential velocity v' if velocity else 'Buoyancy b'}; b0 = {b0:g}",19)
    fig.clip(left,top,width,depth)
    # Grid values are pixel centres, as in MATLAB imagesc. Clip their half-pixel
    # extensions at the physical domain boundary rather than shifting centres.
    dx,dy=width/(data['nr']-1),depth/(data['nz']-1)
    fig.image(left-dx/2,top-dy/2,width+dx,depth+dy,
              data['v_image'] if velocity else data['b_image'])
    paths=data['v_contours'] if velocity else data['b_contours']
    for value,path in paths:fig.line([xy(r,z) for r,z in path],width=.55)
    field=data['fields']
    front=[xy(field.front_r(solz),solz) for solz in
           [field.sol['Zt']+(1-field.sol['Zt'])*i/600 for i in range(601)]]
    fig.line(front,width=1.7)
    fig.endclip()
    fig.line([xy(0,0),xy(4.2,0),xy(4.2,1),xy(0,1),xy(0,0)],width=1.)
    for r in (0,1,2,3,4,4.2):
        x,y=xy(r,0);fig.line([(x,y),(x,y+5)])
        fig.text(x,y+22,f'{r:g}',13)
    for z in (0,.2,.4,.6,.8,1):
        x,y=xy(0,z);fig.line([(x-5,y),(x,y)])
        fig.text(x-10,y+4,f'{z:g}',13,'end')
    fig.text(left+width/2,top+depth+45,'r',18)
    fig.text(left-40,top+depth/2,'z',18)
    barx=left+width+19
    for i in range(120):
        u=1-(i+.5)/120
        y=top+depth*i/120
        fig.polygon([(barx,y),(barx+13,y),(barx+13,y+depth/120+.1),(barx,y+depth/120+.1)],
                    rgb(2*u-1 if velocity else u,velocity))
    for u in (0,.25,.5,.75,1):
        value=(2*u-1)*scale if velocity else u*scale
        fig.text(barx+20,top+(1-u)*depth+4,f'{value:.2g}',12,'start')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf',action='store_true',help='Also create PDFs using ReportLab')
    parser.add_argument('--nr',type=int,default=1600)
    parser.add_argument('--nz',type=int,default=1200)
    parser.add_argument('--nq',type=int,default=8000)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'output/pdf')
    args=parser.parse_args()
    if min(args.nr,args.nz,args.nq)<100:parser.error('Use at least 100 samples per coordinate/family')
    args.output.mkdir(parents=True,exist_ok=True)
    data=[];summary=[]
    for b,bracket in ((2.,(.75,.85)),(4.,(.70,.75)),(7.,(.68,.72))):
        item=prepare(b,bracket,args.nr,args.nz,args.nq);data.append(item)
        summary.append({'b0':b,'alpha':1.,'plot_rmax':4.2,**item['shot'],
                        'tip_radius':item['sol']['rf_t'],'vmax':item['vmax'],
                        'velocity_contours':item['velocity_levels'],'solve_seconds':item['solve_seconds']})
        summary[-1].update(nr=args.nr,nz=args.nz,nq=args.nq,
                           field_min_Q=item['field_min_Q'],field_max_Q=item['field_max_Q'],
                           inner_velocity_min=item['inner_velocity_min'],
                           inner_velocity_max=item['inner_velocity_max'],
                           vmax_with_axis_patch=item['vmax_with_axis_patch'])
        Path(__file__).with_name(f'ou_front_b{b:g}.json').write_text(json.dumps(
            {'shot':item['shot'],'sol':item['sol'],'elapsed_seconds':item['solve_seconds']},indent=2,allow_nan=False)+'\n')
        print(f"ODE b0={b:g}: tip z={item['shot']['Zt']:.9f}; lid r={item['shot']['front_lid_radius']:.9f}",flush=True)
    overview=VectorFigure(args.output/'Ou_ODE_Paper_Cases_b2_b4_b7.svg',1380,1120,args.pdf)
    overview.text(690,30,'Python ODE solver: dry rotational Ou equilibria',24)
    overview.text(690,54,'Physical r and z; black front and analytical field contours; plotting window r = 0 to 4.2',15)
    for row,item in enumerate(data):
        panel(overview,item,65,110+row*330,525,235)
        panel(overview,item,750,110+row*330,525,235,True)
        fig=VectorFigure(args.output/f'Ou_ODE_b{item["b0"]:g}_buoyancy_velocity.svg',1380,520,args.pdf)
        fig.text(690,30,f'Python ODE solver: b0 = {item["b0"]:g}, alpha = 1; dry, incompressible',24)
        panel(fig,item,65,100,525,310)
        panel(fig,item,750,100,525,310,True)
        fig.text(690,495,'Front: heavy black line. Velocity evaluated for every r > 0; v = 0 only at the axis.',14)
        fig.save()
    overview.text(690,1090,'Velocity evaluated for every r > 0; v = 0 only at the axis. These are ODE results, not geometric elements.',14)
    overview.save()
    (args.output/'Ou_ODE_Paper_Cases_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('Saved vector figures to',args.output)


if __name__=='__main__':main()
