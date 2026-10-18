"""Q-surface reconstruction for the dry Ou ODE solution (standard library).

Every surface is evaluated analytically below its lid/front termination point.
Monotone cubic interpolation follows the MATLAB 'pchip' convention. The true
tip and the three warmup points are included, filling the tiny launch interval.
No geometric-model polygons enter this reconstruction.
"""
from bisect import bisect_right
from math import exp, sqrt, log
from ou_front_solver import b_Q, lid_radius


class Pchip:
    """Shape-preserving piecewise cubic Hermite interpolation, no extrapolation."""
    def __init__(self, x, y):
        self.x, self.y = list(x), list(y)
        if len(x) < 2 or len(x) != len(y) or any(b <= a for a, b in zip(x, x[1:])):
            raise ValueError('Interpolation coordinates must strictly increase')
        h = [b-a for a,b in zip(x,x[1:])]
        d = [(b-a)/step for a,b,step in zip(y,y[1:],h)]
        if len(x) == 2:
            self.m = [d[0],d[0]]
            return
        def endpoint(h0,h1,d0,d1):
            m = ((2*h0+h1)*d0-h0*d1)/(h0+h1)
            if m*d0 <= 0:
                return 0.
            if d0*d1 < 0 and abs(m)>3*abs(d0):
                return 3*d0
            return m
        slopes=[endpoint(h[0],h[1],d[0],d[1])]
        for i in range(1,len(x)-1):
            if d[i-1]*d[i] <= 0:
                slopes.append(0.)
            else:
                w1,w2=2*h[i]+h[i-1],h[i]+2*h[i-1]
                slopes.append((w1+w2)/(w1/d[i-1]+w2/d[i]))
        slopes.append(endpoint(h[-1],h[-2],d[-1],d[-2]))
        self.m=slopes

    def __call__(self,x):
        if not self.x[0]-1e-12 <= x <= self.x[-1]+1e-12:
            raise ValueError('Interpolation requested outside data range')
        x=max(self.x[0],min(self.x[-1],x))
        i=min(len(self.x)-2,max(0,bisect_right(self.x,x)-1))
        h=self.x[i+1]-self.x[i]
        t=(x-self.x[i])/h
        return ((2*t**3-3*t**2+1)*self.y[i]+(t**3-2*t**2+t)*h*self.m[i]
                +(-2*t**3+3*t**2)*self.y[i+1]+(t**3-t**2)*h*self.m[i+1])


class OuFields:
    def __init__(self,sol):
        self.sol=sol
        self.alpha,self.b0=sol['alpha'],sol['b0']
        self.Qt=sol['Qt']
        self.qm,self.qp=sol['Qm'][-1],sol['Qp'][-1]
        rows=[(sol['Zt'],self.Qt,self.Qt,sol['rf_t'])]
        rows.extend(tuple(row) for row in sol['warmup'][:-1])
        rows.extend(zip(sol['z'],sol['Qm'],sol['Qp'],sol['rf']))
        if any(b[1]>=a[1] or b[2]<=a[2] for a,b in zip(rows,rows[1:])):
            raise ValueError('Front Q branches are not monotone')
        self.front_rows=rows
        self.z_m=Pchip([row[1] for row in reversed(rows)],[row[0] for row in reversed(rows)])
        self.z_p=Pchip([row[2] for row in rows],[row[0] for row in rows])
        self.front_r=Pchip([row[0] for row in rows],[row[3] for row in rows])
        self.q_m=Pchip([row[0] for row in rows],[row[1] for row in rows])
        self.q_p=Pchip([row[0] for row in rows],[row[2] for row in rows])

    def top(self,Q):
        if Q <= 0:
            return 1.,0.
        if Q < self.qm or Q > self.qp:
            return 1.,lid_radius(Q,self.alpha,self.b0)
        z=self.z_m(Q) if Q <= self.Qt else self.z_p(Q)
        return z,self.front_r(z)

    def radius(self,Q,z):
        zt,rt=self.top(Q)
        if z > zt+1e-12:
            raise ValueError('Q surface terminates below requested height')
        if Q<=0:
            return 0.
        return 1/sqrt(1/rt**2+2*b_Q(Q,self.alpha,self.b0)*(z-zt))

    def curve(self,Q,samples=101):
        zt,rt=self.top(Q)
        if Q<=0:
            return [(0.,zt*i/(samples-1)) for i in range(samples)]
        rb=self.radius(Q,0.)
        points=[]
        for i in range(samples):
            r=rb+(rt-rb)*i/(samples-1)
            z=zt+(1/r**2-1/rt**2)/(2*b_Q(Q,self.alpha,self.b0))
            points.append((r,max(0.,min(zt,z))))
        return points

    def strip(self,qa,qb,samples=81):
        """Polygon between two Q surfaces, closed along lid/front and floor."""
        za,ra=self.top(qa)
        zb,rb=self.top(qb)
        points=self.curve(qa,samples)
        if self.qm <= qa <= self.qp and self.qm <= qb <= self.qp:
            # Sample the connecting front itself, not a straight physical line.
            nfront=max(31,samples//2)
            points.extend((self.front_r(za+(zb-za)*i/nfront),za+(zb-za)*i/nfront)
                          for i in range(1,nfront+1))
        else:
            points.append((rb,zb))
        points.extend(reversed(self.curve(qb,samples)))
        return points

    def labels(self,rmax=4.2,count=900):
        """Dense shading labels: uniform b in core, logarithmic Q outside."""
        qmax=rmax**4/4
        # Exponential b: uniform-in-b samples resolve its visible gradient.
        labels=[0.,self.qm,self.Qt,self.qp,qmax]
        labels.extend(-log(i/count)/self.alpha for i in range(1,count+1))
        labels.extend(exp(log(1e-10)+(log(qmax)-log(1e-10))*i/count)
                      for i in range(count+1))
        return sorted(set(q for q in labels if 0 <= q <= qmax))

    def buoyancy_contours(self,rmin=.2,rmax=3.2):
        levels=[.4*i for i in range(1,int(self.b0/.4+1e-10)+1) if .4*i < self.b0-1e-12]
        paths=[]
        for b in levels:
            Q=-log(b/self.b0)/self.alpha
            curve=self.curve(Q,401)
            keep=[point for point in curve if rmin<=point[0]<=rmax]
            if len(keep)>1:
                paths.append((b,keep))
        return paths

    def grid(self,nr=1600,nz=1200,nq=8000,rmax=4.2):
        """Invert the four monotone Q families onto a physical plotting grid.

        NumPy is used only for this plotting array calculation. The two frontal
        branches are inverted separately so interpolation cannot span the jump.
        """
        import numpy as np
        def interpolate(x,y,newx):
            p=Pchip(x.tolist(),y.tolist())
            slopes=np.asarray(p.m)
            i=np.clip(np.searchsorted(x,newx,side='right')-1,0,len(x)-2)
            h=x[i+1]-x[i];t=(newx-x[i])/h
            return ((2*t**3-3*t**2+1)*y[i]+(t**3-2*t**2+t)*h*slopes[i]
                    +(-2*t**3+3*t**2)*y[i+1]+(t**3-t**2)*h*slopes[i+1])
        if min(nr,nz,nq)<3:
            raise ValueError('Use at least three samples per coordinate/family')
        rgrid=np.linspace(0,rmax,nr);zgrid=np.linspace(0,1,nz)
        Qgrid=np.full((nz,nr),np.nan)
        bounds=(0.,self.qm,self.Qt,self.qp,rmax**4/4)
        families=[]
        for lo,hi in zip(bounds,bounds[1:]):
            q=np.linspace(lo,hi,nq)
            if lo==0:
                q=np.unique(np.concatenate((q,np.geomspace(1e-14,hi,nq))))
            tops=np.asarray([self.top(float(v)) for v in q])
            families.append((q,tops[:,0],tops[:,1]))
        for row,z in enumerate(zgrid):
            for family,(q,zt,rt) in enumerate(families):
                mask=zt>=z-1e-12
                qs=q[mask];zs=zt[mask];rs=rt[mask]
                # Add the exact branch endpoint at this height. No extrapolation
                # across the missing interval of Q labels above the front tip.
                if z>self.sol['Zt'] and family in (1,2):
                    endq=self.q_m(float(z)) if family==1 else self.q_p(float(z))
                    qs=np.append(qs,endq);zs=np.append(zs,z);rs=np.append(rs,self.front_r(float(z)))
                if len(qs)<2:continue
                with np.errstate(divide='ignore',invalid='ignore'):
                    radii=1/np.sqrt(1/rs**2+2*(-self.alpha*self.b0*np.exp(-self.alpha*qs))*(z-zs))
                radii[qs==0]=0.
                ix=np.argsort(radii);radii=radii[ix];qs=qs[ix]
                keep=np.concatenate(([True],np.diff(radii)>1e-13))
                radii=radii[keep];qs=qs[keep]
                if len(qs)<2:continue
                selection=(rgrid>=radii[0]-1e-12)&(rgrid<=radii[-1]+1e-12)
                Qgrid[row,selection]=interpolate(radii,qs,np.clip(rgrid[selection],radii[0],radii[-1]))
            Qgrid[row,0]=0.
        if not np.isfinite(Qgrid).all():
            raise ArithmeticError(f'Field reconstruction left {np.sum(~np.isfinite(Qgrid))} unassigned grid points')
        return rgrid,zgrid,Qgrid

    def velocity_contours(self,levels,rmin=.2,rmax=3.2,count=6000):
        """Find v contours parametrically in Q; never join across a front.

        v=sqrt(Q)/r-r/2 implies r=-v+sqrt(v^2+2*sqrt(Q)).
        Substitute this r into the analytic Q-surface equation to obtain z.
        """
        bounds=(1e-10,self.qm,self.Qt,self.qp,4.2**4/4)
        paths=[]
        for v in levels:
            for lo,hi in zip(bounds,bounds[1:]):
                path=[]
                for i in range(count+1):
                    Q=exp(log(lo)+(log(hi)-log(lo))*i/count)
                    r=-v+sqrt(v*v+2*sqrt(Q))
                    zt,rt=self.top(Q)
                    z=zt+(1/r**2-1/rt**2)/(2*b_Q(Q,self.alpha,self.b0))
                    if rmin <= r <= rmax and -1e-9 <= z <= zt+1e-9:
                        path.append((r,max(0.,min(zt,z))))
                    else:
                        if len(path)>1:
                            paths.append((v,path))
                        path=[]
                if len(path)>1:
                    paths.append((v,path))
        return paths
