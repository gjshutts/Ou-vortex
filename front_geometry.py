"""Extract nonlocal contacts in the ordered-label zero-PV discretisation.

The ordinary family interfaces join consecutive Q labels. Interfaces joining
nonconsecutive labels are candidates for the discrete frontal branch. This
criterion is specific to this one-dimensional, ordered zero-PV family; it is
not a generic front detector for an ocean lens or arbitrary PV distribution.
"""
import math
from coordinate_maps import r_from_R


def frontal_edges(result, tolerance=1e-9):
    planes=result['planes'];polygons=result['polygons']
    if not all(a[0]<b[0] for a,b in zip(planes,planes[1:])):
        raise ValueError('This detector requires strictly increasing Q labels.')
    edges=[]
    for i,polygon in enumerate(polygons):
        for k,start in enumerate(polygon):
            end=polygon[(k+1)%len(polygon)]
            for j in range(i+2,len(planes)):
                if not polygons[j]:continue
                a,b,c=[planes[i][v]-planes[j][v] for v in range(3)]
                normal=math.hypot(a,b)
                if normal and max(abs(math.fsum((a*p[0],b*p[1],c)))/normal for p in (start,end))<=tolerance:
                    edges.append({'i':i,'j':j,'start':start,'end':end,
                                  'buoyancy_jump':abs(planes[i][1]-planes[j][1])})
    return edges


def front_summary(result):
    edges=frontal_edges(result)
    if not edges:return {'segments':0}
    r0=result['parameters']['r_inner']
    points=[p for e in edges for p in (e['start'],e['end'])]
    tip=min(points,key=lambda p:p[1]);roof=max(points,key=lambda p:p[1])
    return {'segments':len(edges),'tip_r':r_from_R(tip[0],r0),'tip_z':tip[1],
            'upper_r':r_from_R(roof[0],r0),'upper_z':roof[1],
            'minimum_label_gap':min(e['j']-e['i'] for e in edges),
            'maximum_label_gap':max(e['j']-e['i'] for e in edges)}


def check_front_chain(edges, tolerance=1e-8):
    """Require a single connected polyline, with two endpoints and no branches."""
    nodes=[];links=[]
    for edge in edges:
        pair=[]
        for point in (edge['start'],edge['end']):
            match=next((i for i,p in enumerate(nodes) if math.dist(p,point)<tolerance),None)
            if match is None:
                match=len(nodes);nodes.append(point)
            pair.append(match)
        if pair[0]!=pair[1]:links.append(pair)
    neighbours=[[] for _ in nodes]
    for a,b in links:
        neighbours[a].append(b);neighbours[b].append(a)
    visited={0};pending=[0]
    while pending:
        for j in neighbours[pending.pop()]:
            if j not in visited:visited.add(j);pending.append(j)
    degrees=[len(v) for v in neighbours]
    ok=len(visited)==len(nodes) and degrees.count(1)==2 and all(d in (1,2) for d in degrees)
    if not ok:raise ValueError('Candidate front is not one unbranched connected chain.')
    return {'nodes':len(nodes),'segments':len(links),'endpoints':2,'connected':True}
