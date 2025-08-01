import numpy as np
import math
import os
#import pydmga

#VORONOI PY#
from shapely.geometry import Polygon, Point


def mirror_data(nodes, radii, maxLim):
    '''Mirror 2D data'''
    delta = 1E-8
    nodesO = np.copy(nodes)
    radiiO = np.copy(radii)
    dist = np.max(radii)*2.
    ind = np.where(nodesO[:,0]<dist)[0]
    n3 = nodesO[ind,:]
    r3 = radiiO[ind]
    n3[:,0] = - n3[:,0] - delta
    ind = np.where(nodesO[:,0]>maxLim[0]-dist)[0]
    n4 = nodesO[ind,:]
    r4 = radiiO[ind]
    n4[:,0] = 2*maxLim[0] - n4[:,0] + delta
    nodes2 = np.vstack((nodesO,n3, n4))
    radii2 = np.hstack((radiiO,r3, r4))
    ind = np.where(nodesO[:,1]<dist)[0]
    n3 = nodesO[ind,:]
    r3 = radiiO[ind]
    n3[:,1] = - n3[:,1] - delta
    ind = np.where(nodesO[:,1]>maxLim[1]-dist)[0]
    n4 = nodesO[ind,:]
    r4 = radiiO[ind]
    n4[:,1] = 2*maxLim[1] - n4[:,1] + delta
    nodes2 = np.vstack((nodes2,n3, n4))
    radii2 = np.hstack((radii2,r3, r4))
    return nodes2, radii2

def power_2d(nodes, radii, dim, maxLim):

    mirnodes, mirradii = mirror_data(nodes, radii, maxLim)

    G = pydmga.geometry.BoxGeometry(-maxLim[0], -maxLim[1], -1., 2*maxLim[0], 2*maxLim[1], 1.)
    C = pydmga.container.Container(G)
    data = [ (i, mirnodes[i,0], mirnodes[i,1], 0, mirradii[i]) for i in range(len(mirradii))]
    C.add(data)
    D = pydmga.diagram.Diagram(C, False)

    polygons = []
    new_regions = []
    areas = []
    centroids = []
    points = []

    polArr = np.array([(0.,0), (0,maxLim[1]), (maxLim[0],maxLim[1]), (maxLim[0],0)])
    pol=Polygon(polArr)

    S = D.size()
    idlist=np.zeros(0)
    for i in range(S):
        print(S,i)
        C = D.get_cell(i)
        print(i)

        for side in C.sides:
            print(side)
            if i<side.neighbour: idsX = i*10000+side.neighbour
            else: idsX = i+side.neighbour*10000
            if not idsX in idlist and side.neighbour>0 and len(side.vertices)==4:
                idlist = np.hstack((idlist, idsX))
                coord = np.array(side.as_coords())
                xxx = np.where(coord[:,2]>0.)[0]
                #beams.append(Beam(coord[xxx[0]][:2],coord[xxx[1]][:2],i, side.neighbour, "regular",nodes, xlim, ylim))
    exit(1)

    return beams




    polArr = np.array([(0.,0), (0,maxLim[1]), (maxLim[0],maxLim[1]), (maxLim[0],0)])
    polArr += shifts

    pol=Polygon(polArr)
  #  pol=Polygon([(0,0), (0,ySize), (xSize,ySize), (xSize,0)])
    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]
        if all(v >= 0 for v in vertices):
            p = Polygon(vor.vertices[vertices])
            if p.is_empty:
                continue
            if p.geom_type == 'Point' or p.geom_type == 'LineString' or p.area < 1e-10:
                continue
            if not pol.intersects(p) or not pol.intersects(Point(vor.points[p1])):
                continue
            new_regions.append(vertices)
            areas.append(p.area)
            polygons.append(p)
            centroids.append(np.asarray(p.centroid))
            points.append(vor.points[p1])
    if np.sum(areas) < .99:
        print('Area errror', np.sum(areas))

    return new_regions, np.asarray(new_vertices), polygons, np.asarray(areas), np.asarray(centroids), np.asarray(points)
