import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import os

from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay

#VORONOI PY#
from shapely.geometry import Polygon, Point


def mirror_data(data):
    '''Cast data array to 3x3 array by mirroring input data'''
    return np.vstack((np.array([0,2]) + np.array([-1,-1]) * data, np.array([0,2]) + data * np.array([1,-1]), np.array([2,2]) + data * np.array([-1,-1]),
                  np.array([-1,1]) * data, data, np.array([2,0]) + data * np.array([-1,1]),
                  data * np.array([-1,-1]), np.array([1,-1]) * data, np.array([2,0]) + data * np.array([-1,-1])))

def copy_data(data):
    '''Cast data array to 3x3 array'''
    return np.vstack((data + np.array([-1, 1]), data + np.array([0,1]), data + np.array([1,1]),
               data + np.array([-1,0]), data, data + np.array([1,0]),
               data + np.array([-1,-1]), data + np.array([0,-1]), data + np.array([1,-1])))

def mirror_dataBeam(data, dim, sizes, shifts=0):
    if (dim == 2):
        dataOut= np.vstack((
        data,
        np.array([0,0]) + data * np.array([-1,1]),
        np.array([sizes[0]*2,0]) + data * np.array([-1,1]),
        np.array([0,sizes[1]*2]) + data * np.array([1,-1]),
        np.array([0,0]) + data * np.array([1,-1])
        ))
    
    if (dim == 3):
        dataOut =  np.vstack((data,
            np.array([0,0,0]) + data * np.array([-1,1,1]),
            np.array([ sizes[0]*2 ,0,0]) + data * np.array([-1,1,1]),
            np.array([ 0 ,0,0]) + data * np.array([1,-1,1]),
            np.array([ 0 ,sizes[1]*2,0]) + data * np.array([1,-1,1]),
            np.array([ 0 ,0,0]) + data * np.array([1,1,-1]),
            np.array([ 0 ,0,sizes[2]*2]) + data * np.array([1,1,-1])
        ))

    dataOut += shifts
    return dataOut



def copy_dataBeam(data, dim, sizes):
    if (dim == 2):
        return np.vstack((
            data + np.array([0, 0]),
            data + np.array([-sizes[0], sizes[1]]),
            data + np.array([0,sizes[1]]),
            data + np.array([sizes[0],sizes[1]]),
            data + np.array([-sizes[0],0]),
            data + np.array([sizes[0],0]),
            data + np.array([-sizes[0],-sizes[1]]),
            data + np.array([0,-sizes[1]]),
            data + np.array([sizes[0],-sizes[1]])
        ))

    if (dim == 3):
        return np.vstack((
            data + np.array([0, 0,0]),
            data + np.array([-sizes[0], sizes[1],0]),
            data + np.array([0,sizes[1],0]),
            data + np.array([sizes[0],sizes[1],0]),
            data + np.array([-sizes[0],0,0]),
            data + np.array([sizes[0],0,0]),
            data + np.array([-sizes[0],-sizes[1],0]),
            data + np.array([0,-sizes[1],0]),
            data + np.array([sizes[0],-sizes[1],0]),
            data + np.array([0, 0,sizes[2]]),
            data + np.array([-sizes[0], sizes[1],sizes[2]]),
            data + np.array([0,sizes[1],sizes[2]]),
            data + np.array([sizes[0],sizes[1],sizes[2]]),
            data + np.array([-sizes[0],0,sizes[2]]),
            data + np.array([sizes[0],0,sizes[2]]),
            data + np.array([-sizes[0],-sizes[1],sizes[2]]),
            data + np.array([0,-sizes[1],sizes[2]]),
            data + np.array([sizes[0],-sizes[1],sizes[2]]),
            data + np.array([0, 0,-sizes[2]]),
            data + np.array([-sizes[0], sizes[1],-sizes[2]]),
            data + np.array([0,sizes[1],-sizes[2]]),
            data + np.array([sizes[0],sizes[1],-sizes[2]]),
            data + np.array([-sizes[0],0,-sizes[2]]),
            data + np.array([sizes[0],0,-sizes[2]]),
            data + np.array([-sizes[0],-sizes[1],-sizes[2]]),
            data + np.array([0,-sizes[1],-sizes[2]]),
            data + np.array([sizes[0],-sizes[1],-sizes[2]])
            ))


def voronoi_2d(vor, sizes):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    polygons = []
    new_regions = []
    areas = []
    centroids = []
    points = []
    new_vertices = vor.vertices.tolist()

    pol=Polygon([(0,0), (0,sizes[1]), (sizes[0],sizes[1]), (sizes[0],0)])
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


def voronoi_3d(vor, sizes):

    if vor.points.shape[1] != 3:
        raise ValueError("Requires 3D input")

    volumes = []
    #
    for i in range (vor.points.shape[0]):
        comp = True
        for p in range (3):
            if (vor.points[i][p] < 0 or vor.points[i][p] > sizes[p] ):
                comp = False

        if (comp):
            vol = cellVolume3d(vor, i)
            volumes.append(vol)
            #print (vol)

   # print (np.sum(volumes))

    return volumes



def cellVolume3d(vor,p):
     #'''Calculate volume of 3d Voronoi cell based on point p. Voronoi diagram is passed in v.'''
    dpoints=[]
    vol=0
    for v in vor.regions[vor.point_region[p]]:
        dpoints.append(list(vor.vertices[v]))
    # ConvexHull .volume
    tri=Delaunay(np.array(dpoints))
    #
    for simplex in tri.simplices:
        vol+=tetravol(np.array(dpoints[simplex[0]]),np.array(dpoints[simplex[1]]),np.array(dpoints[simplex[2]]),np.array(dpoints[simplex[3]]))
    #
    return vol

def tetravol(a,b,c,d):
     #'''Calculates the volume of a tetrahedron, given vertices a,b,c and d (triplets)'''
    tetravol=abs(np.dot((a-d),np.cross((b-d),(c-d))))/6
    return tetravol
