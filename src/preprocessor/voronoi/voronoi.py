import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
import math
import os
import itertools

from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay

#VORONOI PY#
from shapely.geometry import Polygon, Point

SHOW_PLOT = False

def mirror_data(data):
    '''Mirror 2D data'''
    return np.vstack((np.array([0, 2]) + np.array([-1, -1]) * data,
                      np.array([0, 2]) + data * np.array([1, -1]),
                      np.array([2, 2]) + data * np.array([-1, -1]),
                      np.array([-1, 1]) * data,
                      data,
                      np.array([2, 0]) + data * np.array([-1, 1]),
                      data * np.array([-1, -1]), np.array([1, -1]) * data,
                      np.array([2, 0]) + data * np.array([-1, -1])))


def copy_data(data):
    '''Copy 2D data'''
    return np.vstack((data + np.array([-1, 1]),
                      data + np.array([0, 1]),
                      data + np.array([1, 1]),
                      data + np.array([-1, 0]),
                      data,
                      data + np.array([1, 0]),
                      data + np.array([-1, -1]),
                      data + np.array([0, -1]),
                      data + np.array([1, -1])))


def mirror_data_general_full(data):
    '''Mirror data in N dimensions'''
    nvar = data.shape[1]
    new_data = data.copy()
    for d in itertools.product([-1, 0, 1], repeat=nvar):
        if not np.all(np.array(d) == 0):
            mult = np.asanyarray(d)
            mult[mult != 0] = -1
            mult[mult == 0] = 1
            add = np.asanyarray(d)
            add[add != 0] += 1
            new_data = np.vstack((new_data, mult * data + add))
    return new_data


def mirror_data_general_perp(data):
    '''Mirror data in one axis direction'''
    nvar = data.shape[1]
    new_data = data.copy()
    vectors = []
    for i in [-1, 1]:
        row = [i] + (nvar - 1) * [0]
        vectors.extend(itertools.permutations(row))
    for d in vectors:
        if not np.all(np.array(d) == 0):
            mult = np.asanyarray(d)
            mult[mult != 0] = -1
            mult[mult == 0] = 1
            add = np.asanyarray(d)
            add[add != 0] += 1
            new_data = np.vstack((new_data, mult * data + add))
    return new_data


def copy_data_general_full(data):
    '''Copy data'''
    nvar = data.shape[1]
    new_data = data.copy()
    for d in itertools.product([-1, 0, 1], repeat=nvar):
        if not np.all(d == 0):
            new_data = np.vstack((new_data, data + np.asanyarray(d)))
    return new_data


def mirror_dataBeam(data, dim, sizes, shifts=0, weights=None):
    '''Mirror data 2D and 3D'''
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
    if weights is not None:
        if dim == 2:
            weightsOut = np.hstack([weights]*5)
        else:
            weightsOut = np.hstack([weights]*7)
        return dataOut, weightsOut
    return dataOut

def mirror_dataDogBone(data, dim, D, thickness = None):
    '''Mirror data dogbone 2D and 3D'''
    if (dim == 2):
        dataOut= np.vstack((
        data,
        np.array([0,0]) + data * np.array([1,-1]), #nahoru
        np.array([0,6/4*D*2]) + data * np.array([1,-1]), #dolu
        np.array([0,0]) + data * np.array([-1,1]), #doleva
        np.array([D*2,0]) + data * np.array([-1,1]), #doprava

        np.array([0,0]) + data * np.array([-1,-1]), #nahoru doleva
        np.array([2*D,0]) + data * np.array([-1,-1]), #nahoru doprava
        np.array([0,6/4*2*D]) + data * np.array([-1,-1]), #dolu doleva
        np.array([2*D,6/4*2*D]) + data * np.array([-1,-1]), #dolu doprava
        ))

    if(dim==3):
        dataOut= np.vstack((
        data,
        np.array([0,0,0]) + data * np.array([1,-1,1]), #nahoru c
        np.array([0,6/4*D*2,0]) + data * np.array([1,-1,1]), #dolu c
        np.array([0,0,0]) + data * np.array([-1,1,1]), #doleva c
        np.array([D*2,0,0]) + data * np.array([-1,1,1]), #doprava c

        np.array([0,0,0]) + data * np.array([-1,-1,1]), #nahoru doleva c
        np.array([2*D,0,0]) + data * np.array([-1,-1,1]), #nahoru doprava c
        np.array([0,6/4*2*D,0]) + data * np.array([-1,-1,1]), #dolu doleva c
        np.array([2*D,6/4*2*D,0]) + data * np.array([-1,-1,1]), #dolu doprava c

        np.array([0,0,0]) + data * np.array([1,-1,-1]), #nahoru dopredu
        np.array([0,6/4*D*2,0]) + data * np.array([1,-1,-1]), #dolu dopredu
        np.array([0,0,0]) + data * np.array([-1,1,-1]), #doleva dopredu
        np.array([D*2,0,0]) + data * np.array([-1,1,-1]), #doprava dopredu

        np.array([0,0,0]) + data * np.array([-1,-1,-1]), #nahoru doleva dopredu
        np.array([2*D,0,0]) + data * np.array([-1,-1,-1]), #nahoru doprava dopredu
        np.array([0,6/4*2*D,0]) + data * np.array([-1,-1,-1]), #dolu doleva dopredu
        np.array([2*D,6/4*2*D,0]) + data * np.array([-1,-1,-1]), #dolu doprava dopredu

        np.array([0,0,2*thickness]) + data * np.array([1,-1,-1]), #nahoru dozadu
        np.array([0,6/4*D*2,2*thickness]) + data * np.array([1,-1,-1]), #dolu dozadu
        np.array([0,0,2*thickness]) + data * np.array([-1,1,-1]), #doleva dozadu
        np.array([D*2,0,2*thickness]) + data * np.array([-1,1,-1]), #doprava dozadu

        np.array([0,0,2*thickness]) + data * np.array([-1,-1,-1]), #nahoru doleva dozadu
        np.array([2*D,0,2*thickness]) + data * np.array([-1,-1,-1]), #nahoru doprava dozadu
        np.array([0,6/4*2*D,2*thickness]) + data * np.array([-1,-1,-1]), #dolu doleva dozadu
        np.array([2*D,6/4*2*D,2*thickness]) + data * np.array([-1,-1,-1]), #dolu doprava dozadu

        ))

    """
    dataOut = np.asarray(dataOut)
    fig, ax = plt.subplots()
    ax.scatter(dataOut[:,0], dataOut[:,1])
    if SHOW_PLOT:
        plt.show()
    """
    return dataOut

def mirror_dataCylinder(data, center, radius, height, directionDim, quarter = False):
    data = np.asarray(data)
    rad = radius + 1e-5
    print(directionDim)
    if (directionDim == 0):
        print(quarter)
        if quarter == False:
            dataOut =  np.vstack((data,
                np.array([-1e-5,0,0]) + data * np.array([-1,1,1]),
                np.array([ (height +1e-5)*2 ,0,0]) + data * np.array([-1,1, 1])
                ))

            mirroredData = np.zeros( (len(dataOut)*3) )
            mirroredData = np.reshape ( mirroredData, (len(dataOut),3))

            for i in range (len(mirroredData)):
                rad0 = scipy.spatial.distance.cdist( np.reshape(np.array([dataOut[i,0], center[1], center[2]]), (1,3)), np.reshape(dataOut[i,:], (1,3)))

                mirroredData[i,0] = dataOut[i,0]
                mirroredData[i,1] = center[1] + (-center[1]+dataOut[i,1]) * ((2*rad-rad0) / rad0 )
                mirroredData[i,2] = center[2] + (-center[2]+dataOut[i,2]) * ((2*rad-rad0) / rad0 )

            dataOut =  np.vstack((dataOut, mirroredData))


        if quarter == True:
            dataOut =  np.vstack((data,
                np.array([-1e-5,0,0]) + data * np.array([-1,1,1]),
                np.array([ (height +1e-5)*2 ,0,0]) + data * np.array([-1,1, 1])
                ))

            mirroredData = np.zeros( (len(dataOut)*3) )
            mirroredData = np.reshape ( mirroredData, (len(dataOut),3))

            """
            print('first')
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(dataOut[:,0], dataOut[:,1], dataOut[:,2])
            if SHOW_PLOT:
                plt.show()
            """

            for i in range (len(mirroredData)):
                rad0 = scipy.spatial.distance.cdist( np.reshape(np.array([dataOut[i,0], center[1], center[2]]), (1,3)), np.reshape(dataOut[i,:], (1,3)))

                mirroredData[i,0] = dataOut[i,0]
                mirroredData[i,1] = center[1] + (-center[1]+dataOut[i,1]) * ((2*rad-rad0) / rad0 )
                mirroredData[i,2] = center[2] + (-center[2]+dataOut[i,2]) * ((2*rad-rad0) / rad0 )

            dataOut =  np.vstack((dataOut, mirroredData))

            """
            print('first')
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(dataOut[:,0], dataOut[:,1], dataOut[:,2])
            if SHOW_PLOT:
                plt.show()
            """

            dataOutM =  np.vstack((
                np.array([ 0 ,-1e-4,0]) + dataOut * np.array([1,-1,1]),
                np.array([ 0 ,0,-1e-4]) + dataOut * np.array([1,1,-1]),
                np.array([ 0 ,-1e-4,-1e-4]) + dataOut * np.array([1,1,-1]) *  np.array([1,-1,1])
                ))


            for n in range (len(dataOutM)):
                no = dataOutM[n]
                if (no[1]>-1e-2 and no[2]>-1e-2):
                    dataOut =  np.vstack((dataOut,no))

            """
            print('first')
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(dataOut[:,0], dataOut[:,1], dataOut[:,2])
            if SHOW_PLOT:
                plt.show()
            """


    if (directionDim == 1):
        dataOut =  np.vstack((data,
            np.array([0,0,0]) + data * np.array([1,-1,1]),
            np.array([ 0 , height*2,0]) + data * np.array([1,-1, 1])
            ))
        mirroredData = np.zeros( (len(dataOut)*3) )
        mirroredData = np.reshape ( mirroredData, (len(dataOut),3))

        for i in range (len(mirroredData)):
            rad0 = scipy.spatial.distance.cdist( np.reshape(np.array([center[0], dataOut[i,1], center[2]]), (1,3)), np.reshape(dataOut[i,:], (1,3)))


            mirroredData[i,0] = center[0] + (-center[0]+dataOut[i,0]) * ((2*rad-rad0) / rad0 )
            mirroredData[i,1] = dataOut[i,1]
            mirroredData[i,2] = center[2] + (-center[2]+dataOut[i,2]) * ((2*rad-rad0) / rad0 )

        dataOut =  np.vstack((dataOut, mirroredData))


    if (directionDim == 2):
        dataOut =  np.vstack((data,
            np.array([0,0,0]) + data * np.array([1,1,-1]),
            np.array([ 0 , 0, height*2]) + data * np.array([1,1, -1])
            ))
        mirroredData = np.zeros( (len(dataOut)*3) )
        mirroredData = np.reshape ( mirroredData, (len(dataOut),3))

        for i in range (len(mirroredData)):
            rad0 = scipy.spatial.distance.cdist( np.reshape(np.array([center[0], center[1], dataOut[i,2]]), (1,3)), np.reshape(dataOut[i,:], (1,3)))

            mirroredData[i,0] = center[0] + (-center[0]+dataOut[i,0]) * ((2*rad-rad0) / rad0 )
            mirroredData[i,1] = center[1] + (-center[1]+dataOut[i,1]) * ((2*rad-rad0) / rad0 )
            mirroredData[i,2] = dataOut[i,2]

        dataOut =  np.vstack((dataOut, mirroredData))


    """
    print('first')
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(dataOut[:,0], dataOut[:,1], dataOut[:,2])
    if SHOW_PLOT:
        plt.show()
    """

    return dataOut

def mirror_dataTube(data, center, radius, height, thickness, directionDim):
    data = np.asarray(data)
    outerRad = radius + 1e-2
    innerRad = radius - thickness - 1e-2

    if (directionDim == 0):
        dataOut =  np.vstack((data,
            np.array([-1e-5,0,0]) + data * np.array([-1,1,1]),
            np.array([ (height +1e-5)*2 ,0,0]) + data * np.array([-1,1, 1])
            ))

        mirroredOutside = np.zeros( (len(dataOut)*3) )
        mirroredOutside = np.reshape ( mirroredOutside, (len(dataOut),3))

        for i in range (len(mirroredOutside)):
            rad0 = scipy.spatial.distance.cdist( np.reshape(np.array([dataOut[i,0], center[1], center[2]]), (1,3)), np.reshape(dataOut[i,:], (1,3)))
            mirroredOutside[i,0] = dataOut[i,0]
            mirroredOutside[i,1] = center[1] + (dataOut[i,1]-center[1]) * ((2*outerRad-rad0) / rad0 )
            mirroredOutside[i,2] = center[2] + (dataOut[i,2]-center[2]) * ((2*outerRad-rad0) / rad0 )

        mirroredInside = np.zeros( (len(dataOut)*3) )
        mirroredInside = np.reshape ( mirroredInside, (len(dataOut),3))

        for i in range (len(mirroredInside)):
            rad0 = scipy.spatial.distance.cdist( np.reshape(np.array([dataOut[i,0], center[1], center[2]]), (1,3)), np.reshape(dataOut[i,:], (1,3)))
            mirroredInside[i,0] = dataOut[i,0]
            mirroredInside[i,1] = center[1] + (-center[1]+dataOut[i,1]) * ((2*innerRad-rad0) / rad0 )
            mirroredInside[i,2] = center[2] + (-center[2]+dataOut[i,2]) * ((2*innerRad-rad0) / rad0 )
    """
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(dataOut[:,0], dataOut[:,1], dataOut[:,2])
    ax.scatter(mirroredOutside[:,0], mirroredOutside[:,1], mirroredOutside[:,2])
    ax.scatter(mirroredInside[:,0], mirroredInside[:,1], mirroredInside[:,2])
    if SHOW_PLOT:
        plt.show()
    """
    dataOut =  np.vstack((dataOut, mirroredOutside))
    dataOut =  np.vstack((dataOut, mirroredInside))
    return dataOut


def copy_dataBeam(data, dim, sizes):
    '''Copy data 2D and 3D'''
    if (dim == 2):
        return np.vstack((
            data + np.array([0, 0]),
            data + np.array([-sizes[0], sizes[1]]),
            data + np.array([0, sizes[1]]),
            data + np.array([sizes[0], sizes[1]]),
            data + np.array([-sizes[0], 0]),
            data + np.array([sizes[0], 0]),
            data + np.array([-sizes[0], -sizes[1]]),
            data + np.array([0, -sizes[1]]),
            data + np.array([sizes[0],  -sizes[1]])
            ))

    if (dim == 3):
        return np.vstack((
            data + np.array([0, 0, 0]),
            data + np.array([-sizes[0], sizes[1], 0]),
            data + np.array([0, sizes[1], 0]),
            data + np.array([sizes[0], sizes[1], 0]),
            data + np.array([-sizes[0], 0, 0]),
            data + np.array([sizes[0], 0, 0]),
            data + np.array([-sizes[0], -sizes[1], 0]),
            data + np.array([0, -sizes[1], 0]),
            data + np.array([sizes[0], -sizes[1], 0]),
            data + np.array([0, 0, sizes[2]]),
            data + np.array([-sizes[0], sizes[1], sizes[2]]),
            data + np.array([0, sizes[1], sizes[2]]),
            data + np.array([sizes[0], sizes[1], sizes[2]]),
            data + np.array([-sizes[0], 0, sizes[2]]),
            data + np.array([sizes[0], 0, sizes[2]]),
            data + np.array([-sizes[0], -sizes[1], sizes[2]]),
            data + np.array([0, -sizes[1], sizes[2]]),
            data + np.array([sizes[0], -sizes[1], sizes[2]]),
            data + np.array([0, 0, -sizes[2]]),
            data + np.array([-sizes[0], sizes[1], -sizes[2]]),
            data + np.array([0, sizes[1], -sizes[2]]),
            data + np.array([sizes[0], sizes[1], -sizes[2]]),
            data + np.array([-sizes[0], 0, -sizes[2]]),
            data + np.array([sizes[0], 0, -sizes[2]]),
            data + np.array([-sizes[0], -sizes[1], -sizes[2]]),
            data + np.array([0, -sizes[1], -sizes[2]]),
            data + np.array([sizes[0], -sizes[1], -sizes[2]])
            ))


def voronoi_2d(vor, sizes, shifts=0):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    polygons = []
    new_regions = []
    areas = []
    centroids = []
    points = []
    new_vertices = vor.vertices.tolist()

    polArr = np.array([(0.,0), (0,sizes[1]), (sizes[0],sizes[1]), (sizes[0],0)])
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

def voronoi_2d_dogBone(vor, sizes):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    volumes = []
    #
    for i in range (vor.points.shape[0]):
        comp = True
        for p in range (2):
            if (vor.points[i][p] < 0 or vor.points[i][p] > sizes[p] ):
                comp = False

        if (comp):
            vol = cellVolume3d(vor, i)
            volumes.append(vol)

    return volumes

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


def volumesCylinder3d (vor, center, radius, height, directionDim ):
    if vor.points.shape[1] != 3:
        raise ValueError("Requires 3D input")

    volumes = []
    #
    """
    for i in range (vor.points.shape[0]):
        comp = True
        for p in range (3):
            if (vor.points[i][p] < 0 or vor.points[i][p] > sizes[p] ):
                comp = False

        if (comp):
            vol = cellVolume3d(vor, i)
            volumes.append(vol)
    """
    return volumes
