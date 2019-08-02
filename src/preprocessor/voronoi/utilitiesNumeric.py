import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import matplotlib.pylab as plt
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import voronoi

#2d voronoi a teselace
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay

##run voronoi, mirrored data
def runMirroredVoronoi (node_coords, dim, maxLim, shifts=0):
    vor = Voronoi(voronoi.mirror_dataBeam(node_coords, dim, maxLim, shifts))

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim, shifts = shifts)
        return vor, regions, vertices, polygons, areas, centroids, points
    if (dim == 3):
        volumes = voronoi.voronoi_3d(vor, maxLim);
        return vor, volumes

def runCylinderMirroredVoronoi  (node_coords, center, radius, height, directionDim):
    vor = Voronoi(voronoi.mirror_dataCylinder(node_coords, center, radius, height, directionDim))
    volumes = voronoi.volumesCylinder3d (vor, center, radius, height, directionDim )
    return vor, volumes

##################################################
#### General function set by table ####
class generalFunc:
    def __init__(self, table):
        self.table = table

    def getString(self):
        line = 'PWLFunction\t%d'%(len(self.table))

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][0])

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][1])

        return line
####################################################
class PWLFuncFromTxt:
    def __init__(self, filename):
        self.table = np.loadtxt(filename)

    def getString(self):
        line = 'PWLFunction\t%d'%(len(self.table))

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][0])

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][1])

        return line


##################################################
#### General function set by table ####
class constantFunc:
    def __init__(self, val):
        self.val = val

    def getString(self):
        line = 'PWLFunction\t1'
        line += '\t%f'%(0)
        line += '\t%f'%(self.val)

        return line
####################################################

##################################################
#### Sine function ####
class sineFunc:
    def __init__ (self, value, period, shift = None):
        self.value = value
        self.period = period
        self.shift = shift
    def getString(self):
        line = 'SinusFn'
        line += '\tvalue\t%f' %(self.value)
        line += '\tperiod\t%f' %(self.period)
        if (self.shift!=None):
            line += '\tshift\t%f' %(self.shift)

        return line
##################################################

##################################################
#### Saw tooth constant function ####
class sawToothConstFunc:
    def __init__(self, value, period =None, sym=None, lower=None, time=None, num_cycles=None):
        self.value = value
        self.period = period
        self.sym = sym
        self.lower = lower
        self.time = time
        self.num_cycles = num_cycles
    def getString (self):
        line = 'ConstSawToothFn'
        line += '\tvalue\t%f' %(self.value)
        if (self.period!=None):
            line += '\tperiod\t%f' %(self.period)
        if (self.sym!=None):
            line += '\tsym\t%f' %(self.sym)
        if (self.lower!=None):
            line += '\tlower\t%f' %(self.lower)
        if (self.time!=None):
            line += '\ttime\t%f' %(self.time)
        if (self.num_cycles!=None):
            line += '\tlower\t%f' %(self.num_cycles)

        return line
##################################################

##################################################
#### Linear saw tooth constant function ####
class linearSawToothFunc:
    def __init__(self, value, period=None, sym=None, lower=None, time=None, num_cycles=None, multiplier = None):
        self.value = value
        self.period = period
        self.sym = sym
        self.lower = lower
        self.time = time
        self.num_cycles = num_cycles
        self.multiplier = mutiplier
    def getString (self):
        line = 'LinSawToothFn'
        line += '\tvalue\t%f' %(self.value)
        if (self.period!=None):
            line += '\tperiod\t%f' %(self.period)
        if (self.sym!=None):
            line += '\tsym\t%f' %(self.sym)
        if (self.lower!=None):
            line += '\tlower\t%f' %(self.lower)
        if (self.time!=None):
            line += '\ttime\t%f' %(self.time)
        if (self.num_cycles!=None):
            line += '\tlower\t%f' %(self.num_cycles)
        if (self.multiplier!=None):
            line += '\tmultiplier\t%f' %(self.multiplier)
        return line

##################################################

##################################################
#### varying saw tooth function
class varyingSawToothFunction:
     def __init__(self, pwlFn, constSawToothFn ):
         self.constSawToothFn = constSawToothFn
         self.pwlFn = pwlFn
     def getString (self):
         line = 'VaryingSawToothFn\t'
         line += self.pwlFn.getString()[12:]
         line += '\t'
         line += self.constSawToothFn.getString()[16:]
         return line




#reordering of indices
def reorderToDiagonal (node_count, node_coords, vor):
    A = np.zeros( (node_count,node_count) )

    validRidgeIdxs = []
    for i in range (vor.ridge_points.shape[0]):
            pr = False
            if (vor.ridge_points[i][0] < node_count and vor.ridge_points[i][1] < node_count):
                    pr=True
            if (pr):
                #print(vor.ridge_points[i,:])
                validRidgeIdxs.append(i)

    validRidgeIdxs = np.asarray(validRidgeIdxs)

    for i in range (validRidgeIdxs.size):
        pointA = vor.ridge_points[validRidgeIdxs[i]][0]
        pointB = vor.ridge_points[validRidgeIdxs[i]][1]

        coordsA = node_coords[pointA,:]
        coordsB = node_coords[pointB,:]

        dist = np.linalg.norm( coordsB - coordsA )

        A [pointA][pointB] = 1
        A [pointB][pointA] = 1
        #print(dist)


    #print('original connectivity matrix')
    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    #plt.imshow(A)
    #plt.colorbar()
    #plt.show()

    C = np.zeros( (node_count,node_count) )
    C = csr_matrix(A)
    order = reverse_cuthill_mckee(C, symmetric_mode=True)
    order = np.asarray(order)

    B = A[order][:,order]
    #print(B)


    #print('reordered connectivity matrix')
    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    #plt.imshow(B)
    #plt.colorbar()
    #plt.show()

    return order
