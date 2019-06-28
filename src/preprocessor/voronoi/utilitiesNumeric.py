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

#run voronoi, mirrored data
def runMirroredVoronoi (node_coords, dim, maxLim):
    vor = Voronoi(voronoi.mirror_dataBeam(node_coords, dim, maxLim))

    return vor



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
    def __init__ (self, amplitude, freq):
        self.amplitude = amplitude
        self.freq = freq
    def getString(self):
        line = 'SineFunction\t'
        line += '%f\t' %(self.freq)
        line += '%f' %(self.amplitude)
        return line
##################################################




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
