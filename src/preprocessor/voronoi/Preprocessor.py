#!/usr/bin/env python
# coding: utf-8
#
# In[1]:

import os
import sys
import time
import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from IPython.display import clear_output

import math
from sklearn import preprocessing

from scipy.ndimage import rotate

#2d voronoi a teselace
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay

from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix

import utilitiesGeom
import utilitiesMech
import utilitiesModeling
import utilitiesNumeric
import voronoi


print('\n%%%%%%%%% LATTICE PREPROCESSOR STARTED %%%%%%%%%')
start = time.time()

try:
    if not os.path.exists('inpFiles'):
        os.makedirs('inpFiles')
except:
    print('Please create inpFiles directory! Code Exited.')
    sys.exit()

#type of solver. does not matter now
solver = 0

#power tesselation on/off
powerTes = 0

#dimension 2/3
dim = 2
print('Creating a %dd lattice model...' %dim)

#dimensions of a rectangle model
if (dim == 2 ): maxLim = np.array([3,1])
if (dim == 3 ): maxLim = np.array([1,1,1])

#volume of the model (later for check)
volume = np.sum(maxLim)

#size of grains (minimum distance between nodes)
radius = 0.2
minDist = radius

#trials of random node positioning
trials = 5000

#lists for the model
node_coords = []
node_mechBC = []
mechBC_merged = []
transportBC_merged = []
functions = []

#### Defining functions
#0 constant zero
func = []
func.append( np.array([0, 0]) )
fn = utilitiesMech.generalFunc(func)
functions.append (fn)

#1 constant one
func1 = []
func1.append( np.array([0, 1]) )
fn1 = utilitiesMech.generalFunc(func1)
functions.append (fn1)

#2 jen zkusebni funkce pro transportni boundary conditions
transtBC = np.array([0,1])
tsptBC = utilitiesGeom.transportBC(-20, transtBC)
transportBC_merged.append(tsptBC)


#sampling of nodes
if (dim == 2):
    node_coords,node_mechBC, mechBC_merged = utilitiesModeling.create2dCantileverUniTens(maxLim, minDist, trials )
if (dim == 3):
    print('3d model inactive! Exiting.')
    sys.exit()
    #assemble3Dblock()


node_coords = np.asarray(node_coords)
node_count = len(node_coords)
print('Model containing %d nodes successfuly generated.' %(node_count))


print('Conducting Voronoi tesselation...')
####### Mirroring data
#Clean sample
#vor = Voronoi(node_coords)

#mirror_data
vor = Voronoi(voronoi.mirror_dataBeam(node_coords, dim, maxLim))

#copy_data
#vor = Voronoi(voronoi.copy_dataBeam(node_coords, dim, maxLim))


### Conducting Voronoi Tess
if (dim == 2 ):
    regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)

    # fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange', line_width=1, line_alpha=1, point_size=0.5)
    # fig.set_figwidth(maxLim[0]*3)
    # fig.set_figheight(maxLim[1]*3)
    #
    # for i in range(points.shape[0]):
    #     plt.text(points[i,0], points[i,1], str(i))
    #     #print (i)
    #
    # plt.xlim(-maxLim[0]*0.2, maxLim[0]*1.2)
    # plt.ylim(-maxLim[1]*0.2, maxLim[1]*1.2)
    # plt.show()

if (dim == 3 ):
    areas = voronoi.voronoi_3d(vor, maxLim)

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    #
    # d = node_coords
    # xcoords = d[:,0]
    # ycoords = d[:,1]
    # zcoords = d[:,2]
    #
    # ax.scatter3D(xcoords,ycoords,zcoords)
    #
    # ax.set_xlim3d(-maxLim[0]*1,2*maxLim[0])
    # ax.set_ylim3d(-maxLim[1]*1,2*maxLim[1])
    # ax.set_zlim3d(-maxLim[2]*1,2*maxLim[2])
    #
    # plt.show()


#reordering nodes due to their connectivity
order = utilitiesNumeric.reorderToDiagonal(node_count, node_coords, vor)


mechanicalElements = []
transportPaths = []
materials = []

vert_count = -1

young = 30e9
poisson = 0.3
transpC = 11
transpS = 22
density = 2200
linElMaterial = utilitiesMech.linearElasticMaterial(young, poisson, transpC, transpS, density)
materials.append(linElMaterial)



print('Saving model...')
if (dim == 2):
    vert_count = utilitiesGeom.output2D(node_count, dim, maxLim, vor, node_coords, node_mechBC,
    areas, order, mechanicalElements, mechBC_merged, transportPaths,materials, functions, False)
if (dim == 3):
    vert_count = utilitiesGeom.output3D(node_count, dim, maxLim, vor, node_coords, node_mechBC,
    areas, order, mechanicalElements, mechBC_merged, transportPaths, materials, functions, False)


solStep = 1e-4

print('Saving boundary conditions...')
utilitiesGeom.saveMechBC(dim, mechBC_merged)
utilitiesGeom.saveTransportBC(transportBC_merged)

print('Saving master file...')
utilitiesGeom.saveMasterInput(dim, solver, solStep)



end =  time.time() -start
print('All done in %.3f secs.' %end)
print('%%%%%%%%% LATTICE PREPROCESSOR FINISHED %%%%%%%%%\n')
