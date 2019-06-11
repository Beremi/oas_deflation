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

#type of solver. does not matter now
solver = 0

#power tesselation on/off
powerTes = 0

#dimension 2/3
dim = 2
print('Creating a %dd lattice model...' %dim)

#dimensions of a rectangle model
if (dim == 2 ): maxLim = np.array([  5   ,   2  ])
if (dim == 3 ): maxLim = np.array([1,1,1])

#volume of the model (later for check)
volume = np.sum(maxLim)

#size of grains (minimum distance between nodes)
#be cautious with small grains
radius = 0.09
minDist = radius

#trials of random node positioning
trials = 20000

#lists for the model
node_coords = []
node_mechBC = []
mechBC_merged = []
transportBC_merged = []
functions = []

#### Defining functions
#0 constant zero
fn = utilitiesNumeric.constantFunc(0)
functions.append (fn)

#1 loading function, single force top right, bilinear
func1 = []
func1.append( np.array([0,0]) )
func1.append( np.array([50, 50e4]) )
fn1 = utilitiesNumeric.generalFunc(func1)
functions.append (fn1)

#transport function, leftFace, constant
fn2 = utilitiesNumeric.constantFunc(20)
functions.append (fn2)

#transport function, rightFace, bilinear
func3 = []
func3.append( np.array([0,0]) )
func3.append( np.array([50, 500]) )
fn3 = utilitiesNumeric.generalFunc(func3)
functions.append (fn3)

#sampling of nodes
if (dim == 2):
    #cantilever
    #node_coords,node_mechBC, mechBC_merged, transportBC_merged, vor, areas   = utilitiesModeling.create2dCantileverBending(maxLim, minDist, trials )
    #simply supported beam
    node_coords,node_mechBC, mechBC_merged, transportBC_merged, vor, areas   = utilitiesModeling.create2dSSBeamUnifLoad(maxLim, minDist, trials )

if (dim == 3):
    print('3d model inactive! Exiting.')
    sys.exit()
    #assemble3Dblock()

node_coords = np.asarray(node_coords)
node_count = len(node_coords)
print('Model containing %d nodes successfuly generated.' %(node_count))

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



print('\nSaving model...')
if (dim == 2):
    vert_count, verticesIdxDict, vertIdxStart = utilitiesGeom.output2D(node_count, dim, maxLim,
    vor, node_coords, node_mechBC,
    areas, order, mechanicalElements, mechBC_merged, transportPaths, False)
if (dim == 3):
    vert_count = utilitiesGeom.output3D(node_count, dim, maxLim, vor, node_coords, node_mechBC,
    areas, order, mechanicalElements, mechBC_merged, transportPaths, False)


utilitiesGeom.saveMaterials(materials)
utilitiesGeom.saveFunctions(functions)
utilitiesGeom.saveMechBC(dim, mechBC_merged)
utilitiesGeom.saveTransportBC(transportBC_merged, verticesIdxDict, vertIdxStart)
utilitiesGeom.saveExporters()

solStep = 10
utilitiesGeom.saveMasterInput(dim, solver, solStep)


end =  time.time() -start
print('\nAll done in %.3f secs.' %end)
print('%%%%%%%%% LATTICE PREPROCESSOR FINISHED %%%%%%%%%\n')
