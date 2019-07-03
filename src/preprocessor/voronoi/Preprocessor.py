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

#power tesselation on/off  does not matter now
powerTes = 0

#dimension 2 //// dim 3 prohibited now
dim = 2
print('Creating a %dd lattice model...' %dim)

Xdim = 5.
Ydim = 3.
Zdim = 1.

#dimensions of a rectangle model
if (dim == 2 ): maxLim = np.array([  Xdim   ,  Ydim ])
if (dim == 3 ): maxLim = np.array([  Xdim,  Ydim,  Zdim ])


#volume of the model (later for check)
volume = np.sum(maxLim)

#size of grains (minimum distance between nodes)
#be cautious with small grains!
minDist = 0.1
radius = minDist / 2

if (dim == 2):
    dV = 3.141592 * radius **2
if (dim == 3):
    dV = 4/3 * 3.141592 * radius **3

expNodes = volume / dV  * 0.6
print ('Expecting about %d nodes' %expNodes)

#trials of random node positioning
trials = 30000

#lists for the model
node_coords = []
mechBC_merged = []
mechIC_merged = []
trsprtBC_merged = []
trsprtIC_merged = []
functions = []



#creating the model. Select the prepared models.
if (dim == 2):
    #cantilever bending
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dCantileverBending(maxLim, minDist, trials )

    #cantilever uni tension
    #do not use this yet! ...please
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dCantileverUniTens(maxLim, minDist, trials )

    #simply supported beam, uniform load
    node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dSSBeamUnifLoad(maxLim, minDist, trials )

    #single spring test
    #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions = utilitiesModeling.createSingleSpringTestModel( 2 )

    #diamond test
    #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions = utilitiesModeling.createDiamondTestModel(1, .7)

if (dim == 3):
    #cantilever
    node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverBending(maxLim, minDist, trials )

node_coords = np.asarray(node_coords)
node_count = len(node_coords)
print('Model containing %d nodes successfuly generated.' %(node_count))
end =  time.time() -start
print('Model done in %.3f secs.' %end)
end=time.time()
sys.stdout.flush()

#reordering nodes due to their connectivity
#order = utilitiesNumeric.reorderToDiagonal(node_count, node_coords, vor)

materials = []

vert_count = -1

young = 30e9
poisson = 0.3
transpC = 11
transpS = 22
density = 2200
linElMaterial = utilitiesMech.linearElasticMaterial(young, poisson, transpC, transpS, density)
materials.append(linElMaterial)

print('')


#Deconstructing Voronoi diagram and saving the geometry
vert_count, verticesIdxDict, vertIdxStart = utilitiesGeom.extractGeometry(dim, node_count,  maxLim, vor, node_coords, areas)
# saving rest of input
utilitiesGeom.saveMaterials(materials)
utilitiesGeom.saveFunctions(functions)
utilitiesGeom.saveMechBC(dim, mechBC_merged)
if (len(mechIC_merged)>0):  utilitiesGeom.saveMechIC(dim, mechIC_merged)
utilitiesGeom.saveTransportBC(trsprtBC_merged, verticesIdxDict, vertIdxStart)
if (len(trsprtIC_merged)>0):utilitiesGeom.saveTransportIC(trsprtIC_merged)
utilitiesGeom.saveExporters()

solStep = 10
utilitiesGeom.saveMasterInput(dim, solver, solStep)
end =  time.time() -end
print('Saving done in %.3f secs.' %end)

end =  time.time() -start
print('\nAll done in %.3f secs.' %end)
print('%%%%%%%%% LATTICE PREPROCESSOR FINISHED %%%%%%%%%\n')
