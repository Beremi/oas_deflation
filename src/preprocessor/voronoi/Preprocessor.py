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

Xdim = 1.
Ydim = 1.
Zdim = 1.

#dimensions of a rectangle model
if (dim == 2 ): maxLim = np.array([  Xdim   ,  Ydim ])
if (dim == 3 ): maxLim = np.array([  Xdim,  Ydim,  Zdim ])


#volume of the model (later for check)
volume = np.sum(maxLim)

#size of grains (minimum distance between nodes)
#be cautious with small grains!
minDist = 0.02
radius = minDist / 2

elaX = minDist / Xdim * 2

if (dim == 2):
    dV = 3.141592 * radius **2
if (dim == 3):
    dV = 4/3 * 3.141592 * radius **3

expNodes = volume / dV  * 0.6
print ('Expecting about %d nodes' %expNodes)

#trials of random node positioning
trials = 50000

#lists for the model
node_coords = []
mechBC_merged = []
mechIC_merged = []
trsprtBC_merged = []
trsprtIC_merged = []
functions = []


materialZones = []
#matZone 1
matZ = []
if (dim==2):
    boundA = np.array(  [ -1e-8             , -1e-8          ] )
    matZ.append (boundA)
    boundB = np.array(  [ maxLim[0]*elaX    , maxLim[1] + 1e-8] )
    matZ.append (boundB)
    boundA1 = np.array(  [ maxLim[0]-maxLim[0]*elaX , - 1e-8] )
    matZ.append (boundA1)
    boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8]  )
    matZ.append (boundB1)
    materialZones.append(matZ)
if (dim==3):
    boundA = np.array(  [ -1e-8             , -1e-8             , -1e8] )
    matZ.append (boundA)
    boundB = np.array(  [ maxLim[0]*elaX    , maxLim[1] + 1e8   , maxLim[2] + 1e8  ] )
    matZ.append (boundB)
    boundA1 = np.array(  [ maxLim[0]-maxLim[0]*elaX , - 1e-8    , -1e8] )
    matZ.append (boundA1)
    boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8   , maxLim[2] + 1e8 ]  )
    matZ.append (boundB1)
    materialZones.append(matZ)

#creating the model. Select the prepared models.
if (dim == 2):
    #cantilever bending
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dCantileverBending(maxLim, minDist, trials )

    #cantilever  pressure free contraction
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dCantileverUniTens(maxLim, minDist, trials)

    #confined  pressure
    node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions  = utilitiesModeling.create2dbeamConfinedPress(maxLim, minDist, trials )

    #simply supported beam, uniform load
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dSSBeamUnifLoad(maxLim, minDist, trials )

    #single spring test
    #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions = utilitiesModeling.createSingleSpringTestModel( 2 )

    #diamond test
    #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions = utilitiesModeling.createDiamondTestModel(1, 2)

if (dim == 3):
    #cantilever bending
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverBending(maxLim, minDist, trials )

    #cantilever uniform pressure, free contraction
    node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverUniPressFree(maxLim, minDist, trials )

    #cantilever uniform pressure, confined
    #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverUniPressConfined(maxLim, minDist, trials )




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
density = 2200


ft = 2e6
Gt = 500
marsMaterial = utilitiesMech.MarsMaterial(young, poisson, density, ft, Gt)
materials.append(marsMaterial)

#	E0	43.0e9	alpha	0.300000    density 2200.0 tauBar 4.0e6 Kin 0.0 gamma 10.0e6 S 0.0025e6 m 0
fatigueMaterial = utilitiesMech.FatigueMaterial(  43.0e9, 0.300000 , 2200.0, 4.0e6, 0.0, 10.0e6 , 0.0025e6, 0)
#materials.append(fatigueMaterial)



transpC = 11
transpS = 22
transportMaterial = utilitiesMech.TransportMaterial( transpC, transpS)
materials.append(transportMaterial)

linElMaterial = utilitiesMech.linearElasticMaterial(young, poisson, density)
materials.append(linElMaterial)


print('')


#Deconstructing Voronoi diagram and saving the geometry
vert_count, verticesIdxDict, vertIdxStart = utilitiesGeom.extractGeometry(dim, node_count,  maxLim, vor, node_coords, areas, mZ=materialZones)


# saving rest of input
utilitiesGeom.saveMaterials(materials)
utilitiesGeom.saveFunctions(functions)
utilitiesGeom.saveMechBC(dim, mechBC_merged)
if (len(mechIC_merged)>0):  utilitiesGeom.saveMechIC(dim, mechIC_merged)
utilitiesGeom.saveTransportBC(trsprtBC_merged, verticesIdxDict, vertIdxStart)
if (len(trsprtIC_merged)>0):utilitiesGeom.saveTransportIC(trsprtIC_merged)
utilitiesGeom.saveExporters()

solStep = 1e-2
simTime = 10
utilitiesGeom.saveMasterInput(dim, solver, solStep, simTime)
end =  time.time() -end
print('Saving done in %.3f secs.' %end)

end =  time.time() -start
print('\nAll done in %.3f secs.' %end)
print('%%%%%%%%% LATTICE PREPROCESSOR FINISHED %%%%%%%%%\n')
