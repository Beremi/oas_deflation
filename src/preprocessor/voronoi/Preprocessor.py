#!/usr/bin/env python
# coding: utf-8

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
import voronoi


print('\n%%%%%%%%% LATTICE PREPROCESSOR STARTED %%%%%%%%%')
start = time.time()

try:
    if not os.path.exists('/inpFiles/'):
        os.makedirs('/inpFiles/')
except:
    print('Please create inpFiles directory! Code Exited.')
    sys.exit()

#type of solver. does not matter now
solver = 0

#power tesselation on/off
powerTes = 0

#dimension 2/3
dim = 2
print('Creating a %dD lattice model...' %dim)

#dimensions of a rectangle model
if (dim == 2 ): maxLim = np.array([3,1])
if (dim == 3 ): maxLim = np.array([3,2,1])

#volume of the model (later for check)
volume = np.sum(maxLim)

#size of grains (minimum distance between nodes)
radius = 0.15
minDist = radius

#trials of random node positioning
trials = 1000


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

#just a test funtion for trsprtBC
transtBC = np.array([0,1])
tsptBC = utilitiesGeom.transportBC(-20, transtBC)
transportBC_merged.append(tsptBC)



######## FUNCTION FOR CREATING OF A 2D SUPPORTED RECTANGLE MODEL
def assemble2DRectangle ():
    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    ###############generating of nodes, supported line left vertical ###############
    #mech bc
    lineBC = np.array([0,0,0,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of nodes, a single point top right (a line of zero length) ###############
    lineBC = np.array([-1,-1,-1,0,1,0])

    #defining points of the line
    nodeA = np.array([maxLim[0] - indent, maxLim[1] - indent])
    nodeB = np.array([maxLim[0] - indent, maxLim[1] - indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])

    #rect
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, rectBC)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

######## FUNCTION FOR CREATING OF A 3D SUPPORTED RECTANGE
def assemble3Dblock():
    indent = 1e-8

    transtBC = np.array([0,1])
    tsptBC = utilitiesGeom.transportBC(20, transtBC)
    transportBC_merged.append(tsptBC)

    ###############generating of points supported line bottom left ###############
    mechBC = np.array([0,0,0,0,0,0,0,0,0,0,0,0])

    nodeA = np.array([indent , indent, indent])
    nodeB = np.array([indent , maxLim[1] - indent, indent])

    oldLen = len(node_coords)
   # utilitiesGeom.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials, True)
    #
    nrOfPoints =  (len(node_coords)) - oldLen

    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of points supported line top right ###############
    mechBC = np.array([-1,-1,-1, -1,-1,-1,    0,1,0, 0,0,0])

    nodeA = np.array([maxLim[0] - indent , indent, maxLim[2] -indent])
    nodeB = np.array([maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials, True)
    #
    nrOfPoints =  (len(node_coords)) - oldLen
     #print (nrOfPoints)

    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')


    ###############generating of points supported surface  left face ###############
    mechBC = np.array([0,0,0,0,0,0,-1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')


    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    #kvadr
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, mechBC)
    newLen = len(node_coords)-1
    nrOfPoints =  (len(node_coords)) - oldLen
   # for n in range ( nrOfPoints ):
       # mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
       # mechBC_merged.append(mBC)
       # print('adding')
    ####################################################################################################


# In[5]:


#sampling of nodes
if (dim == 2):
    assemble2DRectangle()
if (dim == 3):
   assemble3Dblock()


node_coords = np.asarray(node_coords)
node_count = len(node_coords)

clear_output(True)
print('Model containing %d nodes successfuly generated.' %(node_count))




print('Conducting Voronoi tesselation...')
####### MIRRORING OR COPYING OF POINTS FOR VORONOI TESS

#zrcadlení/rozkopírování bodů

#cisty sampl
#vor = Voronoi(node_coords)

#mirror_data - Vaškovo zrcadlení bodů
#vor = Voronoi(voronoi.mirror_data(node_coords))
#upraveno pro obecny obdelnik
vor = Voronoi(voronoi.mirror_dataBeam(node_coords, dim, maxLim))

#copy_data - Vaškovo rozkopírování bodů
#vor = Voronoi(voronoi.copy_data(node_coords))
#upraveno pro obecny obdelnik
#vor = Voronoi(voronoi.copy_dataBeam(node_coords, dim, maxLim))

if (dim == 2 ):
    regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)

    fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange', line_width=1, line_alpha=1, point_size=0.5)

    fig.set_figwidth(maxLim[0]*3)
    fig.set_figheight(maxLim[1]*3)

    for i in range(points.shape[0]):
        plt.text(points[i,0], points[i,1], str(i))
        #print (i)

    plt.xlim(-maxLim[0]*0.2, maxLim[0]*1.2)
    plt.ylim(-maxLim[1]*0.2, maxLim[1]*1.2)
    #plt.xlim(0, 1)
    #plt.ylim(0, 1)

    #plt.show()

if (dim == 3 ):
    areas = voronoi.voronoi_3d(vor, maxLim)

    fig = plt.figure()
    ax = plt.axes(projection='3d')

   # d = voronoi.copy_dataBeam(node_coords, dim, maxLim)
   # d = voronoi.mirror_dataBeam(node_coords, dim, maxLim)
    d = node_coords
    xcoords = d[:,0]
    ycoords = d[:,1]
    zcoords = d[:,2]

    ax.scatter3D(xcoords,ycoords,zcoords)

   #ax.scatter3D(node_coords[:,0],node_coords[:,1],node_coords[:,2])

    ax.set_xlim3d(-maxLim[0]*1,2*maxLim[0])
    ax.set_ylim3d(-maxLim[1]*1,2*maxLim[1])
    ax.set_zlim3d(-maxLim[2]*1,2*maxLim[2])
   #& plt.axis('equal')

   # ax.set_xlim3d(-maxLim[0],2*maxLim[0])
   # ax.set_ylim3d(-maxLim[0],2*maxLim[0])
   # ax.set_zlim3d(-maxLim[0],2*maxLim[0])

    #plt.show()

#ptA = np.array([0,0])
#ptB = np.array([0.5,0.7])

#id = utilitiesGeom.returnSelectedPts (ptA, ptB, vor.points)

#print (id)

#for i in range (len(id)):#
#    print (vor.points[id[i]])


# In[8]:


# spatial boundary conditions
# najit vertexy na hornim lici
#ptA = np.array([0, 0.99])
#ptB = np.array([3, 1])

#id = utilitiesGeom.returnSelectedPts (ptA, ptB, vor.vertices)

#print (id)

#for i in range (len(id)):#
#    print (vor.vertices[id[i]])


print('Checking connectivitiy...')
# matice "konektivity" nodů
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
plt.imshow(A)
#plt.colorbar()
#plt.show()

print('Diagonalizing...')
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
plt.imshow(B)
#plt.colorbar()
#plt.show()


#print (len(mechBC_merged))
#for d in mechBC_merged:
#    d.printProps()


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


diagonalize = False

#print (len(points))

print('Saving model...')
if (dim == 2):
    vert_count = utilitiesGeom.output2D(node_count, dim, maxLim, vor, node_coords, node_mechBC,
    areas, order, mechanicalElements, mechBC_merged, transportPaths,materials, functions, diagonalize)
if (dim == 3):
    vert_count = utilitiesGeom.output3D(node_count, dim, maxLim, vor, node_coords, node_mechBC,
    areas, order, mechanicalElements, mechBC_merged, transportPaths, materials, functions, diagonalize)

#print (vert_count)

#vert_count = 223
mechSolStatDyn = 0
mechSolImpExp = 0
solStep = 1e-4

print('Saving boundary conditions...')
utilitiesGeom.saveMechBC(dim, mechBC_merged)
utilitiesGeom.saveTransportBC(transportBC_merged)

print('Saving master file...')
utilitiesGeom.saveMasterInput(dim, solver, solStep)



# In[ ]:


#for i in range (len(mechanicalElements)):
 #   print(mechanicalElements[i].getString())

#for item in transportPaths:
   # print (item.getString())

#print(mechBC_merged)

#for d in range (len(mechBC_merged)):
#    mechBC_merged[d].printProps()

#print()
#for item in materials:#
    #print(item.getString())



# In[ ]:


#ptA = np.array([0,0])
#ptB = np.array([0.5,0.7])

#id = utilitiesGeom.returnSelectedPts (ptA, ptB, vor.points)

#print (id)

#for i in range (len(id)):
#    print (vor.points[id[i]])


end =  time.time() -start
print('All done in %.3f secs.' %end)
print('%%%%%%%%% LATTICE PREPROCESSOR FINISHED %%%%%%%%%\n')
