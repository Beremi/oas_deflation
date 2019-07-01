import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os
import utilitiesMech
import utilitiesGeom





def randPointInRectangle(dim, maxLim):
    coords = np.random.random(dim)
    coords *= maxLim
    return coords

def randPointOnLine(dim, nodeA, nodeB):
    coords = np.zeros(dim)
    r = np.random.uniform()
    coords = (nodeB - nodeA)*r +nodeA
    return coords


# generates random points no closer to each other than minDist
# into 2d or 3d block
# maxLim: n-d array of dimensions
def generateNodesRect(maxLim, minDist, dim, trials, node_coords):
    if (dim==2):
        print('Generating 2d block segment of size: %f / %f. This may take few minutes. Do not panic. \n Attempt to use the Cython solution (Vasek) !!!' %(maxLim[0], maxLim[1]) )
    if (dim==3):
        print('Generating 3d block segment of size: %f / %f / %f. This may take long. Keep calm.' %(maxLim[0], maxLim[1], maxLim[2]) )

    tr = 0
    while (tr<trials):
        tr = 0
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInRectangle(dim, maxLim)
            distIsGood = True
            #
            #distIsGood = utilitiesGeom.checkMutDistancesLoops(dim, minDist, node_coords, coords)
            distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, node_coords, coords)
            #distIsGood = cutilitiesGeom.heckMutDistancesCKDTree(dim, minDist, node_coords, coords)

            if (distIsGood == False):
                tr += 1

            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords:
        if (tr < trials):
            node_coords.append(coords)

try:
    from point_generators_cython import generateNodesRect_cython as generateNodesRect
    print('Using Cython version of point generator - generateNodesRect.')
except:
    print('''Using Python version of generator. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')

#generates random points onto a set 3d line. No closer than minDst
#catch corners: samples the boundary points first
def generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners):
    print('Generating 3d line segment from [%f; %f; %f] to [%f; %f; %f] '
     %(nodeA[0], nodeA[1],nodeA[2],nodeB[0], nodeB[1],nodeB[2]) )

    if(catchCorners):
        node_coords.append(np.copy(nodeA))
        node_coords.append(np.copy(nodeB))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointOnLine(dim, nodeA, nodeB)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)



def generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials):
    print('Generating 3d surface segment from [%f; %f; %f] to [%f; %f; %f]'
     %(nodeA[0], nodeA[1],nodeA[2],nodeB[0], nodeB[1],nodeB[2]) )

    tr=0
    while (tr<trials):
        tr = 0;
        distIsGood = False
        while (distIsGood == False):
            coords = np.zeros(dim)
            #
            for c in range (dim):
                if (nodeA[c] == nodeB[c]):
                    coords[c] = nodeA[c]
                else:
                    coords[c] = (nodeB[c] - nodeA[c])*np.random.uniform()  + nodeA[c]

            distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, node_coords, coords)
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)


#generate a single node in 2d or 3d
def generateSingleNode(node, dim, node_coords):
    if (dim == 2):
        print('Generating a single 2d node [%f; %f]' %(node[0], node[1]))
    if (dim == 3):
        print('Generating a single 3d node [%f; %f; %f]' %(node[0], node[1], node[2]))

    node_coords.append(node)



#generates random points onto a set 3d line. No closer than minDst
#catch corners: samples the boundary points first
#equid: you can generate equidistant points on the line
def generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners, equidist):
    print('Generating 2d line segment from [%f; %f] to [%f; %f] '
     %(nodeA[0], nodeA[1], nodeB[0], nodeB[1]) )

    if(catchCorners):
        node_coords.append(np.copy(nodeA))
        node_coords.append(np.copy(nodeB))

    if (not equidist):
        tr=0
        while (tr<trials):
            tr = 0;
            distIsGood = False
            while (distIsGood == False):
                coords = randPointOnLine(dim, nodeA, nodeB)

                distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, node_coords, coords)
                if (distIsGood == False):
                    tr += 1
                if (tr > trials): break
            #
            #Adding node coords
            if (tr < trials):
                node_coords.append(coords)

    else:
        #print('Equid')
        mD = minDist * 1.6
        length = np.linalg.norm(nodeA - nodeB)
        nodeNr = int (length / mD)
        indnt = (length-(nodeNr-1)*mD) / 2
        for i in range (nodeNr):
            coords = np.zeros(2)
            coords[0] = (nodeB[0] - nodeA[0])*indnt/length +(nodeB[0] - nodeA[0])*mD/length*i  + nodeA[0]
            coords[1] = (nodeB[1] - nodeA[1])*indnt/length +(nodeB[1] - nodeA[1])*mD/length*i  + nodeA[1]
            node_coords.append(coords)
