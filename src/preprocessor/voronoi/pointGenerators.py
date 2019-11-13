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
        print('Generating 2d block segment of size: %f / %f. This may take few minutes. Do not panic. \nAlthough attempt to use the Cython solution by Vasek!!!' %(maxLim[0], maxLim[1]) )
    if (dim==3):
        print('Generating 3d block segment of size: %f / %f / %f. This may take long. Keep calm. \nAlthough attempt to use the Cython solution by Vasek!!!' %(maxLim[0], maxLim[1], maxLim[2]) )

    tr = 0
    while (tr<trials):
        tr = 0
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInRectangle(dim, maxLim)
            distIsGood = True
            #
            distIsGood = utilitiesGeom.checkMutDistancesLoops(dim, minDist, node_coords, coords)
            #distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, node_coords, coords)
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


def fuller2D(d, dmax):
    return (1.065*np.sqrt(d/dmax)-0.053*np.power(d/dmax,4)-0.012*np.power(d/dmax,6)-0.0045*np.power(d/dmax,8)-0.0025*np.power(d/dmax,10))

def generateParticlesRect(maxLim, minDiam, maxDiam, volumeRatio, dim, trials, node_coords, radii):
        gap = 0.1
        Volume = np.prod(maxLim)
        d = np.flipud(np.linspace(minDiam,maxDiam,20))  #20 different diameters
        prob = volumeRatio*fuller2D(d, maxDiam)
        num = ((prob[:-1]-prob[1:])*Volume/(np.pi*np.square(d[:-1]/2.))+1.).astype(int)

        point = np.zeros(dim+1)
        iters = 0
        di = 0
        numi= 0
        while (d[di]>minDiam and iters<trials):
            if numi<num[di]:
                point = np.random.rand(2)*(maxLim-d[di]) + d[di]/2.
                radius = d[di]/2.
                if len(node_coords) == 0:                     
                    node_coords = np.vstack((node_coords,point)); 
                    radii = np.hstack((radii, radius));
                    numi += 1
                    continue     
     
                dist = min(np.sum(np.square(node_coords-point),1)-np.square((1+gap)*(radii+radius)))
                if dist>0.:
                    node_coords = np.vstack((node_coords, point));
                    radii = np.hstack((radii, radius));
                    iters = 0
                    numi += 1
                else: iters += 1 
            else:
                di += 1
                numi = 0.

        return node_coords, radii
        
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


def generateNodesOrtoCircle3dRand(center, radius, directionDim, minDist, node_coords, trials):
    print ('Generating a 3d circle surface. Ctr [%f, %f, %f], Rad: %f' %(center[0],center[1],center[2], radius))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInCircle(center, radius, directionDim)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(3, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)

def generateNodesOrtoAnnulus3dRand(center, radius, thickness, directionDim, minDist, node_coords, trials):
    print ('Generating a 3d annulus surface. Ctr [%f, %f, %f], Rad: %f, Thick: %f' %(center[0],center[1],center[2], radius, thickness))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInAnnulus(center, radius, thickness,directionDim)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(3, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)

def generateNodesOrtoTube3dRand(center, radius, height, thickness, directionDim, minDist, node_coords, trials):
    print ('Generating a 3d tube. Ctr [%f, %f, %f], Rad: %f, Thick: %f' %(center[0],center[1],center[2], radius, thickness))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInTube(center, radius, height, thickness, directionDim)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(3, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)


def generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials):
    print ('Generating a 3d circle border. Ctr [%f, %f, %f], Rad: %f' %(center[0],center[1],center[2], radius))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointOnCircle(center, radius, directionDim)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(3, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)

def randPointInCircle(center, radius, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    rn = np.random.uniform()

    if (directionDim == 0 ):
        point[1] = radius * np.cos(angle) * rn
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 1):
        point[0] = radius * np.cos(angle) * rn
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 2):
        point[0] = radius * np.cos(angle) * rn
        point[1] = radius * np.sin(angle) * rn

    return point

def randPointInAnnulus(center, radius, thickness, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    effRadius = (radius - thickness) + thickness *  np.random.uniform()

    if (directionDim == 0 ):
        point[1] =  np.cos(angle) * effRadius
        point[2] =  np.sin(angle) * effRadius
    if (directionDim == 1):
        point[0] =  np.cos(angle) * effRadius
        point[2] =  np.sin(angle) * effRadius
    if (directionDim == 2):
        point[0] =  np.cos(angle) * effRadius
        point[1] =  np.sin(angle) * effRadius

    return point

def randPointInTube(center, radius, height, thickness, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    effRadius = (radius - thickness) + thickness *  np.random.uniform()

    if (directionDim == 0 ):
        point[0] = height * np.random.uniform()
        point[1] =  np.cos(angle) * effRadius
        point[2] =  np.sin(angle) * effRadius
    if (directionDim == 1):
        point[0] =  np.cos(angle) * effRadius
        point[1] = height * np.random.uniform()
        point[2] =  np.sin(angle) * effRadius
    if (directionDim == 2):
        point[0] =  np.cos(angle) * effRadius
        point[1] =  np.sin(angle) * effRadius
        point[2] = height * np.random.uniform()

    return point

def randPointOnCircle(center, radius, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    if (directionDim == 0 ):
        point[1] = radius * np.cos(angle)
        point[2] = radius * np.sin(angle)
    if (directionDim == 1):
        point[0] = radius * np.cos(angle)
        point[2] = radius * np.sin(angle)
    if (directionDim == 2):
        point[0] = radius * np.cos(angle)
        point[1] = radius * np.sin(angle)

    return point

def generateNodesOrtoCilinder3dRand(center, radius, height, directionDim, minDist, node_coords, trials):
    print ('Generating a 3d cylinder segment. Ctr [%f, %f, %f], Rad: %f' %(center[0],center[1],center[2], radius))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInCilinder(center, radius, height, directionDim)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(3, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)

def generateNodesOrtoCilinderSurf3dRand(center, radius, height, directionDim, minDist, node_coords, trials):
    print ('Generating a 3d cylinder surf segment. Ctr [%f, %f, %f], Rad: %f' %(center[0],center[1],center[2], radius))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointOnCilinder(center, radius, height, directionDim)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(3, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)


def randPointInCilinder(center, radius, height, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    rn = np.random.uniform()

    if (directionDim == 0 ):
        point[0] = height * rn
        point[1] = radius * np.cos(angle) * rn
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 1):
        point[0] = radius * np.cos(angle) * rn
        point[1] = height * rn
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 2):
        point[0] = radius * np.cos(angle) * rn
        point[1] = radius * np.sin(angle) * rn
        point[2] = height * rn

    return point

def randPointOnCilinder(center, radius, height, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    if (directionDim == 0 ):
        point[0] = height * np.random.uniform()
        point[1] = radius * np.cos(angle)
        point[2] = radius * np.sin(angle)
    if (directionDim == 1):
        point[0] = radius * np.cos(angle)
        point[1] = height * np.random.uniform()
        point[2] = radius * np.sin(angle)
    if (directionDim == 2):
        point[0] = radius * np.cos(angle)
        point[1] = radius * np.sin(angle)
        point[2] = height * np.random.uniform()

    return point

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


def fuller2D(d, dmax):
    return (1.065 * np.sqrt(d / dmax) - 0.053 * np.power(d / dmax, 4) -
            0.012 * np.power(d / dmax, 6) - 0.0045 * np.power(d/dmax, 8) -
            0.0025 * np.power(d / dmax, 10))


def generateTesC(lx, ly, lz=None, seed=None):

    if lz:
        dim = 3
        size = np.array([lx, ly, lz])
        V = lx*ly*lz
    else:
        V = lx*ly
        dim = 2
        size = np.array([lx, ly])
    zero = np.array([0.] * dim)

    dmax = 2.
    dmin = 0.4
    Pk = 0.75
    d = np.linspace(dmin, dmax, 20)[::-1]
    prob = Pk * fuller2D(d, dmax)
    num = ((prob[:-1] - prob[1:]) * V /
           (np.pi * np.square(d[:-1] / 2.)) + 1.).astype(int)

    nodes = np.empty((0, dim))
    radii = np.empty((0))
    beams = []

    lmin = 1.
    i = 0

    xlim = [0., lx]
    ylim = [0., ly]
    if seed:
        np.random.seed(seed)

    # NODES
    di = 0
    numi = 0
    while d[di] > dmin:
        if numi < num[di]:
            point = np.random.rand(dim)*(size-d[di]) + (zero+d[di]/2.)
            if len(nodes) == 0:
                nodes = np.vstack((nodes, point))
                radii = np.hstack((radii, d[di]/2.))
                numi += 1
                continue

            dist = min(np.sqrt(np.sum(np.square(nodes-point), 1)) -
                       1.1*(radii+d[di]/2.))
            if dist > 0.:
                nodes = np.vstack((nodes, point))
                radii = np.hstack((radii, d[di]/2.))
                numi += 1
                i = 0
            else:
                i += 1
        else:
            di += 1
            numi = 0.

    return nodes, radii
