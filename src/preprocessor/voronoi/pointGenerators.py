import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os
import utilitiesMech
import utilitiesGeom
from mpl_toolkits.mplot3d import Axes3D
from regions import *

SHOW_PLOT = False


def randPointInRectangle(dim, maxLim, useLowBound=False):

    if not useLowBound:
        coords = np.random.random(dim)
        coords *= maxLim
    else:
        topBound = maxLim[0:dim]
        lowBound = maxLim[dim:2*dim]
        coords = lowBound + np.random.random(dim) * (topBound-lowBound)


    return coords

def randPointOnLine(dim, nodeA, nodeB):
    coords = np.zeros(dim)
    r = np.random.uniform()
    coords = (nodeB - nodeA)*r +nodeA
    return coords


# generates random points no closer to each other than minDist
# into 2d or 3d block
# maxLim: n-d array of dimensions
def generateNodesRect(maxLim, minDist, dim, trials, node_coords, useLowBound=False, topMinDist = -1, bottomMinDist=-1):
    if (dim==2):
        print('Generating 2d block segment of size: %f / %f. This may take few minutes. Do not panic. \nAlthough attempt to use the Cython solution by Vasek!!!' %(maxLim[0], maxLim[1]) )
    if (dim==3):
        print('Generating 3d block segment of size: %f / %f / %f. This may take long. Keep calm. \nAlthough attempt to use the Cython solution by Vasek!!!' %(maxLim[0], maxLim[1], maxLim[2]) )

    tr = 0
    while (tr<trials):
        tr = 0
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInRectangle(dim, maxLim,useLowBound=useLowBound)
            distIsGood = True
            #
            distIsGood = utilitiesGeom.checkMutDistancesLoops(dim, minDist, node_coords, list(coords))
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
    print('''Using Python version of generator - generateNodesRect. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')


def fuller2D(d, dmax):
    return (1.065*np.sqrt(d/dmax)-0.053*np.power(d/dmax,4)-0.012*np.power(d/dmax,6)-0.0045*np.power(d/dmax,8)-0.0025*np.power(d/dmax,10))

def fuller3D(d, dmax):
    return np.sqrt(d/dmax)

def generateParticlesOrtoSurface(nodeA, nodeB, minDiam, maxDiam, volumeRatio, dim, trials, node_coords, radii, allow_domain_overlap=False):
    gap = 0.1
    maxLim = np.zeros((3))
    flatDim = -1

    #print (nodeA)
    #print (nodeB)
    for c in range (dim):
        if (nodeA[c] == nodeB[c]):
            maxLim[c] = 0
            flatDim = c
        else:
            maxLim[c] = np.abs(nodeA[c]-nodeB[c])
    #print(maxLim)
    Volume = 0

    maxLim = np.delete(maxLim, flatDim)
    Volume = np.product(maxLim)
    #print(maxLim)
    #print (Volume)


    d = np.flipud(np.linspace(minDiam*0.5,maxDiam,30))
    if(dim==2):
        freq = fuller2D(d, maxDiam)
    elif(dim==3):
        freq = fuller3D(d, maxDiam)

    saturation = 0;
    iters = 0
    di = 0
    radius = maxDiam/2.
    while (2*radius>minDiam and iters<trials):

            point = np.zeros(dim)
            for c in range (dim):
                if (nodeA[c] == nodeB[c]):
                    point[c] = nodeA[c]
                else:
                    if (allow_domain_overlap):
                        point[c] = (nodeB[c] - nodeA[c])*np.random.uniform()  + nodeA[c]
                    else:
                        point[c] = (nodeB[c] - nodeA[c] - radius*2)*np.random.uniform()  + nodeA[c]
                        #np.random.rand(dim)*(maxLim-radius*2) + radius

            #if (allow_domain_overlap): point = np.random.rand(dim)*maxLim
            #else: point = np.random.rand(dim)*(maxLim-radius*2) + radius

            approved = False
            if ( len(node_coords)==0): approved = True
            else:
                delta = np.abs(node_coords-point)

                dist2 = np.sum(np.square(delta),axis=1)
                if (min(dist2 - np.square((1.+gap)*(radii+radius)))>0): approved = True
            if ( approved) :
                node_coords = np.vstack((node_coords, point));
                radii = np.hstack((radii, radius));
                iters = 0

                saturation += np.pi*np.power(radius,2)/Volume


                sys.stdout.write('\r'+'Particles:' +  str(len(node_coords)))
                sys.stdout.flush()


                while ((1. - saturation / volumeRatio ) < freq[di]):
                     di+=1;
                     #print (' di %d' %di)

                radius = (d[di] +(d[di - 1] - d[di]) / (freq[di - 1] - freq[di])*(1. - saturation / volumeRatio - freq[di])) / 2.;
            else: iters += 1

#    print("Saturation of the volume: ", saturation)

    return node_coords, radii


def generateParticlesRect(maxLim, minDiam, maxDiam, volumeRatio, dim, trials, node_coords, radii, allow_domain_overlap = False, periodic_distance=False):
        gap = 0.1
        Volume = np.prod(maxLim)
        d = np.flipud(np.linspace(minDiam*0.5,maxDiam,30))
        if(dim==2):
            freq = fuller2D(d, maxDiam)
        elif(dim==3):
            freq = fuller3D(d, maxDiam)

        #add regular boundary
        """
        point = np.zeros(dim)
        n0 = int(2.*maxLim[0]/(maxDiam+minDiam))
        n1 = int(2.*maxLim[1]/(maxDiam+minDiam))
        s0 = maxLim[0]/n0;
        s1 = maxLim[1]/n1;
        for i in range(n0+1):
            point[0] = i*s0
            point[1] = 0
            node_coords = np.vstack((node_coords, point));
            radii = np.hstack((radii, 0));
            point[1] = maxLim[1]
            node_coords = np.vstack((node_coords, point));
            radii = np.hstack((radii, 0));
        for i in range(n0-1):
            point[1] = (i+1)*s1
            point[0] = 0
            node_coords = np.vstack((node_coords, point));
            radii = np.hstack((radii, 0));
            point[0] = maxLim[0]
            node_coords = np.vstack((node_coords, point));
            radii = np.hstack((radii, 0));

        print(node_coords)
        """

        #option to add external nodes
        """
        print("LOADING EXTERNAL NODES")
        node_coords = np.loadtxt("repnodes.inp")
        radii = node_coords[:,-1]
        node_coords = node_coords[:,0:-1]
        """

        saturation = 0;
        iters = 0
        di = 0
        radius = maxDiam/2.
        while (2*radius>minDiam and iters<trials):
                if (allow_domain_overlap): point = np.random.rand(dim)*maxLim
                else: point = np.random.rand(dim)*(maxLim-radius*2) + radius

                approved = False
                if ( len(node_coords)==0): approved = True
                else:
                    delta = np.abs(node_coords-point)
                    if (periodic_distance):
                        dist2 = np.sum(np.square(np.minimum(delta,maxLim-delta)),axis=1)
                    else:
                        dist2 = np.sum(np.square(delta),axis=1)
                    if (min(dist2 - np.square((1.+gap)*(radii+radius)))>0): approved = True
                if ( approved) :
                    node_coords = np.vstack((node_coords, point));
                    radii = np.hstack((radii, radius));
                    iters = 0
                    if (dim==2):  saturation += np.pi*np.power(radius,2)/Volume
                    elif (dim==3): saturation += np.pi/6*np.power(radius*2.,3)/Volume

                    while ((1. - saturation / volumeRatio ) < freq[di]):
                        di+=1;
                        #print (' di %d' %di)
                    radius = (d[di] +(d[di - 1] - d[di]) / (freq[di - 1] - freq[di])*(1. - saturation / volumeRatio - freq[di])) / 2.;

                    sys.stdout.write('\r'+'Particles:' +  str(len(node_coords)))

                    sys.stdout.flush()

                else: iters += 1

        print("Saturation of the volume: ", saturation)

        """
        print(node_coords.shape, radii.shape)
        node_coords = np.zeros((400,2))
        radii = np.zeros(400)
        for i in range(20):
            for j in range(20):
                node_coords[20*i+j,0] = (i+0.5)/20+np.random.rand()/1000.
                node_coords[20*i+j,1] = (j+0.5)/20+np.random.rand()/1000.
        print(node_coords.shape, radii.shape)
        """

        return node_coords, radii

try:
    from point_generators_cython import generateParticlesRect_cython as generateParticlesRect
    print('Using Cython version of point generator - generateParticlesRect.')
except:
    print('''Using Python version of generator - generateParticlesRect. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')


def generateParticlesDam(maxLim, topsize, minDiam, maxDiam, volumeRatio, dim, trials, node_coords, radii):
        gap = 0.1
        Volume = np.prod(maxLim)- maxLim[1]*maxLim[2]*(maxLim[0]-topsize)/2.
        d = np.flipud(np.linspace(minDiam,maxDiam,20))  #20 different diameters
        prob = volumeRatio*fuller2D(d, maxDiam)
        num = ((prob[:-1]-prob[1:])*Volume/(np.pi*np.square(d[:-1]/2.))+1.).astype(int)

        #option to add external nodes
        """
        node_coords = np.loadtxt("repnodes.inp")
        radii = node_coords[:,-1]
        node_coords = node_coords[:,0:-1]
        """

        alpha = np.arctan( (maxLim[0] - topsize)/maxLim[2] )
        planenorm = np.array([np.cos(alpha), 0., np.sin(alpha)])
        planeconst = -planenorm[0]*maxLim[0] - planenorm[1]*maxLim[1]

        iters = 0
        di = 0
        numi= 0
        while (d[di]>minDiam and iters<trials):
            if numi<num[di]:
                point = np.random.rand(dim)*(maxLim-d[di]) + d[di]/2.
                radius = d[di]/2.
                if (np.dot(point,planenorm) + planeconst > -d[di]/2.):
                    continue #violation of skewed boundary
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

def generateNodesRectPeriodic(maxLim, minDist, dim, trials, node_coords):
    raise NotImplementedError('''Not implemented in Python, only Cython version available. To use the Cython version
          the code has to be build using: python setup.py build_ext --inplace.''')
try:
    from point_generators_cython import generateNodesRectPeriodic_cython as generateNodesRectPeriodic
    print('Using Cython version of point generator - generateNodesRectPeriodic.')
except:
    print('''Using Python version of generator - generateNodesRectPeriodic. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')

#generates random points onto a set 3d line. No closer than minDst
#catch corners: samples the boundary points first
def generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners, equidist):
    print('Generating 3d line segment from [%f; %f; %f] to [%f; %f; %f] '
     %(nodeA[0], nodeA[1],nodeA[2],nodeB[0], nodeB[1],nodeB[2]) )

    if(catchCorners):
        node_coords.append(np.copy(nodeA))
        node_coords.append(np.copy(nodeB))

    if not equidist:
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

    else:
        #print('Equid')
        mD = minDist * 1.6
        length = np.linalg.norm(nodeA - nodeB)
        nodeNr = int (length / mD)
        indnt = (length-(nodeNr-1)*mD) / 2
        for i in range (nodeNr):
            coords = np.zeros(3)
            coords[0] = (nodeB[0] - nodeA[0])*indnt/length +(nodeB[0] - nodeA[0])*mD/length*i  + nodeA[0]
            coords[1] = (nodeB[1] - nodeA[1])*indnt/length +(nodeB[1] - nodeA[1])*mD/length*i  + nodeA[1]
            coords[2] = (nodeB[2] - nodeA[2])*indnt/length +(nodeB[2] - nodeA[2])*mD/length*i  + nodeA[2]
            node_coords.append(coords)


def generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, minDistAmongNewPoints=False):
    print('Generating 3d surface segment from [%f; %f; %f] to [%f; %f; %f]'
     %(nodeA[0], nodeA[1],nodeA[2],nodeB[0], nodeB[1],nodeB[2]) )

    if (minDistAmongNewPoints == True):
        new_points = []
        new_points.append(nodeA)

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
            if (minDistAmongNewPoints == True):
                distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, new_points, coords)
            else:
                distIsGood = utilitiesGeom.checkMutDistancesCdist(dim, minDist, node_coords, coords)

            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            if (minDistAmongNewPoints == True): new_points.append(coords)
            node_coords.append(coords)


def ortho_grid(n, nvar):
    nsim = n ** nvar
    x = np.linspace(0.5 / n, 1 - 0.5 / n, n)
    x_list = [x] * (nvar)
    X = np.meshgrid(*x_list) # format later used in RBF
    ortho_grid = np.array(X).reshape((nvar, nsim)).T

    return X, ortho_grid


def generateOrtogrid(maxLim, minDist, dim, node_coords, size):
    print ('Generating orthogonal grid...', end='')
    if (dim == 2):

        n= int (size/minDist)
        if (n%2>0):
            n+=1


        X, orthogrid = ortho_grid(n, dim)

        xmin = maxLim[0]/2 - (n/2+0.5) * minDist  #maxLim[2]
        xmax = maxLim[0]/2 + (n/2+0.5) * minDist  #maxLim[0]
        ymin = maxLim[3]
        ymax = maxLim[1]

        orthogrid [:,0] *= (xmax-xmin)
        orthogrid [:,1] *= (ymax-ymin)

        orthogrid [:,0] += xmin
        orthogrid [:,1] += ymin - 0.5/n*(ymax-ymin)

    if (dim == 3):

        n= int (size/minDist)
        if (n%2>0):
            n+=1

        X, orthogrid = ortho_grid(n, dim)

        xmin = maxLim[0]/2 - (n/2+0.5) * minDist
        xmax = maxLim[0]/2 + (n/2+0.5) * minDist
        ymin = maxLim[4]
        ymax = maxLim[1]
        zmin = maxLim[5]
        zmax = maxLim[2]

        orthogrid [:,0] *= (xmax-xmin)
        orthogrid [:,1] *= (ymax-ymin)
        orthogrid [:,2] *= (zmax-zmin)

        orthogrid [:,0] += xmin
        orthogrid [:,1] += ymin - 0.5/n*(ymax-ymin)
        orthogrid [:,2] += zmin

        """
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(orthogrid[:,0],orthogrid[:,1],orthogrid[:,2])
        if SHOW_PLOT:
            plt.show()
        #
        #"""

    print('done.')

    for i in range(len(orthogrid)):
        rnd = np.random.uniform(-1e-7, 1e-7)

        node_coords.append(orthogrid[i,:])


def ortho_grid_variable(ns):
     if isinstance(ns, int):
         ns = [ns]
     nvar = len(ns)
     nsim = np.prod(ns)
     x_list = [np.linspace(0.5 / n, 1 - 0.5 / n, n) for n in ns]
     X = np.meshgrid(*x_list) # format later used in RBF
     ortho_grid = np.array(X).reshape((nvar, nsim)).T
     return X, ortho_grid
#X, og = ortho_grid_variable((3, 4))


def generateOrtogrid_variable(maxLim, minDist, node_coords, dimensions):
    print ('Generating orthogonal grid...', end='')
    """
    if (len(dimensions) == 2):

        n= int (size/minDist)
        if (n%2>0):
            n+=1

        X, orthogrid = ortho_grid(n, dim)

        xmin = maxLim[0]/2 - (n/2+0.5) * minDist  #maxLim[2]
        xmax = maxLim[0]/2 + (n/2+0.5) * minDist  #maxLim[0]
        ymin = maxLim[3]
        ymax = maxLim[1]

        orthogrid [:,0] *= (xmax-xmin)
        orthogrid [:,1] *= (ymax-ymin)

        orthogrid [:,0] += xmin
        orthogrid [:,1] += ymin - 0.5/n*(ymax-ymin)
    """
    if (len(dimensions) == 3):
        dimensions += minDist
        n= dimensions/minDist
        ints = []
        for i in range (3):
            ints.append (int (n[i]))

        if (ints[0]%2>0):
            ints[0]+=1

        ints[1] +=1
        ints[2] +=1



        X, orthogrid = ortho_grid_variable((ints[0],ints[1],ints[2]))

        xmin = (maxLim[0]+maxLim[3])/2 - (n[0]/2) * minDist
        #xmax = (maxLim[0]+maxLim[3])/2 + (n[0]/2) * minDist
        ymin = maxLim[4] - minDist /2
        ymax = maxLim[1]
        zmin = maxLim[5] - minDist /2
        zmax = maxLim[2]
        """
        orthogrid [:,0] *= (xmax-xmin)
        orthogrid [:,1] *= (ymax-ymin)
        orthogrid [:,2] *= (zmax-zmin)
        """

        orthogrid *= dimensions

        orthogrid [:,0] += xmin
        orthogrid [:,1] += ymin #- 0.5/n[1]*(ymax-ymin)
        orthogrid [:,2] += zmin

        """
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(orthogrid[:,0],orthogrid[:,1],orthogrid[:,2])
        if SHOW_PLOT:
            plt.show()
        #
        #"""

    print('done.')

    for i in range(len(orthogrid)):
        rnd = np.random.uniform(-1e-7, 1e-7)

        node_coords.append(orthogrid[i,:])


def generateNodesOrtoCircle2dRand(center, radius, minDist, node_coords, trials):
    print ('Generating a 3d circle surface. Ctr [%f, %f], Rad: %f' %(center[0],center[1], radius))

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = randPointInCircle(center, radius, 0)
            coords = np.copy(coords[0:2])
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(2, minDist, node_coords, coords)
            #
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

try:
    from point_generators_cython import generateNodesOrtoTube3dRand_cython as generateNodesOrtoTube3dRand
    print('Using Cython version of point generator - generateNodesOrtoTube3dRand.')
except:
    print('''Using Python version of generator - generateNodesOrtoTube3dRand. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')



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

    point = np.zeros(len(center))

    rn = np.random.uniform()

    if (directionDim == 0 ):
        if len(center) == 2:
            point[0] = radius * np.cos(angle) * rn
            point[1] = radius * np.sin(angle) * rn
        if len(center) == 3:
            point[1] = radius * np.cos(angle) * rn
            point[2] = radius * np.sin(angle) * rn
    """
    if (directionDim == 1):
        point[0] = radius * np.cos(angle) * rn
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 2):
        point[0] = radius * np.cos(angle) * rn
        point[1] = radius * np.sin(angle) * rn
    """

    point += center
    return point[0:len(center)]

def randPointInAnnulus(center, radius, thickness, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)

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

    point += center
    return point

def randPointInSphere(center, radius):
    # works also for circle in xy plane (based on len(center))
    point = np.zeros(len(center))

    angle1 = np.random.uniform() * np.pi * 2
    rn = np.random.uniform()
    # print(rn, rn < 0.2 and "lower than <<<<<<<<<<<<<<<<<<<<<<" or "" )

    if len(point) > 2:
        # 3D case
        angle2 = np.random.uniform() * np.pi
        point[0] = radius * np.cos(angle1) * np.sin(angle2) * rn
        point[1] = radius * np.sin(angle1) * np.sin(angle2) * rn
        point[2] = radius * np.cos(angle2) * rn

    else:
        point[0] = radius * np.cos(angle1) * rn
        point[1] = radius * np.sin(angle1) * rn

    point += center
    return point

def randPointBetweenSpheres(center, radius, thickness):
    # works also both anulus in xy plane
    point = np.zeros(len(center))

    angle1 = np.random.uniform() * np.pi * 2

    effRadius = (radius - thickness) + thickness *  np.random.uniform()

    if len(point) > 2:
        # 3D case
        angle2 = np.random.uniform() * np.pi
        point[0] = effRadius * np.cos(angle1) * np.sin(angle2)
        point[1] = effRadius * np.sin(angle1) * np.sin(angle2)
        point[2] = effRadius * np.cos(angle2)

    else:
        point[0] = effRadius * np.cos(angle1)
        point[1] = effRadius * np.sin(angle1)

    point += center
    return point


def randPointInTube(center, radius, height, thickness, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)

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

    point += center
    return point


def generateNodesCircle2dRand(center, radius, minDist, node_coords, trials, angleLimitA=None, angleLimitB = None, mirrorIndent = 0, radiusSpread = 0):
    #print ('Generating a 2d circle border. Ctr [%s, %s], Rad: %s, Angle limit +-%s' %(center[0],center[1], radius, np.degrees(angleLimitA-angleLimitB)))

    mirroredPoints = []
    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            if (mirrorIndent >0):
                coords, mirrored_coords = randPointOnCircle(center, radius, 2, angleLimitA = angleLimitA, angleLimitB = angleLimitB,mirrorIndent = mirrorIndent)
            else:
                coords = randPointOnCircle(center, radius, 2, angleLimitA = angleLimitA, angleLimitB = angleLimitB,mirrorIndent = mirrorIndent, radiusSpread=radiusSpread)
            #
            distIsGood = utilitiesGeom.checkMutDistancesCdist(2, minDist, node_coords, coords)
            #
            if (distIsGood == False):
                tr += 1
            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)
            if (mirrorIndent >0):
                mirroredPoints.append(mirrored_coords)

    if (mirrorIndent >0):
        return mirroredPoints

def randPointOnCircle(center, radius, directionDim, angleLimitA = None, angleLimitB = None ,mirrorIndent = 0, radiusSpread=0):
    if (angleLimitA == None): angle = np.random.uniform() * np.pi * 2
    else:   angle = np.random.uniform(low=angleLimitA, high=angleLimitB) #* np.pi * 2

    if mirrorIndent > 0:
        radiusSpread = 0

    if (len(center) == 3):
        point = np.zeros(3)
        point += center

        samplingRadius = radius + np.random.uniform(low=-radiusSpread, high = radiusSpread)

        if (directionDim == 0 ):
            point[1] = samplingRadius * np.cos(angle)
            point[2] = samplingRadius * np.sin(angle)
        if (directionDim == 1):
            point[0] = samplingRadius * np.cos(angle)
            point[2] = samplingRadius * np.sin(angle)
        if (directionDim == 2):
            point[0] = samplingRadius * np.cos(angle)
            point[1] = samplingRadius * np.sin(angle)
        return point

    if (len(center) == 2):
        point = np.zeros(2)
        point += center
        mirroredPoint = np.zeros(2)
        mirroredPoint += center

        samplingRadius = radius + np.random.uniform(low=-radiusSpread, high = radiusSpread)
        #print(samplingRadius)
        point[0] += samplingRadius * np.cos(angle)
        point[1] += samplingRadius * np.sin(angle)

        if (mirrorIndent>0):
            mirroredPoint[0] += (radius-2*mirrorIndent) * np.cos(angle)
            mirroredPoint[1] += (radius-2*mirrorIndent) * np.sin(angle)
            return point,mirroredPoint
        else:
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

try:
    from point_generators_cython import generateNodesOrtoCilinder3dRand_cython as generateNodesOrtoCilinder3dRand
    print('Using Cython version of point generator - generateNodesOrtoCilinder3dRand.')
except:
    print('''Using Python version of generator - generateNodesOrtoCilinder3dRand. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')



def generateNodesOrtoCilinderSurf3dRand(center, radius, height, directionDim, minDist, node_coords, trials, angleLimitA=None, angleLimitB=None, mirrorIndent=None):
    print ('Generating a 3d cylinder surf segment. Ctr [%f, %f, %f], Rad: %f' %(center[0],center[1],center[2], radius))

    mirroredPoints = []

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            if (mirrorIndent == None):
                coords = randPointOnCilinder(center, radius, height, directionDim, angleLimitA = angleLimitA, angleLimitB = angleLimitB ,mirrorIndent = mirrorIndent)
            else:
                coords, mirrored_coords = randPointOnCilinder(center, radius, height, directionDim, angleLimitA = angleLimitA, angleLimitB = angleLimitB ,mirrorIndent = mirrorIndent)

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
            if (mirrorIndent != None):
                mirroredPoints.append(mirrored_coords)

    if (mirrorIndent != None):
        return mirroredPoints

try:
    from point_generators_cython import generateNodesOrtoCylinderSurf3dRand_cython as generateNodesOrtoCilinderSurf3dRand
    print('Using Cython version of point generator - generateNodesOrtoCilinderSurf3dRand.')
except:
    print('''Using Python version of generator - generateNodesOrtoCilinderSurf3dRand. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def distance(a, b):
    dst = 0;
    for coord_a, coord_b in zip(a, b):
        dst += (coord_a - coord_b) ** 2.0
    return dst ** 0.5


def minDistTrans(lminR, lmin, dst, rR, rT):
    return lminR + (lmin - lminR) * ( dst - rR ) / (rT - rR)



def generateNodesRemesh(node_coords, trials, maxLim, minDistRemesh, minDist,
                        centersToRemesh, centersPreviouslyRemeshed,
                        radiusRemesh, radiusTransitional,
                        dim, useExistingFineNodes=False,
                        remesherSeed=1):
    np.random.seed(seed=remesherSeed)  # same seed for remesher all the time to be able to replicate the adaptive remesh (test version 26.5.2021)
    PRINT_TEST = False
    print ( 'Generating points to update geometry' )
    if dim == 2:
        border_block = Block(Point(0., 0.),
                             Point(maxLim[0], maxLim[1]))
    elif dim == 3:
        border_block = Block(Point(0., 0., 0),
                             Point(maxLim[0], maxLim[1], maxLim[2]))
    tr = 0
    # remesh fine areas
    ci = 0
    for center in centersToRemesh:
        if not useExistingFineNodes:
            tr = 0
            # print("generating in remesh area %d - " % (ci), end=' ' )
            # print(center)
            ci += 1
            while (tr<trials):
                tr = 0
                distIsGood = False
                while (distIsGood == False):
                    if PRINT_TEST:
                        print("generating node %d, trials: %d/%d" % (len(node_coords), tr, trials), end='\r')
                    distIsGood = True
                    if (tr > trials): break

                    coords = randPointInSphere(center, radiusRemesh)
                    p_coord = (dim == 2) and Point(coords[0], coords[1]) or Point(coords[0], coords[1], coords[2])

                    # check if generated point is not outside the specimen
                    # if rectLims is not None:
                    if not border_block.IsInside(p_coord):
                        distIsGood = False
                        tr += 1
                        continue
                    # check distances to already remeshed centers (if it is not in any previously remeshed circle/sphere)
                    distIsGood = utilitiesGeom.checkMutDistancesLoops(dim,
                                              radiusRemesh,
                                              centersPreviouslyRemeshed,
                                              list(coords))
                    if not distIsGood:
                        tr += 1
                        continue

                    # check distances to other nodes
                    distIsGood = utilitiesGeom.checkMutDistancesLoops(dim,
                                                 minDistRemesh,
                                                 node_coords, list(coords))
                    if not distIsGood:
                        tr +=1
                        continue

                #Adding node coords
                if (tr < trials):
                    if PRINT_TEST:
                        print("appending node in remesh area ---------------------", end=' ' )
                        print(coords)
                    node_coords.append(coords)

        # after remesh of fine area, append this one into centers previously remeshed
        centersPreviouslyRemeshed.append(center)
    ##########################################################################
    # print(centersPreviouslyRemeshed)
    # remesh the transitional areas - same centers, just different radius and only outer ring
    ci = 0
    tr = 0
    for center in centersToRemesh:
        # if PRINT_TEST:
        #     print("generating in transitional area %d" % ci, end='\r' )
        tr = 0
        ci += 1
        while (tr<trials):
            tr = 0
            distIsGood = False
            while (distIsGood == False):
                if PRINT_TEST:
                    print("generating tansitional node %d, trials: %d/%d" % (len(node_coords), tr, trials), end='\r')
                if (tr > trials): break
                distIsGood = True

                coords = randPointBetweenSpheres(center, radiusTransitional,
                            thickness=(radiusTransitional-radiusRemesh))
                p_coord = (dim == 2) and Point(coords[0], coords[1]) or Point(coords[0], coords[1], coords[2])

                # check if generated point is not outside the specimen
                # if rectLims is not None:
                if not border_block.IsInside(p_coord):
                    distIsGood = False
                    tr += 1
                    continue

                # check distances to already remeshed centers - to prevent putting nodes there centers
                distIsGood = utilitiesGeom.checkMutDistancesLoops(dim,
                                        radiusRemesh,
                                        centersPreviouslyRemeshed,
                                        list(coords))
                if not distIsGood:
                    tr += 1
                    continue


                # TODO tady neco nesedí: viz obrázek co to exportuje
                mdt = minDistTrans(minDistRemesh, minDist, distance(coords,
                                                                    center),
                             radiusRemesh, radiusTransitional)
                # mdt = (minDistRemesh + (minDist - minDistRemesh) *
                #                  ( distance(coords, center) - radiusRemesh ) /
                #                  (radiusTransitional - radiusRemesh)
                #                 )

                # print(minDist, minDistRemesh, mdt)
                # print(minDistTrans(minDistRemesh, minDist, radiusTransitional,
                #              radiusRemesh, radiusTransitional))

                # check distances to other nodes
                distIsGood = utilitiesGeom.checkMutDistancesLoops(dim,
                                            mdt,
                                            node_coords, list(coords))
                if not distIsGood:
                    tr += 1
                    continue


            #Adding node coords
            if (tr < trials):
                if PRINT_TEST:
                    print("appending node in transitional area ----------", end=' ' )
                    print(coords)
                node_coords.append(coords)

    return node_coords


try:
    from point_generators_cython import generateNodesRemesh_cython as generateNodesRemesh
    print('Using Cython version of point generator - generateNodesRemesh.')
except Exception as e:
    print(e)
    print('''Cython version of point generator - generateNodesRemesh - not avaliable yet.''')
    # print('''Using Python version of generator. To use the Cython version the
    #       the code has to be build using: python setup.py build_ext --inplace.''')

###############################################################################
###############################################################################
###############################################################################



def randPointInCilinder(center, radius, height, directionDim):
    angle = np.random.uniform() * np.pi * 2

    point = np.zeros(3)
    point += center

    rn = np.random.uniform()

    if (directionDim == 0 ):
        point[0] = height * np.random.uniform()
        point[1] = radius * np.cos(angle) * rn
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 1):
        point[0] = radius * np.cos(angle) * rn
        point[1] = height * np.random.uniform()
        point[2] = radius * np.sin(angle) * rn
    if (directionDim == 2):
        point[0] = radius * np.cos(angle) * rn
        point[1] = radius * np.sin(angle) * rn
        point[2] = height * np.random.uniform()

    return point

def randPointOnCilinder(center, radius, height, directionDim, angleLimitA = None, angleLimitB = None ,mirrorIndent = None):
    #angle = np.random.uniform() * np.pi * 2
    if (angleLimitA == None): angle = np.random.uniform() * np.pi * 2
    else:   angle = np.random.uniform(low=angleLimitA, high=angleLimitB)

    point = np.zeros(3)
    point += center
    mirroredPoint = np.zeros(3)
    mirroredPoint += center

    if (directionDim == 0 ):
        point[0] += height * np.random.uniform()
        point[1] += radius * np.cos(angle)
        point[2] += radius * np.sin(angle)
        if (mirrorIndent!=None):
            mirroredPoint[1] += (radius-2*mirrorIndent) * np.cos(angle)
            mirroredPoint[2] += (radius-2*mirrorIndent) * np.sin(angle)
            mirroredPoint[0] += height * np.random.uniform()
    if (directionDim == 1):
        point[0] += radius * np.cos(angle)
        point[1] += height * np.random.uniform()
        point[2] += radius * np.sin(angle)
        if (mirrorIndent!=None):
            mirroredPoint[0] += (radius-2*mirrorIndent) * np.cos(angle)
            mirroredPoint[2] += (radius-2*mirrorIndent) * np.sin(angle)
            mirroredPoint[1] += height * np.random.uniform()
    if (directionDim == 2):
        point[0] += radius * np.cos(angle)
        point[1] += radius * np.sin(angle)
        point[2] += height * np.random.uniform()
        if (mirrorIndent!=None):
            mirroredPoint[0] += (radius-2*mirrorIndent) * np.cos(angle)
            mirroredPoint[1] += (radius-2*mirrorIndent) * np.sin(angle)
            mirroredPoint[2] += height * np.random.uniform()

    if (mirrorIndent!=None):
        return point,mirroredPoint
    else:
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
        mD = minDist
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
