# distutils: language=c++
cimport cython
import numpy as np
cimport numpy as np
import random
from libc.stdio cimport printf
from libc.stdlib cimport rand, RAND_MAX
from libcpp.vector cimport vector
#from cpprandom cimport mt19937_64, uniform_real_distribution
from libcpp cimport bool
from libc.math cimport sin, cos, abs

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesRect_cython(double[:] maxLim,
                      float minDist,
                      int dim,
                      int trials,
                      list node_coords,
                      bool useLowBound=False):
    print('Generating {:d}d block segment of size: {}'.format(dim, ' / '.join('{:f}'.format(i) for i in maxLim)))
    cdef:
        int generatedPoints = 0
        int p, d, i = 0
        int tr = 0
        double distInt
        vector[double] coords
        int node_coords_input_len = 0
        vector[double] node_coords_temp
        double[:] topBound, lowBound
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(dim):
        coords.push_back(0.0)

    if node_coords:
        node_coords_input_len = len(node_coords)
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i][d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            if (useLowBound==True):
                topBound = maxLim[0:dim]
                lowBound = maxLim[dim:2*dim]
                for d in range(dim):
                    coords[d] = lowBound[d] + np.random.rand() * (topBound[d] - lowBound[d])

            else:
                for d in range(dim):
                    coords[d] = np.random.rand() * maxLim[d]

            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(dim):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(dim):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])










@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoCilinder3dRand_cython(
                      double[:] center,
                      double radius,
                      double height,
                      int directionDim,
                      float minDist,
                      list node_coords,
                      int trials
                      ):
    print ('Generating a 3d cylinder segment. ')
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        int dim = 3
        double distInt
        double angle
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        double[:] rn
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(3):
        coords.push_back(0.0)

    if node_coords:
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i][d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            rn = np.random.rand(5)
            angle = rn[0]  * np.pi * 2
            if (directionDim == 0 ):
                coords[0] = center[0] + height * rn[1]
                coords[1] = center[1] +radius * cos(angle) * rn[4]
                coords[2] = center[2] +radius * sin(angle) * rn[4]
            if (directionDim == 1):
                coords[0] = center[0] +radius * cos(angle) * rn[4]
                coords[1] = center[1] +height * rn[2]
                coords[2] = center[2] +radius * sin(angle) * rn[4]
            if (directionDim == 2):
                coords[0] = center[0] +radius * cos(angle) * rn[4]
                coords[1] = center[1] +radius * sin(angle) * rn[4]
                coords[2] = center[2] +height * rn[3]



            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(3):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])





@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoTube3dRand_cython(
                      double[:] center,
                      double outerRad,
                      double height,
                      double thickness,
                      int directionDim,
                      float minDist,
                      list node_coords,
                      int trials
                      ):
    print ('Generating a 3d tube segment. cython ')
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        int dim = 3
        double distInt
        double angle
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(3):
        coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(3):
                node_coords_temp.push_back(node[d])

    innerRad = outerRad - thickness
    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            angle = random.random()  * np.pi * 2
            rn = random.random()
            if (directionDim == 0 ):
                coords[0] = center[0] + height * random.random()
                coords[1] = center[1] + (innerRad + thickness * rn) * np.cos(angle)
                coords[2] = center[2] + (innerRad + thickness * rn) * np.sin(angle)
            if (directionDim == 1 ):
                coords[0] = center[0] + (innerRad + thickness * rn) * np.cos(angle)
                coords[1] = center[1] + height * random.random()
                coords[2] = center[2] + (innerRad + thickness * rn) * np.sin(angle)
            if (directionDim == 2 ):
                coords[0] = center[0] + (outerRad + thickness * rn) * np.cos(angle)
                coords[1] = center[1] + (outerRad + thickness * rn) * np.sin(angle)
                coords[2] = center[2] + height * random.random()



            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(3):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoCylinderSurf3dRand_cython(
                      double[:] center,
                      double radius,
                      double height,
                      int directionDim,
                      float minDist,
                      list node_coords,
                      int trials,
                       angleLimitA=None, angleLimitB=None, mirrorIndent=None
                      ):
    print ('Generating a 3d cylinder surf cython. ')
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        int dim = 3
        double distInt
        double angle
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(3):
        coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(3):
                node_coords_temp.push_back(node[d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            angle = random.random()  * np.pi * 2
            rn = random.random()
            if (directionDim == 0 ):
                coords[0] = center[0] + height * random.random()
                coords[1] = center[1] + radius * np.cos(angle)
                coords[2] = center[2] + radius * np.sin(angle)
            if (directionDim == 1 ):
                coords[0] = center[0] + radius * np.cos(angle)
                coords[1] = center[1] + height * random.random()
                coords[2] = center[2] + radius * np.sin(angle)
            if (directionDim == 2 ):
                coords[0] = center[0] + radius * np.cos(angle)
                coords[1] = center[1] + radius * np.sin(angle)
                coords[2] = center[2] + height * random.random()



            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(3):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def checkMutDistancesLoops_cython(int dim,
                                  float minDist,
                                  list currentNodes,
                                  list newNode):
    cdef:
        bool distIsGood = True
        int p, d
        int currentNodes_len = len(currentNodes)
        double distInt, dx

    # NOTE JK: currentNodes list must be faltten, but maybe in some older versions it was already sent here as a flatt array (but then does not make sense (len() and range(p*dim + d)))
    if currentNodes_len == 0:
      return distIsGood
    if type(currentNodes[0]) is not float:
      currentNodes = np.concatenate(currentNodes).ravel().tolist()


    for p in range (currentNodes_len):
        distInt = 0
        for d in range(dim):
            dx = currentNodes[p*dim+d] - newNode[d]
            distInt += dx * dx
        distInt = distInt**0.5
        if (distInt < minDist):
            distIsGood = False
            break


    return distIsGood






@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesRectPeriodic_cython(double[:] maxLim,
                      float minDist,
                      int dim,
                      int trials,
                      list node_coords):
    print('Generating {:d}d block segment of size: {}'.format(dim, ' / '.join('{:f}'.format(i) for i in maxLim)))
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        double distInt, dx
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(dim):
        coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(dim):
                node_coords_temp.push_back(node[d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            #for d in range(dim):
            #    coords[d] = dist(gen) * maxLim[d]
            for d in range(dim):
                 coords[d] = np.random.rand() * maxLim[d]
            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(dim):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]

                    if (abs(dx+maxLim[d]) < abs(dx) ):
                        dx += maxLim[d]

                    elif (abs(dx-maxLim[d]) < abs(dx) ):
                        dx -= maxLim[d]


                    distInt += dx * dx
                distInt = distInt**0.5

                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(dim):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1
            #print (generatedPoints)

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])






@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def checkMutDistancesLoopsPeriodic_cython(int dim, float minDist, list currentNodes, list newNode, double[:] maxLim):
    cdef:
        bool distIsGood = True
        int p, d
        int currentNodes_len = len(currentNodes)
        double distInt, dx

    for p in range (currentNodes_len):
        distInt = 0
        for d in range(dim):
            dx = currentNodes[p*dim+d] - newNode[d]

            if (abs(dx+maxLim[d]) < abs(dx) ):
                dx += maxLim[d]
            elif (abs(dx-maxLim[d]) < abs(dx) ):
                dx -= maxLim[d]

            distInt += dx * dx

        distInt = distInt**0.5

        if (distInt < minDist):
            distIsGood = False
            break

    return distIsGood
