# distutils: language=c++
cimport cython
import numpy as np
cimport numpy as np
from libc.stdio cimport printf
from libc.stdlib cimport rand, RAND_MAX
from libcpp.vector cimport vector
from cpprandom cimport mt19937_64, uniform_real_distribution


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesRect_cython(double[:] maxLim,
                      float minDist,
                      int dim,
                      int trials,
                      list node_coords):
    print('Generating {:d}d block segment of size: {}'.format(dim, ' / '.join('{:f}'.format(i) for i in maxLim)))
    cdef:
        int generatedPoints = 0
        cdef int p
        int tr = 0
        double distInt
        int i
        vector[float] coords
        int node_coords_input_len = len(node_coords)
        vector[float] node_coords_temp
        bint distIsGood
        mt19937_64 gen = mt19937_64()
        uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0,1.0)
    coords.push_back(0.0)
    coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(dim):
                node_coords_temp.push_back(node[d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            for i in range(dim):
                coords[i] = dist(gen) * maxLim[i]
            distIsGood = True

            for p in range(generatedPoints):
                distInt = 0
                for i in range(dim):
                    distInt += ((node_coords_temp[p*dim+i] - coords[i])*
                                (node_coords_temp[p*dim+i] - coords[i]))
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for i in range(dim):
                node_coords_temp.push_back(coords[i])
            generatedPoints += 1

    # Copy back to lists
    for i in range(node_coords_input_len, generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[i*dim], node_coords_temp[i*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[i*dim], node_coords_temp[i*dim+1], node_coords_temp[i*dim+2]])

def checkMutDistancesLoops_cython(int dim,
                                  float minDist,
                                  list currentNodes,
                                  list newNode):
    cdef:
        distIsGood = True
        int p, i
        int currentNodes_len = len(currentNodes)
        double distInt
    for p in range (currentNodes_len):
        distInt = 0
        for i in range(dim):
            distInt += ((currentNodes[p*dim+i] - newNode[i])*
                        (currentNodes[p*dim+i] - newNode[i]))
        #distInt = scipy.spatial.distance.euclidean(currentNodes[p], newNode)
        if (distInt < minDist):
            distIsGood = False
            break
    return distIsGood
