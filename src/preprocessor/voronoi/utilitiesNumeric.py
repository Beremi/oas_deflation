import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import matplotlib.pylab as plt
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix

def reorderToDiagonal ():
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


    print('original connectivity matrix')
    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(A)
    #plt.colorbar()
    plt.show()

    C = np.zeros( (node_count,node_count) )
    C = csr_matrix(A)
    order = reverse_cuthill_mckee(C, symmetric_mode=True)
    order = np.asarray(order)

    B = A[order][:,order]
    #print(B)



    print('reordered connectivity matrix')
    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(B)
    #plt.colorbar()
    plt.show()

    return order
