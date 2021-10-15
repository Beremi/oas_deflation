import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os
import time
from mpl_toolkits.mplot3d import Axes3D
import utilitiesMech
import Preprocessor as prepro

masterFile                  = "master.inp"
nodesFile                   = "nodes.inp"
verticesFile                = "vertices.inp"
mechElemsFile               = "mechElems.inp"
trsprtElemsFile             = "trsprtElems.inp"
mechBCFile                  = "mechBC.inp"
mechICFile                  = "mechIC.inp"
trsprtBCFile                = "trsprtBC.inp"
trsprtICFile                = "trsprtIC.inp"
materialsFile               = "materials.inp"
functionsFile               = "functions.inp"
initConditionsMechFile      = "initCondMech.inp"
initConditionsTrsprtFile    = "initCondTrsprt.inp"
auxNodesFile                = "auxNodes.inp"
exportersFile               = "exporters.inp"
blocksFile                  = "blocks.inp"
govNodesFile                = "govNodes.inp"
govNodesTrsptFile           = "govNodesTrspt.inp"
constraintFile              = "constraint.inp"
constraintTrsptFile         = "constraintTrspt.inp"
solverFile                  = "solver.inp"

SHOW_PLOT = False
AXIS_ASPECT_EQUAL = False  # True may cause error using newer matplotlib versions
#
#coplanarity test
def equation_plane(pA, pB, pC, pD):
    """
    >>> pA = np.array([0,1,0])
    >>> pB = np.array([1,1,0])
    >>> pC = np.array([0,1,4])
    >>> pD = np.array([0,0,0])
    >>> equation_plane(pA, pB, pC, pD)
    4
    """

    a1 = pB[0] - pA[0]
    b1 = pB[1] - pA[1]
    c1 = pB[2] - pA[2]
    a2 = pC[0] - pA[0]
    b2 = pC[1] - pA[1]
    c2 = pC[2] - pA[2]

    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * pA[0] - b * pA[1] - c * pA[2])

    # checking if the 4th point satisfies
    # equation of plane a*x + b*y + c*z = 0 #
    condition = a * pD[0] + b * pD[1] + c * pD[2] + d

    #print('%f' %condition )
   # if(condition == 0 ):
        #print("Coplanar")
    #else:
        #print("Not Coplanar err: %.2E" %(a * pD[0] + b * pD[1] + c * pD[2] + d ))
    #
    return condition

#return normalized vector normal to plane determined by 3 points
def getPlaneNormalVector (pA, pB, pC):
    a1 = pB[0] - pA[0]
    b1 = pB[1] - pA[1]
    c1 = pB[2] - pA[2]
    a2 = pC[0] - pA[0]
    b2 = pC[1] - pA[1]
    c2 = pC[2] - pA[2]

    normal = np.zeros(3)
    normal[0] = b1 * c2 - b2 * c1
    normal[1] = a2 * c1 - a1 * c2
    normal[2] = a1 * b2 - b1 * a2

    #print (np.linalg.norm(normal))
    normal = normal / np.linalg.norm(normal)
    #print (np.linalg.norm(normal))

    return  normal

def angleBetweenVectors (vecA, vecB):
    cos = np.dot(vecA,vecB) / np.norm(vecA) / np.norm(vecB)
    angle = np.arccos(np.clip(c, -1, 1))
    print ('Angle is %f' %(ang))


#check if any number in matrix is lower than
def checkLowerThan_old(matrix, minDist):
    '''
    >>> matrix = np.array([[.3, .5], [.1, .2]]).flatten()
    >>> minDist = .6
    >>> checkLowerThan_old (matrix, minDist)
    False
    >>> minDist = .1
    >>> checkLowerThan_old (matrix, minDist)
    True
    '''
    return  all(i >= minDist for i in matrix)

#check if any number in matrix is lower than
def checkLowerThan (matrix, minDist):
    '''
    >>> matrix = np.array([[.3, .5], [.1, .2]])
    >>> minDist = .6
    >>> checkLowerThan (matrix, minDist)
    False
    >>> minDist = .1
    >>> checkLowerThan (matrix, minDist)
    True
    '''
    return np.all(matrix >= minDist)



#check mutual distances between particles using cDist
def checkMutDistancesCdist(dim, minDist, currentNodes, newNode):
    ncrds = np.asarray(currentNodes)
    crds = np.asarray(newNode)
    crds = np.reshape(crds, (-1, dim))
    dists = scipy.spatial.distance.cdist(crds, ncrds , 'euclidean')
    #dists = dists.flatten()
    distIsGood = checkLowerThan(dists, minDist)
    return distIsGood

#check mutual distances between particles using cKDTree
def checkMutDistancesCKDTree (dim, minDist, currentNodes, newNode):
    crds = np.asarray(newNode)
    ncrds = np.asarray (currentNodes)
    tree = scipy.spatial.cKDTree ( ncrds,  leafsize=50 )
    violatingPoints = tree.query_ball_point (x = crds,  r = minDist, n_jobs = -1 )
    distIsGood = True
    if ( len(violatingPoints) != 0):
        distIsGood = False
    return distIsGood

#check mutual distances between particles using loops
def checkMutDistancesLoops (dim, minDist, currentNodes, newNode):
    distIsGood = True
    for p in range (len(currentNodes)):
        distInt = scipy.spatial.distance.euclidean(currentNodes[p], newNode)
        if (distInt < minDist):
            distIsGood = False
            break
    return distIsGood

try:
    from point_generators_cython import checkMutDistancesLoops_cython as checkMutDistancesLoops
    print('Using Cython version of point generator - checkMutDistancesLoops.')
except:
    print('''Using Python version of generator - checkMutDistancesLoops. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')



#check mutual distances between particles using loops in a periodic domain
def checkMutDistancesLoopsPeriodic (dim, minDist, currentNodes, newNode, maxLim):
    distIsGood = True
    return distIsGood

try:
    from point_generators_cython import checkMutDistancesLoopsPeriodic_cython as checkMutDistancesLoopsPeriodic
    print('Using Cython version of point generator - checkMutDistancesLoopsPeriodic.')
except:
    print('''Using Python version of generator - checkMutDistancesLoopsPeriodic. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')


def extractGeometry (master_folder, dim, node_count, maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=None, periodicModel = 0, notches = None, isTube=False, coupled=False, minDist = 0, node_indices_dogbone=[], randomizeMaterial=False):
    if (dim == 2):
        if (periodicModel == 0):
            vert_count, verticesIdxDict, vertIdxStart, totalNodeCount = output2D(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=mZ, notches = notches, coupled=coupled, node_indices_dogbone=node_indices_dogbone, randomizeMaterial=randomizeMaterial)
        if (periodicModel == 1):
            vert_count, verticesIdxDict, vertIdxStart, totalNodeCount = output2DPeriodic(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, minDist, mZ=mZ)
    if (dim == 3):
        if (periodicModel == 0):
            vert_count, verticesIdxDict, vertIdxStart,totalNodeCount = output3D(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=mZ,  notches = notches, isTube=isTube, coupled=coupled, randomizeMaterial=randomizeMaterial)
        if (periodicModel == 1):
            vert_count, verticesIdxDict, vertIdxStart,totalNodeCount = output3Dperiodic(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, minDist, mZ=mZ,  notches = notches, isTube=isTube)
    return vert_count, verticesIdxDict, vertIdxStart, totalNodeCount


#Extract geometry 2d
def output2D(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=None, notches = None, coupled=False, node_indices_dogbone=[], randomizeMaterial=False):
    dim = 2
    print('Extracting the geometry...', end='')
    sys.stdout.flush()
    nodes_out = np.zeros( (node_count, (2 + 1 + 1 )))
    nodes_out[:, 0:2] = vor.points[0:node_count , 0:2]
    nodes_out[:, dim] = 0
    nodes_out[:, dim + 1] = 0

    #relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
    #print ('Area error: %.5e ' %(relAreaError)  )
    ########################################################################################

    node_indices_dogbone = np.asarray(node_indices_dogbone)
    #print((node_indices_dogbone))
    if len(node_indices_dogbone) > 0:
        #cond = np.any((vor.ridge_points[:,:,None] == node_indices_dogbone), axis=2)
        #cond = np.all(cond, axis=1)
        #validRidgeIdxs = np.where(cond)[0]
        validRidgeIdxs = []
        for i in range (vor.ridge_points.shape[0]):
            pr = False
            for p in range (2):
                if (vor.ridge_points[i][p] in node_indices_dogbone):
                    pr=True
            if (pr):
               validRidgeIdxs.append(i)
        validRidgeIdxs = np.asarray(validRidgeIdxs)
    else:
        cond = np.any((vor.ridge_points < node_count) & (vor.ridge_points >= 0), axis=1)
        validRidgeIdxs = np.where(cond)[0]
    #print(validRidgeIdxs.shape)
#    print(validRidgeIdxs)

    #REMOVE
    #validRidgeIdxs = []
    #for i in range (vor.ridge_points.shape[0]):
        #pr = False
        #for p in range (2):
            #if (vor.ridge_points[i][p] < node_count):
                #pr=True

        #if (pr):
           #validRidgeIdxs.append(i)

    #validRidgeIdxs = np.asarray(validRidgeIdxs)
    ########################################################################################
    # vertices: [xA,yA,zA] [origIdx]
    vertices_out = []
    # dictionary of original and new indices of vertices
    verticesIdxDict = {}
    # ridges: nodeAidx, nodeBidx, trsprtBC, vertIdx
    ridges_out = []
    #auxiliary nodes
    aux_nodes = []

    #vertices
    ####################################################
    #list of vertices, list of beams
    print()
    start = time.time()
    vertices_out_set = set()
    for i in range (validRidgeIdxs.size):
        sys.stdout.write('\r'+'Ridge nr. ' + str(i) + '/' +  str(validRidgeIdxs.size)+'  '+          str(int(i/validRidgeIdxs.size*100))+'%')

        sys.stdout.flush()


        #array for two vertices A and B
        vrtxA = np.zeros ( (dim  +1 +1 ) )
        vrtxB = np.zeros ( (dim  +1 +1 ) )

        #original indices of vertices A and B
        vertA = vor.ridge_vertices[validRidgeIdxs[i]][0]
        vertB = vor.ridge_vertices[validRidgeIdxs[i]][1]

        # copying of coordinates of vertices A and B
        vrtxA [:dim] = vor.vertices[vertA]
        vrtxB [:dim] = vor.vertices[vertB]

        #copying of original indices of vertices A and B
        vrtxA[dim] = vertA
        vrtxB[dim] = vertB

        vrtxA[dim+1] = 0
        vrtxB[dim+1] = 0

        #duplicity check
        addVrtxA = True
        addVrtxB = True

        # OLD - SLOW
        #for j in range (len(vertices_out)):
        #    if (vertices_out[j][0] ==  vor.vertices[vertA][0] and vertices_out[j][1] ==  vor.vertices[vertA][1]):
        #        addVrtxA = False
        #    if (vertices_out[j][0] ==  vor.vertices[vertB][0] and vertices_out[j][1] ==  vor.vertices[vertB][1]):
        #        addVrtxB = False

        ## SLOWER ALTERNATIVE
        #vertices_out_arr = np.array(vertices_out)
        #if vertices_out_arr.shape[0] > 0:
        #    vertices_out_arr = vertices_out_arr[:, :dim]
        #    cond = np.all(vertices_out_arr == vor.vertices[vertA], axis=1)
        #    addVrtxA = ~np.any(cond)
        #    cond = np.all(vertices_out_arr == vor.vertices[vertB], axis=1)
        #    addVrtxB = ~np.any(cond)
        addVrtxA = tuple(vor.vertices[vertA].tolist()) not in vertices_out_set
        addVrtxB = tuple(vor.vertices[vertB].tolist()) not in vertices_out_set

        #adding the vertices into the list of vertices if new
        if (addVrtxA == True):
            verticesIdxDict.update( { vertA : len(vertices_out)  } )
            vrtxA [dim] = len(vertices_out)
            vertices_out.append(vrtxA)
            vertices_out_set.add(tuple(vor.vertices[vertA].tolist()))

        if (addVrtxB == True):
            verticesIdxDict.update( { vertB : len(vertices_out)  } )
            vrtxB [dim] = len(vertices_out)
            vertices_out.append(vrtxB)
            vertices_out_set.add(tuple(vor.vertices[vertB].tolist()))

        #ridges
        ########################################################
        #Array for ridge nAidx, nBidx, nrVrt, vertAidx, vertBidx
        rdg = np.zeros ( (2 + 1 +  2) )

        #indices of two nodes that are divided by the ridge
        pointA = vor.ridge_points[validRidgeIdxs[i],0]
        pointB = vor.ridge_points[validRidgeIdxs[i],1]

        #creating auxiliary nodes if one of nodes is outside
        if len(node_indices_dogbone)>0:
            if(pointA in node_indices_dogbone  and pointB not in node_indices_dogbone):
                pA = np.asarray( vor.points[pointA, :]  )
                pB = np.asarray( vor.points[pointB, :]  )
                ptB = (pA + pB)/2

                pointB = node_count + len(aux_nodes)
                aux_nodes.append(ptB)

            if(pointA not in node_indices_dogbone  and pointB in node_indices_dogbone):
                pA = np.asarray( vor.points[pointA, :]  )
                pB = np.asarray( vor.points[pointB, :]  )
                ptA = (pA + pB)/2

                pointA = node_count + len(aux_nodes)
                aux_nodes.append(ptA)



        if(pointA >= node_count and pointB<node_count):
            pA = np.asarray( vor.points[pointA, :]  )
            pB = np.asarray( vor.points[pointB, :]  )
            ptA = (pA + pB)/2

            pointA = node_count + len(aux_nodes)
            aux_nodes.append(ptA)

        if(pointB >= node_count  and pointA<node_count):
            pA = np.asarray( vor.points[pointA, :]  )
            pB = np.asarray( vor.points[pointB, :]  )
            ptB = (pA + pB)/2

            pointB = node_count + len(aux_nodes)
            aux_nodes.append(ptB)
        #
        rdg[0] = pointA
        rdg[1] = pointB

        #number of vertices
        rdg[2] = 2
        #indices of vertices
        rdg[3] = verticesIdxDict[vertA]
        rdg[4] = verticesIdxDict[vertB]
        #adding the ridge into the list of ridges
        ridges_out.append(rdg)
    #
    print(' - time:', time.time()-start)
    v_count = len (vertices_out)
    vertIdxStart = node_count + len(aux_nodes)

    ridges_out = np.array(ridges_out)
    ridges_out[:, 3] += vertIdxStart
    ridges_out[:, 4] += vertIdxStart

    print('done.')
    sys.stdout.flush()
    #output: nodes_out, aux_nodes, vertices_out, ridges_out


    saveNodes(master_folder, aux_nodes, "AuxNode",dim, auxNodesFile)
    if activeMechanics:
        saveNodes(master_folder, nodes_out, "Particle",dim, nodesFile)
        saveMechanicalElements(master_folder, ridges_out, node_count, dim, nodes_out, mZ=mZ, notches = notches, randomizeMaterial=randomizeMaterial)
    else:
        saveNodes(master_folder, nodes_out, "AuxNode",dim, nodesFile)
    if activeTransport:
        saveNodes(master_folder, vertices_out, "TrsprtNode",dim, verticesFile)
        saveTransportElements(master_folder, ridges_out,dim, node_count, v_count, aux_nodes, maxLim, nodes_out, vertices_out, coupled=coupled, mZ=mZ)
    else:
        saveNodes(master_folder, vertices_out, "AuxNode",dim, verticesFile)

    totalPointCount = len(nodes_out) + len(aux_nodes) + len(vertices_out)

    return v_count, verticesIdxDict, vertIdxStart, totalPointCount #, nodes_out, aux_nodes, vertices_out, ridges_out

def findClosest(points, target, dim):
    dist2 = np.zeros(len(points))
    for i in range(dim):
        dist2 += np.square(points[:,i]-target[i])
    index = np.argmin(dist2)
    return index,np.sqrt(dist2[index])

def savePeriodicBlock (master_folder,cpldNds, maxLim, nodes_out):
    cf = open(os.path.join(master_folder,blocksFile),"w")
#    print(cpldNds)
    nblocks = len(cpldNds)
    print("BLOCKS   ", nblocks)
    if (len(maxLim)==2):
        #loads=["\t2\tey\t0\tgxy\t1","\t2\tjy\t0\tjy\t0"]
        loads=["\t2\tey\t0\tgxy\t1","\t2\tgx\t0\tgy\t0"]
        names=["MechanicalPeriodicBC","TransportPeriodicBC"]
        for q in range(nblocks):
            ndepend = len(cpldNds[q])
            #ex ey gxy sx sy sxy
            cf.write("%s\tsize\t2\t%e\t%e\tload\t%s\tpairs\t%d"%(names[q],maxLim[0],maxLim[1],loads[q], ndepend))

            for i in range(len(cpldNds[q])):
                cf.write("\t%d\t%d"%(cpldNds[q][i][0], cpldNds[q][i][1]))

            #plt.plot( [nodes_out[ cpldNds[i][0],0 ], nodes_out[ cpldNds[i][1],0 ]], [nodes_out[ cpldNds[i][0],1 ], nodes_out[ cpldNds[i][1],1 ]],'ro-', color='red')
            #plt.text(nodes_out[ cpldNds[i][0],0 ] , nodes_out[ cpldNds[i][0],1 ], cpldNds[i][0], fontsize=11)
            #plt.text(nodes_out[ cpldNds[i][1],0 ] , nodes_out[ cpldNds[i][1],1 ], cpldNds[i][1], fontsize=11)

            cf.write(os.linesep)
        cf.close()

    if (len(maxLim)==3):
        #loads=["\t2\tey\t0\tgxy\t1","\t2\tjy\t0\tjy\t0"]
        loads=["\t2\tey\t0\tgxy\t1","\t3\tgx\t0\tgy\t0\tgz\t0"]
        names=["MechanicalPeriodicBC","TransportPeriodicBC"]
        for q in range(nblocks):
            ndepend = len(cpldNds[q])
            #ex ey gxy sx sy sxy
            cf.write("%s\tsize\t3\t%e\t%e\t%e\tload\t%s\tpairs\t%d"%(names[q],maxLim[0],maxLim[1],maxLim[2],loads[q], ndepend))

            for i in range(len(cpldNds[q])):
                cf.write("\t%d\t%d"%(cpldNds[q][i][0], cpldNds[q][i][1]))

            #plt.plot( [nodes_out[ cpldNds[i][0],0 ], nodes_out[ cpldNds[i][1],0 ]], [nodes_out[ cpldNds[i][0],1 ], nodes_out[ cpldNds[i][1],1 ]],'ro-', color='red')
            #plt.text(nodes_out[ cpldNds[i][0],0 ] , nodes_out[ cpldNds[i][0],1 ], cpldNds[i][0], fontsize=11)
            #plt.text(nodes_out[ cpldNds[i][1],0 ] , nodes_out[ cpldNds[i][1],1 ], cpldNds[i][1], fontsize=11)

            cf.write(os.linesep)
        cf.close()

    #plt.plot(nodes_out[:,0], nodes_out[:,1], 'o', color='black')

    # if SHOW_PLOT:
    #     plt.show()

def pointWithinCenterBox(point, maxLim):
    if (point[0]>0 and point[0]<maxLim[0]
    and point[1]>0 and point[1]<maxLim[1]
    and point[2]>0 and point[2]<maxLim[2]):
        return True
    else:
        return False

def ridgeWithinCenterBox(vertices, maxLim):
    within = True
    for i in range (len(vertices)):
        if pointWithinCenterBox(vertices[i], maxLim) == False:
            within = False

    return within

def output3D(master_folder, node_count, maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=None, notches=None, isTube=False, coupled=False, randomizeMaterial=False):
    start_time = time.time()
    dim = 3

    print('Extracting the geometry...',  end ='')
    sys.stdout.flush()

    printout = False
    # nody: [x,y,z] [powerR] [area]
    nodes_out = np.zeros( (node_count, (dim + 1 +1)))
    nodes_out[:,  0:dim] = node_coords[0:node_count,  0:dim]

    if ((len(areas) == node_count)):
       nodes_out[:,dim] = areas[:]

    #relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
    #print ('Area Error: %.5E ' %(relAreaError) )

    ########################################################################################################
    # ridges with nodes within sample
    validRidgeIdxs = []

    #print('ridge points')
    #adding ridges with at least one node in sample
    #validRidgeIdxs = np.where(np.any(vor.ridge_points < node_count, axis=1))[0].tolist()
    cond = np.any((vor.ridge_points < node_count) & (vor.ridge_points >= 0), axis=1)
    validRidgeIdxs = np.where(cond)[0]


    validRidgeIdxs = np.asarray(validRidgeIdxs)
    ########################################################################################################
    # vertices: [xA,yA,zA] [origIdx]
    vertices_out = []
    vertices_out_set = set()
    # dictionary of original and new indices of vertices
    verticesIdxDict = {}
    # ridges: nodeAidx, nodeBidx, trsprtBC, vertIdx
    ridges_out = []
    #auxiliary nodes
    aux_nodes = []
    ########################################################################################################
    allCoplanar = True
    for i in range (validRidgeIdxs.size):
        sys.stdout.write('\r'+'Ridge nr. ' + str(i) + '/' +  str(validRidgeIdxs.size)+'  '+          str(int(i/validRidgeIdxs.size*100))+'%')

        sys.stdout.flush()

        rdge = vor.ridge_vertices[validRidgeIdxs[i]]
        #indices of all vertices that form the planar ridge
        for j in range (len(rdge)):
            vrtx = np.zeros ( (dim + 1 +1 +1) )
            #
            #for d in range (dim):
            #    vrtx [d] = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][j]][d]
            vrtx[0:dim] =  vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][j]][0:dim]
            #
            vrtx[dim] = vor.ridge_vertices[validRidgeIdxs[i]][j]
            #
            if vrtx[dim] not in vertices_out_set:
                verticesIdxDict.update( { vrtx[dim] : len(vertices_out)  } )
                vrtx [dim +1] = len(vertices_out)
                vrtx [dim +2] = 0
                vertices_out.append(vrtx)
                vertices_out_set.add(vrtx[dim])

        #ridges
        ########################################################
        #array for the ridge: nodeA, nodeB, trBc, vertCount,newVertIdcs
        nrVertices = len(vor.ridge_vertices[validRidgeIdxs[i]])
        rdg = np.zeros ( (2 + 1 + nrVertices  ) )

        #nodes divided by the ridge
        pointA = vor.ridge_points[validRidgeIdxs[i]][0]
        pointB = vor.ridge_points[validRidgeIdxs[i]][1]

        #auxiliary nodes if one of them is out of sample
        if(pointA >= node_count and pointB<node_count):
            pA = np.asarray( vor.points[pointA, :]  )
            pB = np.asarray( vor.points[pointB, :]  )
            ptA = (pA + pB)/2

            pointA = node_count + len(aux_nodes)
            aux_nodes.append(ptA)


        if(pointB >= node_count  and pointA<node_count):
            pA = np.asarray( vor.points[pointA, :]  )
            pB = np.asarray( vor.points[pointB, :]  )
            ptB = (pA + pB)/2

            pointB = node_count + len(aux_nodes)
            aux_nodes.append(ptB)

        #
        rdg[0] = pointA
        rdg[1] = pointB
        rdg[2] = nrVertices

        #adding vert idcs
        for v in range ( nrVertices ):
            rdg[2+1+v] =  verticesIdxDict[ vor.ridge_vertices[validRidgeIdxs[i]][v] ]


        #coplanarity control
        maxE = 0
        for v in range ( nrVertices-3 ):
            pA = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v]][:]
            pB = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+1]][:]
            pC = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+2]][:]
            pD = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+3]][:]

            tol = 1e-10
            val = equation_plane(pA, pB, pC, pD)
            if (np.abs(val) > maxE): maxE = np.abs(val)
            if ( val > tol):
                allCoplanar = False
                print('Not coplanar!!! Ridge nr. %d, err: %e' %(i, val ))
            #else: print('Coplanar  %d' %i)

        #normal of the ridge surface from first three vertices
        planeNormal = getPlaneNormalVector(vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][0]][:],
                                     vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][1]][:],
                                     vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][2]][:])

        # vector connecting the nodes. Should be identical with plane normal. Otherwise the order of vertices will be swapped.
        pointNormal = vor.points[pointB] - vor.points[pointA]
        pointNormal /= np.linalg.norm(pointNormal)

        #print ('diff:')
        diff = np.linalg.norm(planeNormal - pointNormal)
        if (diff < 1e-10):
            if (printout): print ('Direction of plane normal OK')
        else:
            if (printout): print ('Direction of plane normal REVERSE')
            #rdg[nrVertices:] = rdg[:nrVertices-1:-1]
            #print('pred %s' %rdg)
            rdg[3:len(rdg)]= rdg[3:len(rdg)][::-1]
            #print('po %s \n' %rdg)

        ##############atan2((Vb x Va) . Vn, Va . Vb)##############
        #average point within the ridge surface
        """
        avgPoint = np.zeros(3)
        for d in range (3):
            for l in range ( nrVertices ):
                avgPoint [d] += vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][l]][d]
            avgPoint[d] /= len(vor.ridge_vertices[validRidgeIdxs[i]])
        #

        #mutual angles between vertices and the average point
        angles = np.zeros(nrVertices)
        referenceVector =  vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][0]][:] - avgPoint

        # computing the angles
        # atan2((Vb x Va) . Vn, Va . Vb)
        for l in range ( nrVertices ):
            currVector =  vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][l]][:] - avgPoint
            angles [l] = np.degrees( np.arctan2(  np.dot( np.cross( referenceVector, currVector ), planeNormal),
                                                np.dot(referenceVector,currVector)   ) )
            if (angles [l] < 0):
                angles [l] = 360 - (-angles [l])
        #print(angles)
        """

        ridges_out.append(rdg)
    geom_time = time.time()
    print ('done in %.1f seconds' %(geom_time-start_time))
    if (allCoplanar):
        print ('ALL ridges coplanar OK, maxErr: %e' %maxE)
    else:
        print ('!!! NOT ALL RIDGES COPLANAR, maxErr: %e' %maxE)
    vertIdxStart = node_count + len(aux_nodes)
    v_count = len (vertices_out)

    for i in range (len(ridges_out)):
        ln = len(np.asarray(ridges_out[i]) )
        for l in range (3, ln):
            ridges_out[i][l] += vertIdxStart



    newAuxNodes = 0
    if (activeTransport):
        newAuxNodes = saveTransportElements(master_folder, ridges_out,dim, node_count, v_count, aux_nodes, maxLim, nodes_out, vertices_out, isTube=isTube, coupled=coupled)
    vertIdxStart += newAuxNodes

    for i in range (len(ridges_out)):
        ln = len(np.asarray(ridges_out[i]) )
        for l in range (3, ln):
            ridges_out[i][l] += newAuxNodes


    if activeMechanics:
        saveNodes(master_folder, nodes_out, "Particle",dim, nodesFile)
        saveNodes(master_folder, aux_nodes, "AuxNode",dim, auxNodesFile)
        saveMechanicalElements(master_folder, ridges_out, node_count, dim, nodes_out, mZ=mZ, notches = notches, randomizeMaterial=randomizeMaterial)
    else:
        saveNodes(master_folder, nodes_out, "AuxNode",dim, nodesFile)

    if activeTransport:
        saveNodes(master_folder, vertices_out, "TrsprtNode",dim, verticesFile)
        saveNodes(master_folder, aux_nodes, "AuxNode",dim, auxNodesFile)
        #JM: save transport elements uz je volane drive
        # je potreba, aby bylo volane prvni, protoze jeste generuje nove aux nodes
        #saveTransportElements(master_folder, ridges_out,dim, node_count, aux_nodes, maxLim)
    else:
        saveNodes(master_folder, vertices_out, "AuxNode",dim, verticesFile)

    totalPointCount = len(nodes_out) + len(aux_nodes) + len(vertices_out)

    return v_count, verticesIdxDict, vertIdxStart, totalPointCount


def returnSelectedPtsRadial (innerRad , outerRad, points, axisDim=0, xmin = -1, xmax = -1):
    dim = 3
    #
    selectedPointIdxs = []
    #
    for i in range (len(points)):
        if (axisDim == 0):
            dist = np.sqrt(points[i,1]**2+points[i,2]**2)

        if (dist > innerRad and dist < outerRad):
            if xmin == -1 and xmax == -1:
                selectedPointIdxs.append(i)
            else:
                if xmin <= points[i,0] <= xmax :
                    selectedPointIdxs.append(i)


    return np.array(selectedPointIdxs).astype(int)


def returnSelectedPts_old (boundPtA , boundPtB, points):
    '''
    >>> boundPtA = np.array([0, 0])
    >>> boundPtB = np.array([1, 1])
    >>> points = np.array([[.5, .5], [.1, .9], [2, .1], [2, 3]])
    >>> returnSelectedPts_old(boundPtA, boundPtB, points)
    array([0, 1])
    '''
    dim = len (boundPtA)
    #
    selectedPointIdxs = []
    #
    for i in range (len(points)):
        selected = True
        #
        for d in range (dim):
            if (points[i][d] < boundPtA[d] or points[i][d] > boundPtB[d]):
                selected = False
        #
        if (selected == True):
            selectedPointIdxs.append(i)

    return np.array(selectedPointIdxs).astype(int)


def returnSelectedPts (boundPtA , boundPtB, points):
    '''
    >>> boundPtA = np.array([0, 0])
    >>> boundPtB = np.array([1, 1])
    >>> points = np.array([[.5, .5], [.1, .9], [2, .1], [2, 3]])
    >>> returnSelectedPts(boundPtA, boundPtB, points)
    array([0, 1])
    '''
    return np.where(np.all(~np.logical_or(points < boundPtA, points > boundPtB), axis=1))[0]


def excludeSelectedPts_old (boundPtA , boundPtB, points):
    '''
    >>> boundPtA = np.array([0, 0])
    >>> boundPtB = np.array([1, 1])
    >>> points = np.array([[.5, .5], [.1, .9], [2, .1], [2, 3]])
    >>> excludeSelectedPts_old(boundPtA, boundPtB, points)
    array([2, 3])
    '''
    dim = len (boundPtA)
    #
    selectedPointIdxs = []
    #
    for i in range (len(points)):
        selected = False
        #
        for d in range (dim):
            if (points[i][d] < boundPtA[d] or points[i][d] > boundPtB[d]):
                selected = True
        #
        if (selected == True):
            selectedPointIdxs.append(i)

    return np.array(selectedPointIdxs).astype(int)


def output3Dperiodic(master_folder, node_count, maxLim, vor, node_coords, areas, activeTransport, activeMechanics, minDist, mZ=None, notches=None, isTube=False, coupled=False):
    start_time = time.time()
    dim = len(maxLim)

    print('Extracting the geometry...',  end ='')
    sys.stdout.flush()

    # Mechanical Elements
    #######################################################################################################

    print ('Periodic model, filtering ridges...', end = '')
    is_inside = np.all(abs(vor.points-maxLim/2.)<=maxLim/2., axis=1)
    valid_node_idcs = np.where(is_inside)[0]
    inside_coords = vor.points[valid_node_idcs]
    is_positive = np.all( vor.points>=0, axis=1)
    #cross terms
    is_plus_only = np.zeros((len(vor.points), dim))
    for v in range(dim): is_plus_only[:,v] = ( vor.points[:,v]>maxLim[v])

    mechElemPoints = np.zeros((0,2)).astype(int)
    mechElemVerts = []

    coupledNodesMech = np.zeros((0,2)).astype(int)


    print('\nPocet nodu, se kterymi pocital voronoj: %d' %len(vor.points))

    maxIdx = -1
    for ir,r in enumerate(vor.ridge_points):
         if int(r[0]) > maxIdx:
             maxIdx = int(r[0])
         if int(r[1]) > maxIdx:
            maxIdx = int(r[1])
    print('Nejvyssi index nodu v ridges: %d' %maxIdx)

    if (len(vor.points)-1 == (maxIdx)):
        print ('Ridge spojuji jen nody v samplu.')
        print(' Tohle dela powerTes 0')
        print ('Export probehne v poradku.')
    elif (len(vor.points)-1 < (maxIdx)):
        print('Ridge se odkazuji na nejake nody s indexy, ktere nejsou v samplu. !!!')
        print(' Tohle dela powerTes 1')
        print('Nastane chyba index out of bounds...')
    print('Stiskni enter!\n')
    a = input('').split(" ")[0]


    coupledNodesMech = np.zeros((0,2)).astype(int)

    actual_node_count = len(valid_node_idcs)
    print ('actual node count: %d' %actual_node_count)


    for ir,r in enumerate(vor.ridge_points):
        sys.stdout.write('\rRidge nr.'+str(ir)+' / '+str(len(vor.ridge_points))+' ')
        sys.stdout.flush()

        nAidx = int(r[0])
        nBidx = int(r[1])

        is_diagonal = (is_positive[nAidx] and is_positive[nBidx]) and np.any(is_plus_only[nAidx]) and np.any(is_plus_only[nBidx])
        v = 0
        while is_diagonal and v<dim:
            is_diagonal =  ((not is_plus_only[nAidx,v] and not is_plus_only[nBidx,v]) or (is_plus_only[nAidx,v] != is_plus_only[nBidx,v]) )
            v += 1

        if  ( (is_inside[nBidx] and is_positive[nAidx]) or (is_inside[nAidx] and is_positive[nBidx]) or is_diagonal):
            for nXidx in [nAidx,nBidx]:
                foundnewslave = False
                if not nXidx in valid_node_idcs:
                    valid_node_idcs = np.hstack((valid_node_idcs, nXidx))
                    foundnewslave = True
                    newslave = nXidx
                if foundnewslave:
                    match = np.zeros(dim)
                    for i in range(dim):
                        if vor.points[newslave,i]>maxLim[i]:
                            match[i] = vor.points[newslave,i]-maxLim[i]
                        else: match[i] = vor.points[newslave,i]

                    dist = np.sum(np.square(inside_coords-match),axis=1)
                    master = np.argmin(dist)
                    if (dist[master]>1e-15):
                        print("Mechanical master not found, min square dist ", dist[master], vor.points[newslave] )
                        exit(1)
                    else:
                        coupledNodesMech = np.vstack( (coupledNodesMech, np.array([newslave, valid_node_idcs[master]]).astype(int) ))

            mechElemPoints = np.vstack((mechElemPoints, r))
            mechElemVerts.append(vor.ridge_vertices[int(ir)])

    # Transport Elements
    ########################################################################################################
    is_inside = np.all(abs(vor.vertices-maxLim/2.)<=maxLim/2., axis=1)
    valid_vert_idcs = np.where(is_inside)[0]
    inside_coords = vor.vertices[valid_vert_idcs]
    is_positive = np.all(np.column_stack((np.all( vor.vertices>=0, axis=1), np.all( vor.vertices<maxLim+5*minDist, axis=1))),axis=1)
    #cross terms
    is_plus_only = np.zeros((len(vor.vertices), dim))
    for v in range(dim): is_plus_only[:,v] = ( vor.vertices[:,v]>maxLim[v])

    trsprtElemVerts = np.zeros((0,2)).astype(int)
    trsprtElemNodes = []

    coupledNodesTrsp = np.zeros((0,2)).astype(int)

    for ir,r in enumerate(vor.ridge_vertices):
        vAidx = r[-1]
        for vBidx in r:
            is_diagonal = (is_positive[vAidx] and is_positive[vBidx]) and np.any(is_plus_only[vAidx]) and np.any(is_plus_only[vBidx])
            v = 0
            while is_diagonal and v<dim:
                is_diagonal =  ((not is_plus_only[vAidx,v] and not is_plus_only[vBidx,v]) or (is_plus_only[vAidx,v] != is_plus_only[vBidx,v]) )
                v += 1

            if ( (is_inside[vBidx] and is_positive[vAidx]) or (is_inside[vAidx] and is_positive[vBidx]) or is_diagonal):
                for vXidx in [vAidx,vBidx]:
                    foundnewslave = False
                    if not vXidx in valid_vert_idcs:
                        valid_vert_idcs = np.hstack((valid_vert_idcs, vXidx))
                        foundnewslave = True
                        newslave = vXidx
                    if foundnewslave:
                        match = np.zeros(dim)
                        for i in range(dim):
                            if vor.vertices[newslave,i]>maxLim[i]:
                                match[i] = vor.vertices[newslave,i]-maxLim[i]
                            else: match[i] = vor.vertices[newslave,i]

                        dist = np.sum(np.square(inside_coords-match),axis=1)
                        master = np.argmin(dist)
                        if (dist[master]>1e-6):
                            print("Trasnport master not found, min square dist ", dist[master], vor.vertices[newslave], match )
                            exit(1)
                        else:
                            coupledNodesTrsp = np.vstack( (coupledNodesTrsp, np.array([newslave, valid_vert_idcs[master]]).astype(int) ))


                trspVert = np.sort(np.array([vAidx,vBidx]).astype(int))
                fA = np.where(trsprtElemVerts[:,0]==trspVert[0])[0]
                fA = fA[ np.where(trsprtElemVerts[fA,1]==trspVert[1])[0] ]
                if (len(fA)==0):
                    trsprtElemVerts = np.vstack((trsprtElemVerts, trspVert))
                    trsprtElemNodes.append(vor.ridge_points[ir])
                else:
                   trsprtElemNodes[fA[0]] = np.hstack((trsprtElemNodes[fA[0]], vor.ridge_points[ir]))
            vAidx = vBidx


    #build surfaces of contact elements
    for i,p in enumerate(trsprtElemNodes):
        q = [p[0],p[1]]
        p = p[2:]

        while len(p)>0:
            idx = np.where(p==q[-1])[0]
            if (not len(idx)==1):
                cout << "Error in continuity of the transport facet" << endl;
                exit();
            if idx[0]%2==0:
                q.append(p[idx[0]+1])
                p = np.hstack((p[:idx[0]],p[idx[0]+2:]))
            else:
                q.append(p[idx[0]-1])
                p = np.hstack((p[:idx[0]-1],p[idx[0]+1:]))
        trsprtElemNodes[i] = q[:-1]

    ######################
    print('renumbering everything...')

    valid_node_idcs = np.array(valid_node_idcs).astype(int) #nodes with degrees of freedom

    numnodes = len(valid_node_idcs)
    mechauxnodes = np.unique(np.hstack(trsprtElemNodes))
    mapping = np.arange(len(mechauxnodes))
    trueauxindicators = np.isin(mechauxnodes,valid_node_idcs)
    trueaux  = np.where(trueauxindicators==False)[0]
    falseaux = np.where(trueauxindicators==True)[0]
    sort_idx = valid_node_idcs.argsort()
    mapping[falseaux] = sort_idx[np.searchsorted(valid_node_idcs,mechauxnodes[falseaux],sorter = sort_idx)]
    mapping[trueaux] = np.arange(len(trueaux))+numnodes


    sort_idx = mechauxnodes.argsort()
    for i in range(len(trsprtElemNodes)):
        m = sort_idx[np.searchsorted(mechauxnodes,trsprtElemNodes[i],sorter = sort_idx)]
        trsprtElemNodes[i] = mapping[m]
    mechauxnodes = mechauxnodes[trueaux]


    sort_idx = valid_node_idcs.argsort()
    mechElemPoints = sort_idx[np.searchsorted(valid_node_idcs,mechElemPoints,sorter = sort_idx)]

    coupledNodesMech = np.array(coupledNodesMech).astype(int)
    coupledNodesMech = sort_idx[np.searchsorted(valid_node_idcs,coupledNodesMech,sorter = sort_idx)]


    valid_vert_idcs = np.array(valid_vert_idcs).astype(int) #vertices with degrees of freedom
    numverts = len(valid_vert_idcs)
    trspauxnodes = np.unique(np.hstack(mechElemVerts))
    mapping = np.arange(len(trspauxnodes))
    trueauxindicators = np.isin(trspauxnodes,valid_vert_idcs)
    trueaux  = np.where(trueauxindicators==False)[0]
    falseaux = np.where(trueauxindicators==True)[0]
    sort_idx = valid_vert_idcs.argsort()
    mapping[trueaux] = np.arange(len(trueaux))+numnodes+len(mechauxnodes)
    mapping[falseaux] = sort_idx[np.searchsorted(valid_vert_idcs,trspauxnodes[falseaux],sorter = sort_idx)]+numnodes+len(mechauxnodes)+len(trueaux)

    sort_idx = trspauxnodes.argsort()
    for i in range(len(mechElemVerts)):
        m = sort_idx[np.searchsorted(trspauxnodes,mechElemVerts[i],sorter = sort_idx)]
        mechElemVerts[i] = mapping[m]
    trspauxnodes = trspauxnodes[trueaux]


    sort_idx = valid_vert_idcs.argsort()
    trsprtElemVerts = sort_idx[np.searchsorted(valid_vert_idcs,trsprtElemVerts,sorter = sort_idx)]+numnodes+len(mechauxnodes)+len(trueaux)


    coupledNodesTrsp = np.array(coupledNodesTrsp).astype(int)
    coupledNodesTrsp = sort_idx[np.searchsorted(valid_vert_idcs,coupledNodesTrsp,sorter = sort_idx)] + numnodes+len(mechauxnodes)+len(trueaux)

    """
    fullnodes = np.vstack((vor.points[valid_node_idcs], vor.points[mechauxnodes], vor.vertices[trspauxnodes], vor.vertices[valid_vert_idcs]))

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for p in trsprtElemVerts:
        ax.plot(fullnodes[p,0],fullnodes[p,1],fullnodes[p,2],color="k")

    for p in mechElemVerts:
        ax.plot(vor.vertices[p,0],vor.vertices[p,1],vor.vertices[p,2],color="b", ls=":")

    plt.show()
    """

    savePeriodicBlock(master_folder,[coupledNodesMech,coupledNodesTrsp],maxLim, vor.points[mechauxnodes])

    if activeMechanics:
        particles = np.column_stack((vor.points[valid_node_idcs], np.zeros(len(valid_node_idcs))))
        saveNodes(master_folder, particles, "Particle",dim, nodesFile)
    else:
        saveNodes(master_folder, vor.points[valid_node_idcs], "AuxNode",dim, nodesFile)
    if activeTransport:
        saveNodes(master_folder, vor.vertices[valid_vert_idcs], "TrsprtNode",dim, verticesFile)
    else:
        saveNodes(master_folder, vor.vertices[valid_vert_idcs], "AuxNode",dim, verticesFile)
    saveNodes(master_folder, np.vstack(( vor.points[mechauxnodes], vor.vertices[trspauxnodes])), "AuxNode",dim, auxNodesFile)


    if activeMechanics:
        inpf = open(os.path.join(master_folder,mechElemsFile),"w")
        inpf.write("#ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tverticesIdxs\tMaterial\n")
        for k in range(len(mechElemPoints)):
            inpf.write("LTCBEAM\t%d\t%d\t%d"%(mechElemPoints[k,0],mechElemPoints[k,1], len(mechElemVerts[k]) ))
            for p in mechElemVerts[k]:
                inpf.write("\t%d"%(p))
            inpf.write("\t0\n")
        inpf.close()

    if activeTransport:
        inpf = open(os.path.join(master_folder,trsprtElemsFile),"w")
        inpf.write("#ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tverticesIdxs\tMaterial\n")
        for k in range(len(trsprtElemVerts)):
            inpf.write("LTCTRSP\t%d\t%d\t%d"%(trsprtElemVerts[k,0],trsprtElemVerts[k,1], len(trsprtElemNodes[k]) ))
            for p in trsprtElemNodes[k]:
                inpf.write("\t%d"%(p))
            inpf.write("\t1\n")
        inpf.close()


    totalPointCount = len(valid_node_idcs) + len(mechauxnodes) + len(trspauxnodes) + len(valid_vert_idcs)
    v_count = len(valid_vert_idcs)
    vertIdxStart = totalPointCount-v_count

    checkSavedModel(master_folder, dim, activeMechanics, activeTransport)

    return v_count, [], vertIdxStart, totalPointCount


def output2DPeriodic(master_folder, node_count,  maxLim, vor, node_coords, areas, activeMechanics, activeTransport, minDist, mZ=None):
    return output3Dperiodic(master_folder, node_count,  maxLim, vor, node_coords, areas, activeMechanics, activeTransport, minDist, mZ=mZ)









































def excludeSelectedPts (boundPtA , boundPtB, points):
    '''
    >>> boundPtA = np.array([0, 0])
    >>> boundPtB = np.array([1, 1])
    >>> points = np.array([[.5, .5], [.1, .9], [2, .1], [2, 3]])
    >>> excludeSelectedPts(boundPtA, boundPtB, points)
    array([2, 3])
    '''
    return np.where(np.any(np.logical_or(points < boundPtA, points > boundPtB), axis=1))[0]


def saveTransportIC(master_folder,transportIC_merged):
    print('Saving TRSPRT initial conditions...', end ='')
    sys.stdout.flush()
    trsprtIC_out = []

    for i in range (len(transportIC_merged)):
        ic = np.zeros (2)
        ic[0] = transportIC_merged[i].getVrtxIdx()
        ic[1] = transportIC_merged[i].getPressure()
        trsprtIC_out.append(ic)

    headerLine = 'vrtxIdx\tpressure'
    fl=open(os.path.join(master_folder,trsprtICFile) ,'w')
    np.savetxt(fl, trsprtIC_out, delimiter='\t', fmt='%d\t%f', header = headerLine)
    fl.close()
    print('done.')



def saveMechIC(master_folder,dim, nodes_mechICmerged):
    print('Saving MECH initial conditions...', end='')
    sys.stdout.flush()
    mechIC_out = []

    for i in range (len(nodes_mechICmerged)):
        if (dim == 2):
            ic = np.zeros (1 + 3)
            ic[0] = nodes_mechICmerged[i].getNodeIdx()
            ic[1:] = nodes_mechICmerged[i].getMechIC()
        elif (dim == 3):
            ic = np.zeros (1 + 6)
            ic[0] = nodes_mechICmerged[i].getNodeIdx()
            ic[1:] = nodes_mechICmerged[i].getMechIC()
        mechIC_out.append(ic)

    #
    if (dim == 2):
        headerLine = 'nodeIdx\tTrVelocX\tTrVelocY\tRotVelocZ'
        fl=open(os.path.join(master_folder,mechICFile) ,'w')
        np.savetxt(fl, mechIC_out, delimiter='\t', fmt='%d\t%f\t%f\t%f', header = headerLine)
        fl.close()
    elif (dim == 3):
        headerLine = 'nodeIdx\tTrVelocX\tTrVelocY\tTrVelocZ\tRotVelocX\tRotVelocY\tRotVelocZ'
        fl=open(os.path.join(master_folder,mechICFile) ,'w')
        np.savetxt(fl, mechIC_out, delimiter='\t', fmt='%d\t%f\t%f\t%f\t%f\t%f\t%f', header = headerLine)
        fl.close()

    print('done.')

def saveRadii(master_folder,radii):
    headerLine = 'radii'
    fl=open(os.path.join(master_folder,"radii.out") ,'w')
    np.savetxt(fl, radii, fmt='%e', header = headerLine)
    fl.close()

def saveMechBC(master_folder,dim, nodes_mechBCmerged, govNodesBC = False):
    print('Saving MECH boundary conditions...', end='')
    sys.stdout.flush()

    #print (len(nodes_mechBCmerged))
    mechBC_out = []

    for i in range (len(nodes_mechBCmerged)):
        if (dim == 2):
            bc = np.zeros (1 + 6)
            bc[0] = nodes_mechBCmerged[i].getNodeIdx()
            bc[1:7] = nodes_mechBCmerged[i].getMechBC()[0:6]
            if (len(nodes_mechBCmerged[i].getMechBC()) ==  7):
                bc[0] += nodes_mechBCmerged[i].getMechBC()[6]

        elif (dim == 3):
            bc = np.zeros (1 + 12)
            bc[0] = nodes_mechBCmerged[i].getNodeIdx()
            bc[1:] = nodes_mechBCmerged[i].getMechBC()

        mechBC_out.append(bc)

    #
    if (govNodesBC == False and len(mechBC_out)>0 ):
        print('Saving node BC...')
        if (dim == 2):
            headerLine = 'bcType\tnodeIdx\tKinTrX\tKinTrY\tKinRotZ\tStTrX\tStTrY\tStRotZ'
            fl=open(os.path.join(master_folder,mechBCFile) ,'w')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='NodalBC\t%d\t%d\t%d\t%d\t%d\t%d\t%d', header = headerLine)
            fl.close()
        elif (dim == 3):
            headerLine = 'bcType\tnodeIdx\tKinTrX\tKinTrY\tKinTrZ\tKinRotX\tKinRotY\tKinRotZ\tStTrX\tStTrY\tStTrZ\tStRotX\tStRotY\tStRotZ'
            fl=open(os.path.join(master_folder,mechBCFile) ,'w')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='NodalBC\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d', header = headerLine)
            fl.close()
    elif(len(mechBC_out)>0 ):
        print('Saving rigid plate BC...')
        #print(mechBC_out)
        if (dim == 2):
            fl=open(os.path.join(master_folder,mechBCFile) ,'a')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='NodalBC\t%d\t%d\t%d\t%d\t%d\t%d\t%d')
            fl.close()
        elif (dim == 3):
            fl=open(os.path.join(master_folder,mechBCFile) ,'a')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='NodalBC\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d')
            fl.close()

    print('done.')



def saveTransportBC(master_folder,transportBCmerged, verticesDict, vertIdxStart, govNodesBC = False, totalNodeCount=-1):
    print('Saving TRSPRT boundary conditions...', end = '')

    sys.stdout.flush()
    trsptBC_out = []

    ttlNdcnt = totalNodeCount

    for i in range (len(transportBCmerged)):
        idx = transportBCmerged[i].getNodeIdx()

        if (govNodesBC==True):
            bc = np.zeros ((1 + 1 + 1))
            bc[0] = transportBCmerged[i].getNodeIdx()
            bc[1:] = transportBCmerged[i].getTrsprtBC()
            trsptBC_out.append(bc)

        if (govNodesBC==False and idx in verticesDict):
            bc = np.zeros ((1 + 1 + 1))
            bc[0] = verticesDict[transportBCmerged[i].getNodeIdx()] + vertIdxStart
            bc[1:] = transportBCmerged[i].getTrsprtBC()
            trsptBC_out.append(bc)



    print('len len %d' %len(trsptBC_out))

    if (govNodesBC==False and len(trsptBC_out)>0):
        headerLine = 'bcType\tvrtxIdx\tTrsptP\tTrsptJ'
        fl=open(os.path.join(master_folder,trsprtBCFile) ,'w')
        if (len(trsptBC_out)>0):
            np.savetxt(fl, trsptBC_out, delimiter='\t', fmt='NodalBC\t%d\t%d\t%d', header = headerLine)
        fl.close()
    elif(len(trsptBC_out)>0):
        print('Saving trspt rigid plate bc...')
        fl=open(os.path.join(master_folder,trsprtBCFile) ,'w')
        if (len(trsptBC_out)>0):
            np.savetxt(fl, trsptBC_out, delimiter='\t', fmt='NodalBC\t%d\t%d\t%d')
        fl.close()

    print('done.')





def saveSolver(master_folder, solver, solStep, minStep, maxStep, simTime, limitTolerance= 1e-1, maxIt=20, tolerance = 1e-3):
    f=open(os.path.join(master_folder,solverFile),'w')
    f.write('%s\n'%solver)
    f.write('time_step\t%e\n'%solStep)
    f.write('min_time_step\t%e\n'%minStep)
    f.write('max_time_step\t%e\n'%maxStep)
    f.write('total_time\t%f\n'%simTime)
    f.write('limit_tolerance\t%e\n'%limitTolerance)
    f.write('maxIt\t%d\n'%maxIt)
    f.write('tolerance\t%e\n'%tolerance)

    f.close()

def saveMasterInput(master_folder,dim, solver, solStep, minStep, maxStep, simTime, activeTransport, activeMechanics, periodic=False, constraint=False, constraintTrspt=False, limitTolerance= 1e-1, maxIt=20, tolerance = 1e-3):
     print('Saving master file...', end='')
     sys.stdout.flush()
     fl=open(os.path.join(master_folder,masterFile),'w')

     fl.write("Dimension\t%d\n"%dim)
     fl.write("Solver\t%s\n"%(solverFile))
     saveSolver(master_folder, solver, solStep, minStep, maxStep, simTime, limitTolerance, maxIt, tolerance=tolerance)

     if not periodic:

         """
         if (solver == "SteadyStateLinearSolver"):
                fl.write('Solver\tSteadyStateLinearSolver\ttime_step\t%e\ttotal_time\t%e\tconj_grad_precission\t1e-18
conj_grad_relative_maxit\t2\n' %(solStep, simTime))
         if (solver == "SteadyStateNonLinearSolver"):
                fl.write('Solver\tSteadyStateNonLinearSolver\ttime_step\t%e\tmax_time_step\t%e\tmin_time_step\t%e\ttotal_time\t%e\tlimit_tolerance\t%e\tmaxIt\t%d\tconj_grad_precission\t1e-18
conj_grad_relative_maxit\t2\n' %(solStep,  maxStep, minStep, simTime, limitTolerance, maxIt))
         """

         if not constraint:
             fl.write("NodeFiles\t3\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile))
         else:
             if (constraint and not constraintTrspt):
                 fl.write("NodeFiles\t4\t%s\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile, govNodesFile))
                 fl.write('PBlockFiles\t1\t%s\n' %(constraintFile))
             elif (not constraint and constraintTrspt):
                 fl.write("NodeFiles\t4\t%s\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile, govNodesTrsptFile))
                 fl.write('PBlockFiles\t1\t%s\n' %(constraintTrsptFile))
             elif (constraint and constraintTrspt):
                 fl.write("NodeFiles\t5\t%s\t%s\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile, govNodesFile, govNodesTrsptFile))
                 fl.write('PBlockFiles\t2\t%s\t%s\n' %(constraintFile, constraintTrsptFile))

         fl.write("MatFiles\t1\t%s\n"%materialsFile)
         if (activeTransport and activeMechanics):
             fl.write("ElemFiles\t2\t%s\t%s\n"%(mechElemsFile,trsprtElemsFile))
             fl.write("BCFiles\t2\t%s\t%s\n"%(mechBCFile,trsprtBCFile))
         elif  (activeTransport):
             fl.write("ElemFiles\t1\t%s\n"%(trsprtElemsFile))
             fl.write("BCFiles\t1\t%s\n"%(trsprtBCFile))
         elif  (activeMechanics):
             fl.write("ElemFiles\t1\t%s\n"%(mechElemsFile))
             fl.write("BCFiles\t1\t%s\n"%(mechBCFile))

     else:
         if (os.path.isfile(os.path.join(master_folder,auxNodesFile))):
             fl.write("NodeFiles\t3\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile))
         else:
             fl.write("NodeFiles\t2\t%s\t%s\n"%(nodesFile,verticesFile))
         fl.write("MatFiles\t1\t%s\n"%materialsFile)
         if (activeTransport and activeMechanics):
             fl.write("ElemFiles\t2\t%s\t%s\n"%(mechElemsFile,trsprtElemsFile))
         elif  (activeTransport):
             fl.write("ElemFiles\t1\t%s\n"%(trsprtElemsFile))
         elif  (activeMechanics):
             fl.write("ElemFiles\t1\t%s\n"%(mechElemsFile))
         fl.write("PBlockFiles\t1\t%s\n" %blocksFile)
     fl.write("FunctionFiles\t1\t%s\n"%functionsFile)
     fl.write("ExporterFiles\t1\t%s"%exportersFile)
     #fl.write('INITMECH:\t%s\n' % initConditionsMechFile   )
     #fl.write('INITTRSPRT:\t%s\n' % initConditionsTrsprtFile   )

     fl.close()
     print('done.')

def saveMaterials (master_folder,materials, regime = 'w'):
    print()
    print ('Saving materials...', end='')
    sys.stdout.flush()

    if not os.path.exists(os.path.join(master_folder,materialsFile)):
        with open(os.path.join(master_folder,materialsFile), regime) as f:
            for item in materials:
                f.write("%s\n" % item.getString() )
    else:
        print ('materials already saved, skipping...', end='')

    print('done')

def saveFunctions (master_folder,functions):
    print ('Saving functions...', end='')
    sys.stdout.flush()
    ### FUNCTIONS
    with open(os.path.join(master_folder,functionsFile), 'w') as f:
        headerLine = '#FuncTyped'
        f.write("%s\n" % headerLine )
        for item in functions:
            f.write("%s\n" % item.getString() )
          # print (item.getString())
    print('done.')





def saveExporters(master_folder,activeTransport, activeMechanics, exporters=[]):
    if len(exporters)==0:
        print('Saving default exporters...', end='')
        sys.stdout.flush()
        fl=open(os.path.join(master_folder,exportersFile),'w')
        if activeMechanics:
            fl.write('#TXTNodalExporter translations 2 ux uy\n')
            fl.write('#TXTNodalExporter pressure 1 pressure\n')
            if not activeTransport:
                fl.write('VTKElementExporter out  saveEvery 1e-4 cellData 2 damage crack_opening pointData 1 nodal_stress\n')
            fl.write('#VTKRCExporter faces  saveEvery 1e-4 cellData 1 damage\n')
            fl.write('#TXTGaussPointExporter damageT 11 x y z normal_x normal_y normal_z damage strainTY strainTZ strainPLTY strainPLTZ\n')
        if activeTransport:
            fl.write('TXTNodalExporter pressure 1 pressure\n')
            fl.write('VTKElementExporter elems saveEvery 0.0001 cellData 2 damage crack_opening pointData 1 pressure\n')

        fl.close()

    else:
        print('Saving exporters specified in prep_master...')
        sys.stdout.flush()
        fl=open(os.path.join(master_folder,exportersFile),'w')
        for exporter in exporters:
            print(exporter)
            for e in exporter:
                fl.write('%s\t'%e)

            fl.write('\n')
        fl.close()

    print('done.')



def saveNodes (master_folder,nodes_out, nodetype, dim, filename, virtualDoF=0):
    print('Saving nodes: %s...' %nodetype, end='')
    sys.stdout.flush()

    nodes_out = np.array(nodes_out)
    #writing nodes
    #print(len(nodes_out))
    num = dim

    if (dim == 2):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY"
        fmt= nodetype + '\t%.15e\t%.15e'
        if nodetype=="GovParticle":
            fmt= 'Particle' + '\t%.15e\t%.15e'

    elif (dim == 3):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY\tnodeCrdZ"
        fmt= nodetype + '\t%.15e\t%.15e\t%.15e'
        if nodetype=="GovParticle":
            fmt= 'Particle' + '\t%.15e\t%.15e\t%.15e'
    if nodetype=="Particle":
        headerLine = headerLine + "\tpowRadius"
        fmt = fmt + '\t%.15e'
        num = num + 1


    fl=open(os.path.join(master_folder,filename),'w')
    if len(nodes_out) > 0:
        np.savetxt(fl,  nodes_out[:,:num], delimiter='\t',   fmt=fmt,  header = headerLine)

    if virtualDoF !=0:
        for i in range (virtualDoF):
            fl.write('MechDoF\n')

    fl.close()
    print('done.')
    sys.stdout.flush()



def saveMechanicalElements (master_folder,ridges_out, node_count, dim, nodes, mZ=None, notches = None, randomizeMaterial = False):
    print('Saving MECH elements...', end ='')
    sys.stdout.flush()
    #filtering ridges to ridges with both nodes in sample -> mech elements
    mechElemRidges = []
    for m in range (len(ridges_out)):
        if (ridges_out[m][0] < node_count and ridges_out[m][1] < node_count and ridges_out[m][0] >=0  and ridges_out[m][1] >= 0):
            mechElemRidges.append( ridges_out[m].copy() )
    print ('Mech elements: %d' %len(mechElemRidges))

    onlyMechNodesConnected = True
    elaElems = []
    fig, ax = plt.subplots()


    if (mZ!=None and len(mZ)>0):
        print('Material zones recognized.')


        for i in range (len(mechElemRidges)):
            nodeA = nodes[int(mechElemRidges[i][0])]
            nodeB = nodes[int(mechElemRidges[i][1])]

            if (int(mechElemRidges[i][0]) >= node_count or int(mechElemRidges[i][1]) >= node_count):
                onlyMechNodesConnected = False

            if (dim==2):

                #rebars
                if (mZ[0][0]=='circle'):
                    inRebar = False
                    for rebar in range (mZ[0][3]):
                        if (np.linalg.norm(nodeA[0:2]-mZ[rebar][2]) < mZ[rebar][1] ) and (np.linalg.norm(nodeB[0:2]-mZ[rebar][2]) < mZ[rebar][1] ):
                            inRebar = True
                    if inRebar:
                        mechElemRidges[i] = np.hstack( (mechElemRidges[i], np.array([2])) )
                    else:
                        mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )

                elif ( (mZ[0][0][0] < nodeA[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeA[1] < mZ[0][1][1] and
                      mZ[0][0][0] < nodeB[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeB[1] < mZ[0][1][1]) or
                      (mZ[0][2][0] < nodeA[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeA[1] < mZ[0][3][1] and
                      mZ[0][2][0] < nodeB[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeB[1] < mZ[0][3][1])   ):
                    mechElemRidges[i] = np.hstack( (mechElemRidges[i], np.array([2])) )

                else:
                    mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )



            if (dim==3):
                triangle = False
                if (mZ[0][0][0] > mZ[0][1][0]):
                    triangle = True

                #print (mZ[0][0][0])
                #print (mZ[0][1][0])
                #print (mZ[0][0][0] - mZ[0][1][0])
                #print (triangle)

                if ( triangle == False and
                      ((mZ[0][0][0] < nodeA[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeA[1] < mZ[0][1][1] and
                      mZ[0][0][2] < nodeA[2] < mZ[0][1][2] and
                      mZ[0][0][0] < nodeB[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeB[1] < mZ[0][1][1] and
                      mZ[0][0][2] < nodeB[2] < mZ[0][1][2] )
                      or
                      mZ[0][2][0] < nodeA[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeA[1] < mZ[0][3][1] and
                      mZ[0][2][2] < nodeA[2] < mZ[0][3][2] and
                      mZ[0][2][0] < nodeB[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeB[1] < mZ[0][3][1] and
                      mZ[0][2][2] < nodeB[2] < mZ[0][3][2]) ):
                      mechElemRidges[i] =  np.hstack( (mechElemRidges[i], np.array([2])) )
                #print('len mz %d' %len(mZ))

                if (  len(mZ)==2 and
                      ((mZ[1][0][0] < nodeA[0] < mZ[1][1][0] and
                      mZ[1][0][1] < nodeA[1] < mZ[1][1][1] and
                      mZ[1][0][2] < nodeA[2] < mZ[1][1][2] and
                      mZ[1][0][0] < nodeB[0] < mZ[1][1][0] and
                      mZ[1][0][1] < nodeB[1] < mZ[1][1][1] and
                      mZ[1][0][2] < nodeB[2] < mZ[1][1][2] )) ):
                      mechElemRidges[i] =  np.hstack( (mechElemRidges[i], np.array([3])) )
                      print('node in notch')

                if ( triangle == True ):
                      isPresentA = False
                      isPresentB = False

                      xmin = mZ[0][1][0]
                      xmax = mZ[0][0][0]
                      ymin = mZ[0][0][1]
                      ymax = mZ[0][1][1]
                      zmin = mZ[0][0][2]
                      zmax = mZ[0][1][2]

                      xTop = (xmin+xmax)/2 - xmin
                      yTop = ymax - ymin


                      x = nodeA[0] - xmin
                      y = ymax - nodeA[1]
                      z = nodeA[2]

                      #print ('cx %f' %(xmin+xTop*0.5/yTop))
                      #print ('cy %f' %(xmax-yTop*0.05/xTop))
                      if (nodeA[0] > xmin and nodeA[0]<(xmin+xmax)/2):
                          xMINlim = xmin + xTop * y / yTop
                          yMINlim = ymax - yTop * x / xTop

                          if (  nodeA[0]>xMINlim and nodeA[1]>yMINlim and zmin<z<zmax):
                              isPresentA = True
                              #ax.scatter(nodeA[0], nodeA[1])
                              #print ('x minlim %f' %xMINlim)
                              #print ('y minlim %f' %yMINlim)

                      if (nodeA[0] > (xmin+xmax)/2 and nodeA[0] < xmax ):
                          xMAXlim = xmax - xTop * y / yTop
                          yMINlim = ymax - yTop * (xmax-nodeA[0]) / xTop

                          if (  nodeA[0]<xMAXlim and nodeA[1]>yMINlim and zmin<z<zmax):
                              isPresentA = True
                              #ax.scatter(nodeA[0], nodeA[1])
                              #print ('xtop %f' %xTop)
                              #print ('ytop %f' %yTop)
                              #print ('x %f < %f' %(nodeA[0], xMAXlim))
                              #print ('y %f > %f' %(nodeA[1], yMINlim))

                              #a = input('').split(" ")[0]


                      x = nodeB[0] - xmin
                      y = ymax - nodeB[1]
                      z = nodeB[2]

                      if (nodeB[0] > xmin and nodeB[0]<xTop):
                          xMINlim = xmin + xTop * y / yTop
                          yMINlim = ymax - yTop * x / xTop

                          if (  nodeB[0]>xMINlim and nodeB[1]>yMINlim and zmin<z<zmax):
                              isPresentB = True

                      if (nodeB[0] > xTop and nodeB[0] < xmax ):
                          xMAXlim = xmax - xTop * y / yTop
                          yMINlim = ymax - yTop * (xmax-nodeA[0]) / xTop

                          if (  nodeB[0]<xMAXlim and nodeB[1]>yMINlim and zmin<z<zmax):
                              isPresentB = True

                      if isPresentA and isPresentB:
                          mechElemRidges[i] =  np.hstack( (mechElemRidges[i], np.array([2])) )
                          elaElems.append(nodeA[0:3])
                          elaElems.append(nodeB[0:3])

                          #ax.scatter(nodeB[0], nodeB[1])
                          #print('node in triangle')


                else:
                      mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )
                """
                elif (mZ[0][2][0] < nodeA[0] < mZ[0][3][0] and
                          mZ[0][2][1] < nodeA[1] < mZ[0][3][1] and
                          mZ[0][2][2] < nodeA[2] < mZ[0][3][2] and
                          mZ[0][2][0] < nodeB[0] < mZ[0][3][0] and
                          mZ[0][2][1] < nodeB[1] < mZ[0][3][1] and
                          mZ[0][2][2] < nodeB[2] < mZ[0][3][2]) :
                          mechElemRidges[i] =  np.hstack( (mechElemRidges[i], np.array([2])) )
                          print('node in B')
                        #print('found ela element')
                """

    else:
        for i in range (len(mechElemRidges)):
            mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )

    if (onlyMechNodesConnected):
        print ('MechElems connect only MechNodes. That is ok.')
    else:
        print ('MechElems CONNECT WRONG NODES !!!')

    np.set_printoptions(threshold=np.inf)
    #print()
    np.asarray(elaElems)


    #plt.show()


    if (notches!=None ):
        print('Filtering out elements connecting notches...' )
        elementsWithoutNotch = []
        for i in range (len(mechElemRidges)):
            nA = mechElemRidges[i][0]
            nB = mechElemRidges[i][1]

            addElem = True
            for notch in notches:
                if ((nA in notch[0] and nB in notch[1]) or ((nB in notch[0] and nA in notch[1]))):
                    addElem = False
                    break

            if addElem == True:
                elementsWithoutNotch.append(mechElemRidges[i])

        mechElemRidges = elementsWithoutNotch

        print('done.')


    if randomizeMaterial == True:
        randomizeMechMaterial(master_folder, mechElemRidges, nodes, materialType='mars')

    if (dim == 2):
        headerLine = 'ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tvrtxAIdx\tvrtxBIdx\tMaterial'
        fl=open(os.path.join(master_folder,mechElemsFile),'w')
        #print(mechElemRidges[0])
        np.savetxt(fl, mechElemRidges, delimiter='\t',fmt='LTCBEAM\t%d\t%d\t%d\t%d\t%d\t%d', header = headerLine )
        fl.close()

    if (dim == 3):
        headerLine = '#ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tverticesIdxs\tMaterial\n'
        fl=open(os.path.join(master_folder,mechElemsFile),'w')
        ro = np.asarray(mechElemRidges[0])
        fl.write(headerLine)
        for i in range (len(mechElemRidges)):
            ro = np.array(mechElemRidges[i], ndmin=2)
            fmt='LTCBEAM\t%d\t%d\t%d'
            np.savetxt(fl,  ro, delimiter='\t', fmt=fmt+'\t%d'*(ro.shape[1]-3)+ '\t0')

    sys.stdout.flush()


def randomizeMechMaterial (master_folder, mechanicalRidges, nodes, materialType='mars'):
    print('Randomizing materials...')
    materials = []
    for i in range(len(mechanicalRidges)):
        sys.stdout.write('\r'+'Mech elem. nr. ' + str(i) + '/' + str(len(mechanicalRidges)))
        sys.stdout.flush()

        element = mechanicalRidges[i]
    #    print(element)

        nodeA = nodes[int(element[0])]
        nodeB = nodes[int(element[1])]
        integrationPoint = (nodeA-nodeB) /2

        material = getRandomizedMaterialProperties(integrationPoint, materialType)
        materials.append(material)

        #zmena indexu materialu prvku
        element[-1] = i



    saveMaterials(master_folder, materials)



def getRandomizedMaterialProperties(integrationPoint, materialType='mars'):
    #initial material properties
    if materialType == 'mars':
        initYoung = 30e9
        initAlpha = 0.25
        initDensity = 2200
        initFt = 2.5e6
        initGt = 50

        #tady se zrandomizuji parametry materialu
        randCoef = getRandCoef(integrationPoint)
        myYoung = initYoung * randCoef
        myAlpha = initAlpha * randCoef
        myDensity = initDensity * randCoef
        myFt = initFt * randCoef
        myGt = initGt * randCoef

        material =  utilitiesMech.MarsMaterial(myYoung, myAlpha, myDensity, myFt, myGt)

    return material

def getRandCoef(integrationPoint):

    randCoef = np.random.uniform(low=0.9, high=1.1)

    return randCoef






def saveTransportElements(master_folder,ridges_out, dim, node_count, vertCount, aux_nodes, maxLim, nodes_out, vertices_out, isTube=False, coupled=False, mZ=[]):
    print('Creating TRSPRT elements...', end='')
    sys.stdout.flush()
    transportElements = []
    transportElements_dict = {}
    ridges_out = np.asarray(ridges_out)

    onlyVerticesConnected = True

    nds = node_count
    aux = len(aux_nodes)
    vrt = vertCount


    if (dim == 2):
        for i in range (len(ridges_out)):
            connNds = []
            connNds.append (ridges_out[i,0])
            connNds.append (ridges_out[i,1])
            vrtA = ridges_out[i,3]
            vrtB = ridges_out[i,4]
            trp = utilitiesMech.transportPath (vrtA, vrtB, connNds, 1)
            transportElements.append (trp)

            if (vrtA <node_count or vrtA>=(nds+aux+vrt)) or  (vrtB <node_count or vrtB>=(nds+aux+vrt)):
                onlyVerticesConnected = False

    if (dim==3):
        for i in range (len(ridges_out)):
            ro = np.asarray(ridges_out[i])

            vrtA = ro[3]
            vrtB = ro[4]
            if (vrtA <node_count or vrtA>=(nds+aux+vrt)) or  (vrtB <node_count or vrtB>=(nds+aux+vrt)):
                onlyVerticesConnected = False

            #print (ro)
            for n in range (3, len(ro)):
                newPath = True
                m = n+1
                if (n==len(ro)-1):
                    m = 3
                path_ends = frozenset((ro[m], ro[n]))
                #print('%d ; %d = %d ; %d' %(n,m, ro[n], ro[m]))
                elem = transportElements_dict.get((path_ends), None)
                if elem:
                    elem.addConnectedNodes(ro)
                else:
                    connNds = []
                    connNds.clear()
                    connNds.append (ro[0])
                    connNds.append (ro[1])
                    #if (ro[0]>node_count and ro[1]>node_count):
                    #    print ('both aux')
                    #transportElements.append (utilitiesMech.transportPath (ro[n], ro[m], connNds.copy(), 1))
                    transportElements_dict[path_ends] = utilitiesMech.transportPath (ro[n], ro[m], connNds.copy(), 1)
            transportElements = transportElements_dict.values()

    if (onlyVerticesConnected):
        print('Transport elements connect only vertices. That is ok.')
    else:
        print('Transport elements CONNECT WRONG POINTS !!!!')
    sys.stdout.flush()

    #  i[b], i[a] = i[a], i[b]
    print ('Reordering connected nodes...', end='')
    allReorderedFine = True
    for elem in transportElements:
        #print('\n\n\n')
        restart = True
        while restart:
            restart = False
            for i in range (1, len( elem.connectedNodes )-2, 2) :
                #print(i)
                if (elem.connectedNodes[i] != elem.connectedNodes[i+1]):
                    #print(elem.connectedNodes, end='')
                    #print(' size %d' % len(elem.connectedNodes))
                    simIdx = -1
                    for j in range (i+1, len(elem.connectedNodes)):
                        #print('search for %d ' %elem.connectedNodes[i] , end='')
                        #print(elem.connectedNodes[j])
                        if (elem.connectedNodes[i] == elem.connectedNodes[j]):
                            simIdx = j
                            #print('simIdx: %d ' %simIdx)
                            break
                    if (simIdx == -1):
                        #print ("___ Not found ", end='')
                        #swapnout prvni dva nody
                        #print ('First couple swapping %d for %d' %(0,1))
                        elem.connectedNodes[0], elem.connectedNodes[1] =  elem.connectedNodes[1], elem.connectedNodes[0]
                        elem.connectedNodes[2], elem.connectedNodes[3] =  elem.connectedNodes[3], elem.connectedNodes[2]
                        #
                        elem.connectedNodes[0], elem.connectedNodes[2] =  elem.connectedNodes[2], elem.connectedNodes[0]
                        elem.connectedNodes[1], elem.connectedNodes[3] =  elem.connectedNodes[3], elem.connectedNodes[1]
                        #
                        #print(elem.connectedNodes)

                        for j in range (i+1, len(elem.connectedNodes)):
                            #print('search for %d ' %elem.connectedNodes[i] , end='')
                            #print(elem.connectedNodes[j])
                            if (elem.connectedNodes[i] == elem.connectedNodes[j]):
                                simIdx = j
                                #print('simIdx: %d ' %simIdx)
                                break
                        #if (simIdx == -1):
                        #    elem.connectedNodes[0], elem.connectedNodes[1] =  elem.connectedNodes[1], elem.connectedNodes[0]

                    isFstNode = True
                    if (int(simIdx) % 2):
                        isFstNode = False

                    if (isFstNode):
                        #swapping linked nodes
                        #print ('1st swapping %d for %d and %d for %d' %(i+1,simIdx, i+2, simIdx+1))
                        elem.connectedNodes[i+1], elem.connectedNodes[simIdx] =  elem.connectedNodes[simIdx], elem.connectedNodes[i+1]
                        elem.connectedNodes[i+2], elem.connectedNodes[simIdx+1] =  elem.connectedNodes[simIdx+1], elem.connectedNodes[i+2]
                    else:
                        #swapping and reversing linked nodes
                        #print ('2nd swapping %d for %d and %d for %d' %(i+1,simIdx, i+2, simIdx-1))
                        elem.connectedNodes[i+1], elem.connectedNodes[simIdx] =  elem.connectedNodes[simIdx], elem.connectedNodes[i+1]
                        elem.connectedNodes[i+2], elem.connectedNodes[simIdx-1] =  elem.connectedNodes[simIdx-1], elem.connectedNodes[i+2]

                        if (np.abs(i+1  - simIdx) == 1 or np.abs(i  - simIdx)==1):
                            #print('swapping order')
                            elem.connectedNodes[i+1], elem.connectedNodes[simIdx] =  elem.connectedNodes[simIdx], elem.connectedNodes[i+1]

                    #print(elem.connectedNodes)
                    #print()
                    if (i == len( elem.connectedNodes )-3):
                        for j in range (1, len( elem.connectedNodes )-2, 2) :
                            if (elem.connectedNodes[j] != elem.connectedNodes[j+1]):
                                restart = True
            reorderOk = True
            for i in range (1, len( elem.connectedNodes )-2, 2) :
                if (elem.connectedNodes[i] != elem.connectedNodes[i+1]):
                    reorderOk = False
                    #print ('!!! %d Reorder not ok %s \n\n\n' %((i), elem.connectedNodes))
                    allReorderedFine = False
                    break
            #if (reorderOk == True):
            #    print ('Reordered fine: %s' %elem.connectedNodes)


    if (allReorderedFine == True):
        print('reordered fine...', end = '')
    else:
        print('NOT ALL REORDERED FINE...', end ='')
    print('done.')
    sys.stdout.flush()


    auxNodesInitLength = len (aux_nodes)
    #print(auxNodesInitLength)
    updatedElems = 0
    newAuxNodesA = []
    if (auxNodesInitLength != 0 and dim==3):
        print('Generating additional aux_nodes (if 0th and last node are aux nodes, creating another auxNode in a corner)...', end='')
        beamMidpoint = maxLim /2
        for elem in transportElements:
            #print('\n\n new element')
            nds = node_count
            #aux = len(aux_nodes)
            vertexA = vertices_out[int(elem.vertexA-node_count-aux)][0:dim]
            vertexB = vertices_out[int(elem.vertexB-node_count-aux)][0:dim]

            diffIdx = -1
            equalCoords = 0
            for d in range(3):
                if np.abs(vertexA[d]-vertexB[d])<1e-12:
                    #print(np.abs(vertexA[d]-vertexB[d]))
                    equalCoords +=1
                else:
                    diffIdx = d

            #if (len(elem.connectedNodes)==4 and elem.connectedNodes[len(elem.connectedNodes)-1]==elem.connectedNodes[0]):
            #    print (elem.connectedNodes)

            if (elem.connectedNodes[0]>=len(nodes_out) and elem.connectedNodes[len(elem.connectedNodes)-1]>=len(nodes_out) and
            (equalCoords == 1 or isTube==True) ):
                anodeA = np.asarray (aux_nodes[ int(elem.connectedNodes[0]-node_count) ][:])
                anodeB = np.asarray (aux_nodes[ int(elem.connectedNodes[len(elem.connectedNodes)-1]-node_count) ][:])

                #nanode = (anodeA + anodeB) /2

                #print(elem.getString())
                elem.addSingleConnectedNode( elem.connectedNodes[len(elem.connectedNodes)-1])
                #elem.addSingleConnectedNode( node_count + len(aux_nodes) )
                #elem.addSingleConnectedNode( node_count + len(aux_nodes) )saveTransportElements
                elem.addSingleConnectedNode( elem.connectedNodes[0] )

                #aux_nodes.append(nanode)
                #newAuxNodesA.append(nanode)
                #print(elem.getString())
                #a = input('').split(" ")[0]


            if (elem.connectedNodes[0]>=len(nodes_out) and elem.connectedNodes[len(elem.connectedNodes)-1]>=len(nodes_out) and equalCoords == 2  and isTube==False):
                #print('tube %s' %isTube)
                #print('Adding corner node')
                #print((elem.connectedNodes[0]/node_count))
                #print((elem.connectedNodes[len(elem.connectedNodes)-1]/node_count))
                updatedElems +=1
                #print(elem.connectedNodes)
                #print('old elem: %s' %elem.getStringyString(len(nodes_out), auxNodesInitLength))
                anodeA = np.asarray (aux_nodes[ int(elem.connectedNodes[0]-node_count) ][:])
                anodeB = np.asarray (aux_nodes[ int(elem.connectedNodes[len(elem.connectedNodes)-1]-node_count) ][:])
                #beamMidpoint = (anodeA + anodeB)/2
                #print( anodeA )
                #print( anodeB )
                nanode = np.zeros(3)

            #    print('vertexA %d array: %d' %(elem.vertexA, elem.vertexA-node_count-aux))
                vertexA = vertices_out[int(elem.vertexA-node_count-aux)][0:dim]
                vertexB = vertices_out[int(elem.vertexB-node_count-aux)][0:dim]

                #print(elem.connectedNodes)

                ridgeCoords = []
                lasti = -1
                for i in ((elem.connectedNodes)):
                    if (i!=lasti):
                        #print(i)
                        if (i<len(nodes_out)):
                            #print(nodes_out[int(i),0:3])
                            ridgeCoords.append(nodes_out[int(i),0:3])
                        if (i>=len(nodes_out)):
                            #print(aux_nodes[int(i-len(nodes_out))][0:3])
                            ridgeCoords.append(aux_nodes[int(i-len(nodes_out))][0:3])
                    lasti=i
                ridgeCoords = np.asarray(ridgeCoords)

                #face plane from first three points (nodes)
                pA = ridgeCoords[0]
                pB = ridgeCoords[1]
                pC = ridgeCoords[2]
                vecA = pB - pA
                vecB = pC - pA
                #planeN = np.dot(vecA,vecB)

                #lineDir = vertexB - vertexA

                #ndotu = np.dot(planeN, lineDir)
                #w = vertexB - pA
                #si = np.dot(-planeN, w) / ndotu
                #psi = w + si * lineDir + pA
                #print ('psi: %s' %psi)
                #print ('pa: %s' %pA)

                #print(int(elem.connectedNodes[0]/node_count))
                #print(int(elem.connectedNodes[len(elem.connectedNodes)-1]/node_count))
                elemMidpoint = (vertexA + vertexB) /2
                ridgeMidpoint = (anodeA + anodeB) /2
                #vector from elem midpoint to ridge midpoint
                vecV =  elemMidpoint - ridgeMidpoint


                #nanode = ridgeMidpoint + 1e-6
                #print('\nMidpoint aux node: %s' %nanode)
                nanode = elemMidpoint
                nanode[diffIdx] = anodeA[diffIdx]
                #print('New aux node: %s' %nanode)
                #
                #nanode = psi
                #adding new aux node to connected nodes
                #print('old elem: %s' %elem.connectedNodes)
                elem.addSingleConnectedNode( elem.connectedNodes[len(elem.connectedNodes)-1])
                elem.addSingleConnectedNode( node_count + len(aux_nodes) )
                elem.addSingleConnectedNode( node_count + len(aux_nodes) )
                elem.addSingleConnectedNode( elem.connectedNodes[0] )

                aux_nodes.append(nanode)
                newAuxNodesA.append(nanode)


                #print(ridgeCoords)
                ridgeCoords = []
                lasti = -1
                for i in ((elem.connectedNodes)):
                    if (i!=lasti):
                        #print(i)
                        if (i<len(nodes_out)):
                            #print(nodes_out[int(i),0:3])
                            ridgeCoords.append(nodes_out[int(i),0:3])
                        if (i>=len(nodes_out)):
                            #print(aux_nodes[int(i-len(nodes_out))][0:3])
                            ridgeCoords.append(aux_nodes[int(i-len(nodes_out))][0:3])
                    lasti=i
                ridgeCoords = np.asarray(ridgeCoords)
                if SHOW_PLOT:
                    fig = plt.figure()
                    ax = Axes3D(fig)
                    if AXIS_ASPECT_EQUAL:
                        ax.set_aspect('equal')
                    ax.plot3D([vertexA[0], vertexB[0]], [vertexA[1], vertexB[1]], [vertexA[2], vertexB[2]], marker='o')
                    #ax.plot3D([anodeA[0], anodeB[0]], [anodeA[1], anodeB[1]], [anodeA[2], anodeB[2]], marker='o')
                    ax.scatter3D(nanode[0], nanode[1], nanode[2])
                    #ax.scatter3D(anodeA[0], anodeA[1], anodeA[2])
                    #ax.scatter3D(anodeB[0], anodeB[1], anodeB[2])
                    #ax.scatter3D(ridgeMidpoint[0], ridgeMidpoint[1], ridgeMidpoint[2])
                    for r in range (len(ridgeCoords)-1):
                        ax.plot3D([ridgeCoords[r,0], ridgeCoords[r+1,0]], [ridgeCoords[r,1], ridgeCoords[r+1,1]], [ridgeCoords[r,2], ridgeCoords[r+1,2]], marker='x')
                    ax.plot3D([ridgeCoords[0,0], ridgeCoords[len(ridgeCoords)-1,0]], [ridgeCoords[0,1], ridgeCoords[len(ridgeCoords)-1,1]], [ridgeCoords[0,2], ridgeCoords[len(ridgeCoords)-1,2]], marker='x')

                    """
                    ax.plot3D([0,maxLim[0]], [0,0], [0,0], color='black')
                    ax.plot3D([0,maxLim[0]], [maxLim[1],maxLim[1]], [0,0], color='black')
                    ax.plot3D([0,maxLim[0]], [0,0], [maxLim[2],maxLim[2]], color='black')
                    ax.plot3D([0,maxLim[0]], [maxLim[1],maxLim[1]], [maxLim[2],maxLim[2]], color='black')

                    ax.plot3D([0,0], [0,maxLim[1]], [0,0], color='black')
                    ax.plot3D([0,0], [0,maxLim[1]], [maxLim[2],maxLim[2]], color='black')
                    ax.plot3D([0,0], [0,0], [0,maxLim[2]], color='black')
                    ax.plot3D([0,0], [maxLim[1],maxLim[1]], [0,maxLim[2]], color='black')

                    ax.plot3D([maxLim[0],maxLim[0]], [0,maxLim[1]], [0,0], color='black')
                    ax.plot3D([maxLim[0],maxLim[0]], [0,maxLim[1]], [maxLim[2],maxLim[2]], color='black')
                    ax.plot3D([maxLim[0],maxLim[0]], [0,0], [0,maxLim[2]], color='black')
                    ax.plot3D([maxLim[0],maxLim[0]], [maxLim[1],maxLim[1]], [0,maxLim[2]], color='black')
                    """
                    plt.show()

                    #print('new elem: %s' %elem.getString())
                    #print('new elem: %s' %elem.getStringyString(len(nodes_out), auxNodesInitLength))
                    #print(elem.connectedNodes)
                    #print()

            if (len(elem.connectedNodes)==2):
                print ('pruser')


        print('done.')
        sys.stdout.flush()

    print('Another renumbering of vertices...', end='')
    auxNodesFinalLength = len (aux_nodes)
    newAuxNodes = auxNodesFinalLength - auxNodesInitLength
    print('new aux nodes: %d' %newAuxNodes)
    for elem in transportElements:
        elem.vertexA += newAuxNodes
        elem.vertexB += newAuxNodes
    sys.stdout.flush()


    print('Coplanarity control before saving...', end='')
    wrongRidges = 0
    nodes_out = np.asarray(nodes_out)
    aux_nodes = np.asarray(aux_nodes)
    #print(ndCrds[0])
    #print(ndCrds[0, 0:dim])
    for element in transportElements:

        nodesCoords = []
        #print(element.connectedNodes)
        for nIdx in range(0,len(element.connectedNodes),2):
            nIdx = int(element.connectedNodes[nIdx])

            if (nIdx<len(nodes_out)):
                nodesCoords.append(nodes_out[nIdx, 0:dim])
            elif (nIdx<len(nodes_out)+auxNodesInitLength):
                nodesCoords.append(aux_nodes[nIdx-node_count, 0:dim])
            elif (nIdx<len(nodes_out)+auxNodesFinalLength):
                nodesCoords.append(aux_nodes[nIdx-node_count, 0:dim])
            else:
                print('Unknown node nr %d ' %nIdx, end='')
                print('Max idx: %d ' %(len(nodes_out)+len(aux_nodes)-1), end='')

        nodesCoords = np.asarray(nodesCoords)
        allCoplanar, val = checkCoplanarity(nodesCoords, 1e-15)

        if ( not allCoplanar):
            wrongRidges +=1
            print('Wrong element!')
            print(element.getStringyString(len(nodes_out), auxNodesInitLength))
            print('Not coplanar!!! Ridge nr. %d, err: %.10e' %(transportElements.index(element), val ))
            print(element.getReducedString())

            vrtcs = []
            ndss = []
            for n in element.connectedNodes:
                if (n<len(nodes_out)):
                    print(nodes_out[int(n),0:3])
                    vrtcs.append(nodes_out[int(n),0:3])
                if (n>=len(nodes_out)):
                    vrtcs.append(aux_nodes[int(n-len(nodes_out)), 0:3])
            vrtcs = np.asarray(vrtcs)

            vo = np.asarray(vertices_out)
            ndss.append(vo[int(element.vertexA - len(nodes_out)-len(aux_nodes)),0:3])
            ndss.append(vo[int(element.vertexB - len(nodes_out)-len(aux_nodes)),0:3])
            ndss = np.asarray(ndss)

            if SHOW_PLOT:
                fig = plt.figure()
                ax = Axes3D(fig)
                #ax.auto_scale_xyz([0, maxLim[0]], [0, maxLim[1]], [0, maxLim[2]])
                ax.scatter3D(vrtcs[:,0], vrtcs[:,1], vrtcs[:,2])
                ax.scatter3D(ndss[:,0], ndss[:,1], ndss[:,2])

                newAuxNodesA =np.asarray(newAuxNodesA)
                ax.scatter3D(newAuxNodesA[:,0], newAuxNodesA[:,1], newAuxNodesA[:,2])
                #ax.plot(ndss[0,:], ndss[1,:])

                ax.plot3D([0,0,0], [maxLim[0],0,0])
                ax.plot3D([0,maxLim[1],0], [maxLim[0],maxLim[1],0])
                #ax.plot(np.array([0,0,maxLim[2]]), np.array([maxLim[0],0,maxLim[2]]))
                #ax.plot(np.array([0,maxLim[1],maxLim[2]]), np.array([maxLim[0],maxLim[1],maxLim[2]]))

                ax.scatter(vrtcs[len(vrtcs)-1,0], vrtcs[len(vrtcs)-1,1], vrtcs[len(vrtcs)-1,2])
                plt.show()


            allCoplanar = False



    print('done. ', end='')
    #print('Updated elems: %d' %updatedElems)
    print('Wrong elems: %d' %wrongRidges)


    if len(mZ)>0:
        print ('Material zones detected for transport')
        if dim == 2:
            for elem in transportElements:
                if mZ[0][0] =='circle':
                    inRebar = False
                    vo = np.asarray(vertices_out)
                    vertA = vo[int(elem.vertexA - len(nodes_out)-len(aux_nodes)),0:2]
                    vertB = vo[int(elem.vertexB - len(nodes_out)-len(aux_nodes)),0:2]

                    for rebar in range (mZ[0][3]):
                        if ((np.linalg.norm(vertA-mZ[rebar][2]) < mZ[rebar][1] ) and
                        (np.linalg.norm(vertB-mZ[rebar][2]) < mZ[rebar][1] )):
                            inRebar = True

                    if inRebar:
                        #print('inrebar')
                        elem.material = 3
                        #print(elem.material)





    print('Saving TRSPRT elements...', end='')
    sys.stdout.flush()
    print(coupled)
    print('Trsprt elements: %d' %len(transportElements))
    with open(os.path.join(master_folder,trsprtElemsFile), 'w') as f:
        headerLine = '#ElemType\tvrtxAIdx\tvrtxBIdx\tnrOfNodes\tnodesIdx\tMaterial'
        f.write("%s\n" % headerLine )
        for element in transportElements:
            #print ("%s\n" % element.getString() )
            if (dim==2):f.write("%s\n" % element.getString(coupled=coupled) )
            if (dim==3): f.write("%s\n" % element.getReducedString(coupled=coupled) )

    return newAuxNodes

def saveRigidPlates(master_folder, dim, rigidPlates, totalNodeCount, trspt=False, expansionRingsProps=[]):
    print('Saving rigid plates...', end='')

    if trspt == False:
        file = constraintFile
    else:
        file = constraintTrsptFile

    with open(os.path.join(master_folder,file), 'w') as f:
        totNodeC = totalNodeCount
        headerLine = '#ConstraintType\tGovNodeIdx\tXmin\tXmax\tYmin\tYmax'
        if (dim==3):
            headerLine += '\tZmin\tZmax'
        f.write("%s\n" % headerLine )
        for i in range(len(rigidPlates)):
            rplt = rigidPlates[i]
            rplt.govNodeIdx = totalNodeCount + i

            f.write("%s\n" % rplt.getString() )

            totNodeC += 1

        if len(expansionRingsProps)>0 and dim == 2:
            rebarCount = expansionRingsProps[0]
            rebarDepth = expansionRingsProps[1]
            rebarDiameter = expansionRingsProps[2]
            maxLim = expansionRingsProps[3]

            print ('rebarcount %d' %rebarCount)

            for r in range (rebarCount):
                if (rebarCount==1):
                    #puvodni poloha rebars polovina od kraje
                    centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
                else:
                    #poloha rebars presne jak je ve clanku
                    centre = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth  ])
                #f.write( 'ExpansionRing	%d %e %e %e %e volExpFn %d\n' %(totNodeC+r, centre[0],centre[1],  0, rebarDiameter/2, 2 ) )
                #f.write( 'ExpansionRingDoFLoad %d %e %e %e %e expansionMaster %d\n' %(totNodeC+r, centre[0],centre[1],  0, rebarDiameter/2, totNodeC+rebarCount ) )
                f.write( 'ExpansionRingSingleDoFLoad %d %e %e 0 %e \n' %( totNodeC+r, centre[0],centre[1], rebarDiameter/2*1.01 ) )



    print('done.')


## ForceGauge file_name gauge_name which_force how_many node_ids
#ForceGauge LD reactF fx 1 5678
def saveForceGauge(master_folder, dir,columnName, nodeIdx):
    print('Saving force %s gauge for node %d' %(dir, nodeIdx) )
    fl=open(os.path.join(master_folder,exportersFile),'a')
    fl.write('ForceGauge LD %s %s 1 %d\n' %(dir, columnName, nodeIdx))
    fl.close()

### displacement gauge uy x1coord y1coord x2coord y2coord
#DisplacementGauge LD phiX rotx 0 0 0 1.0 0 0
#DisplacementGauge LD uX ux 0 0 0 1.0 0 0
def saveDisplacementGauge(master_folder, columnName, dir, coordsA, coordsB):
    print('Saving displacement %s gauge between point %s and %s' %(dir, coordsA, coordsB) )
    fl=open(os.path.join(master_folder,exportersFile),'a')
    if (len(coordsA)==2):
        fl.write('DisplacementGauge LD %s\t%s\t%e\t%e\t%e\t%e\n' %(columnName, dir, coordsA[0],coordsA[1],coordsB[0],coordsB[1]))
    if (len(coordsA)==3):
        fl.write('DisplacementGauge LD %s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\n' %(columnName, dir, coordsA[0],coordsA[1], coordsA[2],coordsB[0],coordsB[1], coordsB[2]))
    fl.close()

def saveConstraint(master_folder, dim, govNodes, govNodesMechBC, rigidPlates, totalNodeCount, nodes, expansionRingsProps=[],virtualDoF=0, nodesMechBC = []):
    #saving gov nodes
    saveNodes (master_folder,govNodes, "GovParticle", dim, govNodesFile, virtualDoF=virtualDoF)
    #saving gov nodes mech BC
    for i in range (len(govNodesMechBC)):
        m = govNodesMechBC[i]
        m.nodeIdx = totalNodeCount + i
    saveMechBC(master_folder,dim, govNodesMechBC, govNodesBC = True)

    #saving force gauges for rigid plates
    for i in range (len(govNodesMechBC)):
        if virtualDoF == 0:
            saveForceGauges(master_folder, dim, govNodesMechBC[i].nodeIdx, name='MechPLT%d'%i)
        else:
            if i < len(govNodesMechBC)-virtualDoF :
                saveForceGauges(master_folder, dim, govNodesMechBC[i].nodeIdx, name='MechPLT%d'%i)
            else:
                saveForceGauges(master_folder, dim, govNodesMechBC[i].nodeIdx, name='VirtualDOF%d'%i, virtualDoF=virtualDoF)


    #saving rigid plates
    for rp in rigidPlates:
        rp.getNodesAffected(nodes)
    saveRigidPlates(master_folder, dim, rigidPlates, totalNodeCount, expansionRingsProps=expansionRingsProps)

    for rp in range (len(rigidPlates)):
        print('Nodes affected by Rigid plate #%d:' %rp)
        print(rigidPlates[rp].getNodesAffected(nodes))

    print('mechbc')
    #if (len(nodesMechBC)!=0):
    #    for i in range(len(nodesMechBC)):
    #s        saveForceGauges(master_folder, dim, nodesMechBC[i].nodeIdx, name='Node%d'%nodesMechBC[i].nodeIdx, moments=False)


def saveConstraintTransport(master_folder, dimension, govNodesTrspt, govNodesTrsptBC, rigidPlatesTrspt, totalNodeCount, node_coords, vert_count, verticesIdxDict, vertIdxStart):
    print ('Saving Transport constraint...')
    print(govNodesTrspt)
    saveNodes (master_folder,govNodesTrspt, "TrsprtNode", dimension, govNodesTrsptFile)

    for i in range (len(govNodesTrsptBC)):
        m = govNodesTrsptBC[i]
        m.nodeIdx = totalNodeCount + i


    for rp in rigidPlatesTrspt:
        rp.renumberVertices(verticesIdxDict, vertIdxStart)


    #print (govNodesTrsptBC)
    saveTransportBC(master_folder,govNodesTrsptBC, verticesIdxDict, vertIdxStart, govNodesBC=True, totalNodeCount=totalNodeCount)


    for i in range (len(govNodesTrsptBC)):
        saveForceGauges(master_folder, dimension, govNodesTrsptBC[i].nodeIdx, moments=False, name='TrsptPLT%d'%i, transport=True)


    saveRigidPlates(master_folder, dimension, rigidPlatesTrspt, totalNodeCount, trspt=True)


def saveForceGauges(master_folder, dimension, nodeIdx, moments=True, name='', transport = False, virtualDoF=-1):
    if transport == False:
        if (name==''):
            saveForceGauge(master_folder, 'fx#%d'%nodeIdx , 'fx', nodeIdx )
            saveForceGauge(master_folder, 'fy#%d'%nodeIdx, 'fy', nodeIdx )
            if (dimension==3): saveForceGauge(master_folder, 'fz#%d'%nodeIdx, 'fz', nodeIdx )
            if moments == True:
                if (dimension==3): saveForceGauge(master_folder, 'mx#%d'%nodeIdx , 'mx', nodeIdx )
                if (dimension==3): saveForceGauge(master_folder, 'my#%d'%nodeIdx, 'my', nodeIdx )
                saveForceGauge(master_folder, 'mz#%d'%nodeIdx, 'mz', nodeIdx )
        else:
            saveForceGauge(master_folder, 'fx#%s'%name , 'fx', nodeIdx )
            saveForceGauge(master_folder, 'fy#%s'%name, 'fy', nodeIdx )
            if (dimension==3): saveForceGauge(master_folder, 'fz#%s'%name, 'fz', nodeIdx )
            if moments == True:
                if (dimension==3): saveForceGauge(master_folder, 'mx#%s'%name , 'mx', nodeIdx )
                if (dimension==3): saveForceGauge(master_folder, 'my#%s'%name, 'my', nodeIdx )
                if not (virtualDoF>0):
                    saveForceGauge(master_folder, 'mz#%s'%name, 'mz', nodeIdx )
    else:
        if (name==''):
            #saveForceGauge(master_folder, 'fx#%d'%nodeIdx , 'fx', nodeIdx )
            saveForceGauge(master_folder, 'flux#%d'%nodeIdx, 'flux', nodeIdx )

        else:
            #saveForceGauge(master_folder, 'fx#%s'%name , 'fx', nodeIdx )
            saveForceGauge(master_folder, 'flux#%s'%name, 'flux', nodeIdx )



#saveDisplacementGauge(master_folder, columnName, dir, coordsA, coordsB):
def saveDisplacementGauges(master_folder, dimension, name, coordsA, coordsB, rotations = False):
    saveDisplacementGauge(master_folder, 'ux#%s'%name, 'ux', coordsA, coordsB )
    saveDisplacementGauge(master_folder, 'uy#%s'%name, 'uy', coordsA, coordsB )
    if (dimension==3): saveDisplacementGauge(master_folder, 'uz#%s'%name, 'uz', coordsA, coordsB )
    if rotations == True:
        if (dimension==3): saveDisplacementGauge(master_folder, 'rotx#%s'%name , 'rotx', coordsA, coordsB )
        if (dimension==3): saveDisplacementGauge(master_folder, 'roty#%s'%name, 'roty', coordsA, coordsB )
        saveDisplacementGauge(master_folder, 'rotz#%s'%name, 'rotz', coordsA, coordsB )


def saveMeasuringGauges(master_folder, dimension, measuringGauges):
    print('Saving measuring gauges...')
    print ('%d gauges' %len(measuringGauges))
    for mg in measuringGauges:
        saveDisplacementGauges(master_folder, dimension, mg.name, mg.coordsA, mg.coordsB, rotations = mg.rotation)


def checkSavedModel(master_folder, dim, activeMechanics, activeTransport):
    print('\nChecking  generated files...')
    allOK = True
    cols = np.arange(1,dim+1)

    print('Loading back node coords...', end='')
    test_nodeCoords = np.genfromtxt(os.path.join(master_folder,nodesFile),  dtype= None, encoding='ascii', usecols=cols)
    print('\t\t %d nodes loaded.' %len(test_nodeCoords))


    if (os.path.exists(os.path.join(master_folder,auxNodesFile))):
        print('Loading back aux node coords...', end='')
        test_auxNodeCoords = np.genfromtxt(os.path.join(master_folder,auxNodesFile),  dtype= None, encoding='ascii', usecols=cols)
        print('\t\t %d nodes loaded.' %len(test_auxNodeCoords))
    else:
        test_auxNodeCoords = np.zeros(3)

    print('Loading back vertices coords...', end='')
    test_verticesCoords = np.genfromtxt(os.path.join(master_folder,verticesFile),  dtype= None, encoding='ascii', usecols=cols)
    print('\t\t %d vertices loaded.' %len(test_verticesCoords))

    test_solverNodeArray = np.vstack((test_nodeCoords, test_auxNodeCoords, test_verticesCoords))

    if (activeMechanics):
        print('Loading back mechanical elements...', end='')
        test_mechElems = []
        with open(os.path.join(master_folder,mechElemsFile)) as f:
            for line in f:
                test_mechElems.append(line.split())
        test_mechElems.pop(0)
        print('\t %d MechElems loaded.' %len(test_mechElems))

        print('Reassembling MechElems, checking face coplanarity... ')
        maxErr = 0
        wrongElems = 0
        for mechElem in test_mechElems:
            #LTCBEAM	200	75	4	449	450	458	451	0	0
            name = mechElem[0]
            nA = test_solverNodeArray[int(mechElem[1])]
            nB = test_solverNodeArray[int(mechElem[2])]
            verticesNr = int(mechElem[3])
            vertices = []
            for v in range (verticesNr):
                vertices.append(test_solverNodeArray [int(mechElem[4+v])] )
            material =  int (mechElem[4+verticesNr])

            #checking coplanarity
            allCoplanar, val = checkCoplanarity(vertices, 1e-15)
            if (val > maxErr): maxErr = val
            if (allCoplanar == False):
                wrongElems +=1
                #a = input('').split(" ")[0]
                #print(mechElem)
                #print(val)
                #print(mechElem[3:3+verticesNr])
                #print(vertices)
        if (wrongElems==0):
            print('All faces coplanar. Mech Elems OK. MaxErr: %e' %maxErr)
        else:
            print ('Wrong faces: %d/%d !!!!' %(wrongElems, len(test_mechElems)))
            print('MaxErr: %e' %maxErr)
            allOK = False


    if (activeTransport):
        print('Loading back transport elements...', end='')
        test_trsprtElems = []
        with open(os.path.join(master_folder,trsprtElemsFile)) as f:
            for line in f:
                test_trsprtElems.append(line.split())
        test_trsprtElems.pop(0)
        print('\t %d TrsprtElems loaded.' %len(test_trsprtElems))

        print('Reassembling TrsprtElems, checking face coplanarity... ')
        maxErr = 0
        wrongElems = 0
        for trsprtElem in test_trsprtElems:
            #LTCTRSP	436	435	3	213	154	196	1
            #print(trsprtElem)
            name = trsprtElem[0]
            nA = test_solverNodeArray[int(trsprtElem[1])]
            nB = test_solverNodeArray[int(trsprtElem[2])]
            verticesNr = int(trsprtElem[3])
            vertices = []
            for v in range (verticesNr):
                vertices.append(test_solverNodeArray [int(trsprtElem[4+v])] )
            material =  int (trsprtElem[4+verticesNr])

            #checking coplanarity
            allCoplanar, val = checkCoplanarity(vertices, 1e-8)
            if (val > maxErr): maxErr = val
            if (allCoplanar == False):
                wrongElems +=1
                #a = input('').split(" ")[0]
                #print(mechElem)
                #print(mechElem[3:3+verticesNr])
                #print(vertices)
        if (wrongElems==0):
            print('All faces coplanar. Trsprt Elems OK. MaxErr: %e' %maxErr)
        else:
            print ('Wrong faces: %d !!!!' %wrongElems)
            allOK = False

        if (allOK == True):
            print('Model seems ok. All fine.')
        else:
            print('!!!!! Saved with problems !!!!!')


def checkCoplanarity(points, maxError):
    allCoplanar = True
    nodesCoords = np.asarray(points)
    #print('checking face: %s' %nodesCoords)
    val = 0
    for v in range ( len(nodesCoords)-3 ):
        val = 0
        pA = nodesCoords [v]
        pB = nodesCoords [v+1]
        pC = nodesCoords [v+2]
        pD = nodesCoords [v+3]

        val = equation_plane(pA, pB, pC, pD)
        if (np.abs(val) > maxError):
            #print(val)
            allCoplanar = False

    return allCoplanar, np.abs(val)



if __name__ == '__main__':
    import doctest
    doctest.testmod()







#
