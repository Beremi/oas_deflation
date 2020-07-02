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
constraintFile              = "constraint.inp"
solverFile                  = "solver.inp"

#
#coplanarity test
def equation_plane(pA, pB, pC, pD):
    """
    pA = np.array([0,1,0])
    pB = np.array([1,1,0])
    pC = np.array([0,1,4])
    pD = np.array([0,0,0])
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
    cos = dot(vecA,vecB)/norm(vecA)/norm(vecB)
    angle = arccos(clip(c, -1, 1))
    print ('Angle is %f' %(ang))


#check if any number in matrix is lower than
def checkLowerThan (matrix, minDist):
    return   all(i >= minDist for i in matrix)



#check mutual distances between particles using cDist
def checkMutDistancesCdist(dim, minDist, currentNodes, newNode):
    ncrds = np.asarray(currentNodes)
    crds = np.asarray(newNode)
    crds = np.reshape(crds, (-1, dim))
    dists = scipy.spatial.distance.cdist(crds, ncrds , 'euclidean')
    dists = dists.flatten()
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
    print('''Using Python version of generator. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')



#check mutual distances between particles using loops in a periodic domain
def checkMutDistancesLoopsPeriodic (dim, minDist, currentNodes, newNode, maxLim):
    distIsGood = True
    return distIsGood

try:
    from point_generators_cython import checkMutDistancesLoopsPeriodic_cython as checkMutDistancesLoopsPeriodic
    print('Using Cython version of point generator - checkMutDistancesLoopsPeriodic.')
except:
    print('''Using Python version of generator. To use the Cython version the
          the code has to be build using: python setup.py build_ext --inplace.''')


def extractGeometry (master_folder, dim, node_count, maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=None, periodicModel = 0, nodePositions = None, coupledNodes = None, mirtype = None, notches = None):
    if (dim == 2):
        if (periodicModel == 0):
            vert_count, verticesIdxDict, vertIdxStart, totalNodeCount = output2D(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=mZ, notches = notches)
        if (periodicModel == 1):
            vert_count, verticesIdxDict, vertIdxStart, totalNodeCount = output2DPeriodic(master_folder, node_count,  maxLim, vor, node_coords, areas, nodePositions, coupledNodes, mirtype, mZ=mZ )
    if (dim == 3):
        vert_count, verticesIdxDict, vertIdxStart,totalNodeCount = output3D(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=mZ,  notches = notches)
    return vert_count, verticesIdxDict, vertIdxStart, totalNodeCount


#Extract geometry 2d
def output2D(master_folder, node_count,  maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=None, notches = None):
    dim = 2
    print('Extracting the geometry...', end='')
    sys.stdout.flush()
    nodes_out = np.zeros( (node_count, (2 + 1 + 1 )))
    nodes_out[:, 0:2] = vor.points[0:node_count , 0:2]
    nodes_out[:, dim] = 0
    nodes_out[:, dim + 1] = areas[:]

    #relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
    #print ('Area error: %.5e ' %(relAreaError)  )
    ########################################################################################
    cond = np.any((vor.ridge_points < node_count) & (vor.ridge_points >= 0), axis=1)
    validRidgeIdxs = np.where(cond)[0]
    #print(validRidgeIdxs.shape)
    #print(validRidgeIdxs)

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
    for i in range (validRidgeIdxs.size):
        #array for two vertices A and B
        vrtxA = np.zeros ( (dim  +1 +1 ) )
        vrtxB = np.zeros ( (dim  +1 +1 ) )

        #original indices of vertices A and B
        vertA = vor.ridge_vertices[validRidgeIdxs[i]][0]
        vertB = vor.ridge_vertices[validRidgeIdxs[i]][1]

        # copying of coordinates of vertices A and B
        for d in range (dim):
            vrtxA [d] = vor.vertices[vertA][d]
            vrtxB [d] = vor.vertices[vertB][d]

        #copying of original indices of vertices A and B
        vrtxA[dim] = vertA
        vrtxB[dim] = vertB

        vrtxA[dim+1] = 0
        vrtxB[dim+1] = 0

        #duplicity check
        addVrtxA = True
        addVrtxB = True
        for j in range (len(vertices_out)):
            if (vertices_out[j][0] ==  vor.vertices[vertA][0] and vertices_out[j][1] ==  vor.vertices[vertA][1]):
                addVrtxA = False
            if (vertices_out[j][0] ==  vor.vertices[vertB][0] and vertices_out[j][1] ==  vor.vertices[vertB][1]):
                addVrtxB = False

        #adding the vertices into the list of vertices if new
        if (addVrtxA == True):
            verticesIdxDict.update( { vertA : len(vertices_out)  } )
            vrtxA [dim] = len(vertices_out)
            vertices_out.append(vrtxA)

        if (addVrtxB == True):
            verticesIdxDict.update( { vertB : len(vertices_out)  } )
            vrtxB [dim] = len(vertices_out)
            vertices_out.append(vrtxB)

        #ridges
        ########################################################
        #Array for ridge nAidx, nBidx, nrVrt, vertAidx, vertBidx
        rdg = np.zeros ( (2 + 1 +  2) )

        #indices of two nodes that are divided by the ridge
        pointA = vor.ridge_points[validRidgeIdxs[i],0]
        pointB = vor.ridge_points[validRidgeIdxs[i],1]

        #creating auxiliary nodes if one of nodes is outside
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
    v_count = len (vertices_out)
    vertIdxStart = node_count + len(aux_nodes)

    for i in range (len(ridges_out)):
        ridges_out[i][3] += vertIdxStart
        ridges_out[i][4] += vertIdxStart

    print('done.')
    sys.stdout.flush()
    #output: nodes_out, aux_nodes, vertices_out, ridges_out

    saveNodes(master_folder, aux_nodes, "AuxNode",dim, auxNodesFile)
    if activeMechanics:
        saveNodes(master_folder, nodes_out, "Particle",dim, nodesFile)
        saveMechanicalElements(master_folder, ridges_out, node_count, dim, nodes_out, mZ=mZ, notches = notches)
    else:
        saveNodes(master_folder, nodes_out, "AuxNode",dim, nodesFile)
    if activeTransport:
        saveNodes(master_folder, vertices_out, "TrsprtNode",dim, verticesFile)
        saveTransportElements(master_folder, ridges_out,dim, node_count, v_count, aux_nodes, maxLim, nodes_out, vertices_out)
    else:
        saveNodes(master_folder, vertices_out, "AuxNode",dim, verticesFile)

    totalPointCount = len(nodes_out) + len(aux_nodes) + len(vertices_out)

    return v_count, verticesIdxDict, vertIdxStart, totalPointCount #, nodes_out, aux_nodes, vertices_out, ridges_out



#Extract geometry 2d periodic torus
def output2DPeriodic(master_folder, node_count,  maxLim, vor, node_coords, areas, nodePositions, coupledNodes, mirtype, mZ=None):
    dim = 2



    print('Filtering valid ridges of 2d periodic model...')
    print ('Filtering valid ridges...', end ='')
    # arrays for nodes of the resulting model
    valid_ridges = np.empty((0,1)).astype(int)
    valid_ridge_nodes = np.empty((0,2)).astype(int)
    valid_ridge_vertices = np.empty((0,2)).astype(int)

    """
    # selecting points
    for ir,r in enumerate(vor.ridge_points):
        if ((nodePositions[r[0]]>0 and nodePositions[r[1]]>0)):
            plt.plot( [vor.points[ r[0],0 ] , vor.points[r[1],0 ]], [vor.points[ r[0],1 ] , vor.points[r[1],1 ]] , 'ro-', color='green', alpha = 0.5)
            plt.text(vor.points[ r[0],0 ] , vor.points[ r[0],1 ] , nodePositions[r[0]], fontsize=11)
            plt.text(vor.points[ r[1],0 ] , vor.points[ r[1],1 ] , nodePositions[r[1]], fontsize=11)

            valid_ridges = np.vstack((valid_ridges, ir))
            valid_ridge_nodes = np.vstack((valid_ridge_nodes, r))
            valid_ridge_vertices = np.vstack((valid_ridge_vertices, vor.ridge_vertices[ir]))


        if ((nodePositions[r[0]]*nodePositions[r[1]]<0)):
            plt.plot( [vor.points[ r[0],0 ] , vor.points[r[1],0 ]], [vor.points[ r[0],1 ] , vor.points[r[1],1 ]] , 'ro-', color='red', alpha = 0.5)
            plt.text(vor.points[ r[0],0 ], vor.points[ r[0],1 ] , nodePositions[r[0]], fontsize=11)
            plt.text(vor.points[ r[1],0 ], vor.points[ r[1],1 ] , nodePositions[r[1]], fontsize=11)

            valid_ridges = np.vstack((valid_ridges, ir))
            valid_ridge_nodes = np.vstack((valid_ridge_nodes, r))
            valid_ridge_vertices = np.vstack((valid_ridge_vertices, vor.ridge_vertices[ir]))


    plt.show()
    print('done.')
    #"""


    #if (mirtype[r[0]]==0 and mirtype[r[1]]>=0) or (mirtype[r[1]]==0 and mirtype[r[0]]>=0) or mirtype[r[0]]*mirtype[r[1]]==2:
    # selecting points
    for ir,r in enumerate(vor.ridge_points):
        if ((mirtype[r[1]]>=0 and mirtype[r[0]]==0) or (mirtype[r[1]]==0 and mirtype[r[0]]>=0)):
            #plt.plot( [vor.points[ r[0],0 ] , vor.points[r[1],0 ]], [vor.points[ r[0],1 ] , vor.points[r[1],1 ]] , 'ro-', color='green', alpha = 0.5)
            #plt.text(vor.points[ r[0],0 ] , vor.points[ r[0],1 ] , nodePositions[r[0]], fontsize=11)
            #plt.text(vor.points[ r[1],0 ] , vor.points[ r[1],1 ] , nodePositions[r[1]], fontsize=11)

            valid_ridges = np.vstack((valid_ridges, ir))
            valid_ridge_nodes = np.vstack((valid_ridge_nodes, r))
            valid_ridge_vertices = np.vstack((valid_ridge_vertices, vor.ridge_vertices[ir]))

        if ( mirtype[r[0]]*mirtype[r[1]]==2 ):
            #plt.plot( [vor.points[ r[0],0 ] , vor.points[r[1],0 ]], [vor.points[ r[0],1 ] , vor.points[r[1],1 ]] , 'ro-', color='red', alpha = 0.5)
            #plt.text(vor.points[ r[0],0 ], vor.points[ r[0],1 ] , nodePositions[r[0]], fontsize=11)
            #plt.text(vor.points[ r[1],0 ], vor.points[ r[1],1 ] , nodePositions[r[1]], fontsize=11)

            valid_ridges = np.vstack((valid_ridges, ir))
            valid_ridge_nodes = np.vstack((valid_ridge_nodes, r))
            valid_ridge_vertices = np.vstack((valid_ridge_vertices, vor.ridge_vertices[ir]))

    #plt.show()
    #"""
    print('done.')


    #input("Press Enter to continue...")

    print ('Renumbering nodes...', end = '')
    valid_ridges = valid_ridges.flatten()
    valid_ridge_nodes = valid_ridge_nodes.flatten()
    valid_ridge_nodes = np.unique(valid_ridge_nodes)
    validNodesDict = {}
    for i in range (len(valid_ridge_nodes)):
        validNodesDict.update( {  int(valid_ridge_nodes[i]) : int(i)  } )

    node_count = len(valid_ridge_nodes)
    print('done.')


    cpldNds = []
    #coupledNodes = []
    for i in range (len(coupledNodes)):
        if ( coupledNodes[i][0] in valid_ridge_nodes and coupledNodes[i][1] in valid_ridge_nodes ):
            #plt.plot( [vor.points[ coupledNodes[i][0],0 ] , vor.points[ coupledNodes[i][1],0 ]], [vor.points[ coupledNodes[i][0],1 ] , vor.points[ coupledNodes[i][1],1 ]] ,'ro-', color='red')
            #plt.text(vor.points[ coupledNodes[i][0],0 ] , vor.points[ coupledNodes[i][0],1 ]  , nodePositions[coupledNodes[i][0]], fontsize=11)
            #plt.text(vor.points[ coupledNodes[i][1],0 ] , vor.points[ coupledNodes[i][1],1 ]  , nodePositions[coupledNodes[i][1]], fontsize=11)

            cpldNds.append ( np.array( [ validNodesDict[int(coupledNodes[i][0])]    ,   validNodesDict[int(coupledNodes[i][1])]     ]    )   )

    #plt.show()
    print ('done.')
    print('filtering done.')

    print('Extracting the geometry of 2d periodic torus...')
    sys.stdout.flush()
    nodes_out = np.zeros( (node_count, (2 + 1 + 1 )))

    for i in range (node_count):
        nodes_out [i, 0:2] = vor.points[ valid_ridge_nodes[i], : ]

    """
    print ('Node valid nr %s' %nd)
    print ('in vor idx %s' % valid_ridge_nodes[nd])
    print ('in out dict %s' %validNodesDict[valid_ridge_nodes[nd]])

    print ('crds orig %s' %vor.points[ valid_ridge_nodes[nd] ,: ])
    print ('crds save %s' %nodes_out[ validNodesDict[valid_ridge_nodes[nd]] ,:] )

    #plt.plot(vor.points[:,0], vor.points[:,1], 'x', color='red')
    #plt.plot(nodes_out[:,0], nodes_out[:,1], 'o', color='black')
    #plt.plot(vor.points[0,0], vor.points[0,1], '*', color='green')

    #print (len(valid_ridge_nodes))
    #print (len(nodes_out))
    #print (len(vor.points))
    #"""

    nodes_out[:, dim] = 0
    nodes_out[:, dim + 1] = 0

    ########################################################################################
    validRidgeIdxs = []
    for i in range (vor.ridge_points.shape[0]):
        pr = False
        if ( vor.ridge_points[i][0] in validNodesDict and vor.ridge_points[i][1] in validNodesDict):
                pr=True
        if (pr):
           validRidgeIdxs.append(i)
          # print (vor.points[vor.ridge_points[i][p]])


    validRidgeIdxs = np.asarray(validRidgeIdxs)
    ########################################################################################
    # vertices: [xA,yA,zA] [origIdx]
    #vertices_out = []
    # dictionary of original and new indices of vertices
    verticesIdxDict = {}
    # ridges: nodeAidx, nodeBidx, trsprtBC, vertIdx
    ridges_out = []
    #auxiliary nodes
    aux_nodes = []

    #vertices
    vertices_out = []
    ####################################################
    #list of vertices, list of beams
    for i in valid_ridges:
        #input("Press Enter to continue...")
        #print()
        #print(i)
    #for i in range (validRidgeIdxs.size):
        #array for two vertices A and B
        vrtxA = np.zeros ( (dim  +1 +1 ) )
        vrtxB = np.zeros ( (dim  +1 +1 ) )

        #original indices of vertices A and B
        vertA = vor.ridge_vertices[i][0]
        vertB = vor.ridge_vertices[i][1]

        # copying of coordinates of vertices A and B
        for d in range (dim):
            vrtxA [d] = vor.vertices[vertA][d]
            vrtxB [d] = vor.vertices[vertB][d]

        #copying of original indices of vertices A and B
        vrtxA[dim] = vertA
        vrtxB[dim] = vertB

        vrtxA[dim+1] = 0
        vrtxB[dim+1] = 0

        #duplicity check
        addVrtxA = True
        addVrtxB = True
        for j in range (len(vertices_out)):
            if (vertices_out[j][0] ==  vor.vertices[vertA][0] and vertices_out[j][1] ==  vor.vertices[vertA][1]):
                addVrtxA = False
            if (vertices_out[j][0] ==  vor.vertices[vertB][0] and vertices_out[j][1] ==  vor.vertices[vertB][1]):
                addVrtxB = False

        #adding the vertices into the list of vertices if new
        if (addVrtxA == True):
            verticesIdxDict.update( { vertA : len(vertices_out)  } )
            vrtxA [dim] = len(vertices_out)
            vertices_out.append(vrtxA)

        if (addVrtxB == True):
            verticesIdxDict.update( { vertB : len(vertices_out)  } )
            vrtxB [dim] = len(vertices_out)
            vertices_out.append(vrtxB)
        #"""
        #ridges
        ########################################################
        #Array for ridge nAidx, nBidx, nrVrt, vertAidx, vertBidx
        rdg = np.zeros ( (2 + 1 +  2) )

        #indices of two nodes that are divided by the ridge
        pointA = vor.ridge_points[i,0]
        pointB = vor.ridge_points[i,1]

        #
        #new_validNodesDict = dict([(value, key) for key, value in validNodesDict.items()])
        #
        #print (new_validNodesDict)
        key = pointA
        #print ('node %d to dict %d' %(key, validNodesDict[key]))
        rdg[0] = int(validNodesDict[key])
        key = pointB
        #print ('node %d to dict %d' %(key, validNodesDict[key]) )
        rdg[1] = int(validNodesDict[key])

        #number of vertices
        rdg[2] = 2

        #indices of vertices
        key = vertA
        #print ('vert %d to dict %d' %(key, verticesIdxDict[key]) )
        rdg[3] = verticesIdxDict[vertA]
        key = vertB
        #print ('vert %d to dict %d' %(key, verticesIdxDict[key]) )
        rdg[4] = verticesIdxDict[vertB]
        #adding the ridge into the list of ridges
        ridges_out.append(rdg)

        """
        print ('rdg %s' %rdg)
        #print ('crds orig %s' %vor.points[ pointA ,: ])
        print ('crds %s' %nodes_out[ int(rdg[0]) ,: ])
        print ('crds %s' %nodes_out[ int(rdg[1]) ,: ])
        print ('vr %s' %vor.vertices[vertA])
        print ('vr %s' %vor.vertices[vertB])
        print ('vor %s' %vertices_out[ int(rdg[3])])
        print ('vor %s' %vertices_out[ int(rdg[4])])
        #print ('vrt %s' %)
        #"""
    #
    aux_nodes = []
    v_count = len (vertices_out)
    vertIdxStart = node_count + len(aux_nodes)
    print ('vertIdxStart %d' %vertIdxStart)

    for i in range (len(ridges_out)):
        ridges_out[i][3] += vertIdxStart
        ridges_out[i][4] += vertIdxStart

    print('done.')
    sys.stdout.flush()
    #output: nodes_out, aux_nodes, vertices_out, ridges_out

    saveNodes(master_folder, nodes_out, "Particle",dim, nodesFile)
    #saveNodes(master_folder, nodes_out, aux_nodes, dim)
    #saveNodes(master_folder, aux_nodes, "AuxNode",dim, verticesFile)
#   saveVertices(master_folder, vertices_out, dim, withoutTransport=True)
    saveNodes(master_folder, vertices_out, "AuxNode",dim, verticesFile)
    saveMechanicalElements(master_folder, ridges_out, node_count, dim, nodes_out, mZ=mZ)
    #saveTransportElements(master_folder, ridges_out,dim, node_count, aux_nodes, maxLim)

    savePeriodicBlock(master_folder,cpldNds,maxLim, nodes_out)
    totalPointCount = len(nodes_out) + len(aux_nodes) + len(vertices_out)

    return v_count, verticesIdxDict, vertIdxStart, totalPointCount#, nodes_out, aux_nodes, vertices_out, ridges_out




def savePeriodicBlock (master_folder,cpldNds, maxLim, nodes_out):
    cf = open(os.path.join(master_folder,blocksFile),"w")

    ndepend = len(cpldNds)
    #ex ey gxy sx sy sxy
    cf.write("BasicPeriodicBC\tsize\t2\t%e\t%e\tload\t%d\tey\t%d\tgxy\t%d\tpairs\t%d"%(maxLim[0],maxLim[1],2,0,1, ndepend))


    for i in range(len(cpldNds)):
        cf.write("\t%d\t%d"%(cpldNds[i][1], cpldNds[i][0]))

        #plt.plot( [nodes_out[ cpldNds[i][0],0 ], nodes_out[ cpldNds[i][1],0 ]], [nodes_out[ cpldNds[i][0],1 ], nodes_out[ cpldNds[i][1],1 ]],'ro-', color='red')
        #plt.text(nodes_out[ cpldNds[i][0],0 ] , nodes_out[ cpldNds[i][0],1 ], cpldNds[i][0], fontsize=11)
        #plt.text(nodes_out[ cpldNds[i][1],0 ] , nodes_out[ cpldNds[i][1],1 ], cpldNds[i][1], fontsize=11)

    cf.write(os.linesep)
    cf.close()

    #plt.plot(nodes_out[:,0], nodes_out[:,1], 'o', color='black')

    #plt.show()












def output3D(master_folder, node_count, maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=None, notches=None):

    dim = 3
    print('Extracting the geometry...',  end ='')
    sys.stdout.flush()

    printout = False
    # nody: [x,y,z] [powerR] [area]
    nodes_out = np.zeros( (node_count, (dim + 1 +1)))
    nodes_out[:,  0:dim] = node_coords[:,  0:dim]

    if ((len(areas) == node_count)):
        nodes_out[:,dim] = areas[:]

    relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
    #print ('Area Error: %.5E ' %(relAreaError) )

    ########################################################################################################
    # ridges with nodes within sample
    validRidgeIdxs = []

    #print('ridge points')
    #adding ridges with at least one node in sample
    for i in range (vor.ridge_points.shape[0]):
        pr = False
        for p in range (2):
            if (vor.ridge_points[i][p] < node_count):
                pr=True
        if (pr):
            validRidgeIdxs.append(i)


    validRidgeIdxs = np.asarray(validRidgeIdxs)
    ########################################################################################################
    # vertices: [xA,yA,zA] [origIdx]
    vertices_out = []
    # dictionary of original and new indices of vertices
    verticesIdxDict = {}
    # ridges: nodeAidx, nodeBidx, trsprtBC, vertIdx
    ridges_out = []
    #auxiliary nodes
    aux_nodes = []
    ########################################################################################################
    allCoplanar = True
    for i in range (validRidgeIdxs.size):

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
            addVrtx = True
            for j in range (len(vertices_out)):
                if (vertices_out[j][dim] == vrtx[dim]):
                    addVrtx = False

            if (addVrtx == True):
                verticesIdxDict.update( { vrtx[dim] : len(vertices_out)  } )
                vrtx [dim +1] = len(vertices_out)
                vrtx [dim +2] = 0
                vertices_out.append(vrtx)

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
        newAuxNodes = saveTransportElements(master_folder, ridges_out,dim, node_count, v_count, aux_nodes, maxLim, nodes_out, vertices_out)
    vertIdxStart += newAuxNodes

    for i in range (len(ridges_out)):
        ln = len(np.asarray(ridges_out[i]) )
        for l in range (3, ln):
            ridges_out[i][l] += newAuxNodes


    if activeMechanics:
        saveNodes(master_folder, nodes_out, "Particle",dim, nodesFile)
        saveNodes(master_folder, aux_nodes, "AuxNode",dim, auxNodesFile)
        saveMechanicalElements(master_folder, ridges_out, node_count, dim, nodes_out, mZ=mZ, notches = notches)
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



def returnSelectedPts (boundPtA , boundPtB, points):
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

def excludeSelectedPts (boundPtA , boundPtB, points):
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
            bc[1:] = nodes_mechBCmerged[i].getMechBC()
        elif (dim == 3):
            bc = np.zeros (1 + 12)
            bc[0] = nodes_mechBCmerged[i].getNodeIdx()
            bc[1:] = nodes_mechBCmerged[i].getMechBC()

        mechBC_out.append(bc)

    #
    if (govNodesBC == False and len(mechBC_out)>0 ):
        print('Saving node BC...')
        if (dim == 2):
            headerLine = 'nodeIdx\tKinTrX\tKinTrY\tKinRotZ\tStTrX\tStTrY\tStRotZ'
            fl=open(os.path.join(master_folder,mechBCFile) ,'w')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='%d\t%d\t%d\t%d\t%d\t%d\t%d', header = headerLine)
            fl.close()
        elif (dim == 3):
            headerLine = 'nodeIdx\tKinTrX\tKinTrY\tKinTrZ\tKinRotX\tKinRotY\tKinRotZ\tStTrX\tStTrY\tStTrZ\tStRotX\tStRotY\tStRotZ'
            fl=open(os.path.join(master_folder,mechBCFile) ,'w')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d', header = headerLine)
            fl.close()
    elif(len(mechBC_out)>0 ):
        print('Saving rigid plate BC...')
        #print(mechBC_out)
        if (dim == 2):
            fl=open(os.path.join(master_folder,mechBCFile) ,'a')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='%d\t%d\t%d\t%d\t%d\t%d\t%d')
            fl.close()
        elif (dim == 3):
            fl=open(os.path.join(master_folder,mechBCFile) ,'a')
            np.savetxt(fl, mechBC_out, delimiter='\t', fmt='%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d')
            fl.close()

    print('done.')

def saveSolver(master_folder, solver, solStep, minStep, maxStep, simTime, limitTolerance= 1e-1, maxIt=20):
    f=open(os.path.join(master_folder,solverFile),'w')
    f.write('%s\n'%solver)
    f.write('time_step\t%e\n'%solStep)
    f.write('min_time_step\t%e\n'%minStep)
    f.write('max_time_step\t%e\n'%maxStep)
    f.write('total_time\t%f\n'%simTime)
    f.write('limit_tolerance\t%e\n'%limitTolerance)
    f.write('maxIt\t%d\n'%maxIt)

    f.close()

def saveMasterInput(master_folder,dim, solver, solStep, minStep, maxStep, simTime, activeTransport, activeMechanics, periodic=False, constraint=False, limitTolerance= 1e-1, maxIt=20):
     print('Saving master file...', end='')
     sys.stdout.flush()
     fl=open(os.path.join(master_folder,masterFile),'w')

     if not periodic:
         fl.write("Dimension\t%d\n"%dim)
         """
         if (solver == "SteadyStateLinearSolver"):
                fl.write('Solver\tSteadyStateLinearSolver\ttime_step\t%e\ttotal_time\t%e\n' %(solStep, simTime))
         if (solver == "SteadyStateNonLinearSolver"):
                fl.write('Solver\tSteadyStateNonLinearSolver\ttime_step\t%e\tmax_time_step\t%e\tmin_time_step\t%e\ttotal_time\t%e\tlimit_tolerance\t%e\tmaxIt\t%d\n' %(solStep,  maxStep, minStep, simTime, limitTolerance, maxIt))
         """
         fl.write("Solver\t%s\n"%(solverFile))
         saveSolver(master_folder, solver, solStep, minStep, maxStep, simTime, limitTolerance, maxIt)

         if not constraint:
             fl.write("NodeFiles\t3\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile))
         else:
             fl.write("NodeFiles\t4\t%s\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile, govNodesFile))
             fl.write('PBlockFiles\t1\t%s\n' %(constraintFile))
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
         fl.write("Dimension\t%d\n"%dim)
         if (solver == 0):
                fl.write('Solver\tSteadyStateNonLinearSolver\ttime_step\t%e\tmax_time_step\t%e\tmin_time_step\t%e\ttotal_time\t%f\n' %(solStep,  maxStep, minStep, simTime))

         fl.write("NodeFiles\t2\t%s\t%s\n"%(nodesFile,verticesFile))
         fl.write("MatFiles\t1\t%s\n"%materialsFile)
         fl.write("ElemFiles\t1\t%s\n"%(mechElemsFile))
         fl.write("PBlockFiles\t1\t%s\n" %blocksFile)
     fl.write("FunctionFiles\t1\t%s\n"%functionsFile)
     fl.write("ExporterFiles\t1\t%s"%exportersFile)
     #fl.write('INITMECH:\t%s\n' % initConditionsMechFile   )
     #fl.write('INITTRSPRT:\t%s\n' % initConditionsTrsprtFile   )

     fl.close()
     print('done.')

def saveMaterials (master_folder,materials):
    print ('Saving materials...', end='')
    sys.stdout.flush()
    ### MATERIALS
    with open(os.path.join(master_folder,materialsFile), 'w') as f:
        headerLine = 'matType\tYoungM\tPoisson\tTranspC\tTranspS\tDensity'
       # f.write("%s\n" % headerLine )
        for item in materials:
            f.write("%s\n" % item.getString() )
          # print (item.getString())
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



def saveTransportBC(master_folder,transportBCmerged, verticesDict, vertIdxStart):
    print('Saving TRSPRT boundary conditions...', end = '')
    sys.stdout.flush()
    trsptBC_out = []

    for i in range (len(transportBCmerged)):
        bc = np.zeros ((1 + 1 + 1))
        bc[0] = verticesDict[transportBCmerged[i].getNodeIdx()] + vertIdxStart
        bc[1:] = transportBCmerged[i].getTrsprtBC()
        trsptBC_out.append(bc)

    headerLine = 'vrtxIdx\tTrsptP\tTrsptJ'
    fl=open(os.path.join(master_folder,trsprtBCFile) ,'w')
    if (len(trsptBC_out)>0):
        np.savetxt(fl, trsptBC_out, delimiter='\t', fmt='%d\t%d\t%d', header = headerLine)
    fl.close()

    print('done.')


def saveExporters(master_folder,activeTransport, activeMechanics):
    print('Saving exporters...', end='')
    sys.stdout.flush()
    fl=open(os.path.join(master_folder,exportersFile),'w')
    if activeMechanics:
        fl.write('#TXTNodalExporter translations 2 ux uy\n')
        fl.write('#TXTNodalExporter pressure 1 pressure\n')
        fl.write('VTKElementExporter out  saveEvery 1e-4 cellData 1 damage\n')
        fl.write('#VTKRCExporter faces  saveEvery 1e-1 cellData 1 damage\n')
        fl.write('#TXTGaussPointExporter damageT 11 x y z normal_x normal_y normal_z damage strainTY strainTZ strainPLTY strainPLTZ\n')
    if activeTransport:
        fl.write('TXTNodalExporter pressure 1 pressure\n')

    fl.close()

    print('done.')



def saveNodes (master_folder,nodes_out, nodetype, dim, filename):
    print('Saving nodes: %s...' %nodetype, end='')
    sys.stdout.flush()
    nodes_out = np.array(nodes_out)
    #writing nodes
    #print(len(nodes_out))
    num = dim
    if (dim == 2):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY"
        fmt= nodetype + '\t%.15e\t%.15e'
    elif (dim == 3):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY\tnodeCrdZ"
        fmt= nodetype + '\t%.15e\t%.15e\t%.15e'
    if nodetype=="Particle":
        headerLine = headerLine + "\tpowRadius"
        fmt = fmt + '\t%.15e'
        num = num + 1

    fl=open(os.path.join(master_folder,filename),'w')
    np.savetxt(fl,  nodes_out[:,:num], delimiter='\t',   fmt=fmt,  header = headerLine)
    fl.close()
    print('done.')
    sys.stdout.flush()


def saveMechanicalElements (master_folder,ridges_out, node_count, dim, nodes, mZ=None, notches = None):
    print('Saving MECH elements...', end ='')
    sys.stdout.flush()
    #filtering ridges to ridges with both nodes in sample -> mech elements
    mechElemRidges = []
    for m in range (len(ridges_out)):
        if (ridges_out[m][0] < node_count and ridges_out[m][1] < node_count and ridges_out[m][0] >=0  and ridges_out[m][1] >= 0):
            mechElemRidges.append( ridges_out[m].copy() )
    print ('Mech elements: %d' %len(mechElemRidges))

    onlyMechNodesConnected = True

    if (mZ!=None and len(mZ)>0):
        print('Material zones recognized.')
        for i in range (len(mechElemRidges)):
            nodeA = nodes[int(mechElemRidges[i][0])]
            nodeB = nodes[int(mechElemRidges[i][1])]

            if (int(mechElemRidges[i][0]) >= node_count or int(mechElemRidges[i][1]) >= node_count):
                onlyMechNodesConnected = False

            if (dim==2):
                if ( (mZ[0][0][0] < nodeA[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeA[1] < mZ[0][1][1] and
                      mZ[0][0][0] < nodeB[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeB[1] < mZ[0][1][1]) or
                      (mZ[0][2][0] < nodeA[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeA[1] < mZ[0][3][1] and
                      mZ[0][2][0] < nodeB[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeB[1] < mZ[0][3][1])   ):
                    mechElemRidges[i] = np.hstack( (mechElemRidges[i], np.array([2])) )
                    #print('found ela element')
                else:
                    mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )

            if (dim==3):
                if ( mZ[0][0][0] < nodeA[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeA[1] < mZ[0][1][1] and
                      mZ[0][0][2] < nodeA[2] < mZ[0][1][2] and
                      mZ[0][0][0] < nodeB[0] < mZ[0][1][0] and
                      mZ[0][0][1] < nodeB[1] < mZ[0][1][1] and
                      mZ[0][0][2] < nodeB[2] < mZ[0][1][2] ):
                      mechElemRidges[i] = np.hstack( (mechElemRidges[i], np.array([2])) )
                      #print('node in A')
                if (mZ[0][2][0] < nodeA[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeA[1] < mZ[0][3][1] and
                      mZ[0][2][2] < nodeA[2] < mZ[0][3][2] and
                      mZ[0][2][0] < nodeB[0] < mZ[0][3][0] and
                      mZ[0][2][1] < nodeB[1] < mZ[0][3][1] and
                      mZ[0][2][2] < nodeB[2] < mZ[0][3][2]) :
                      mechElemRidges[i] = np.hstack( (mechElemRidges[i], np.array([2])) )
                      #print('node in B')
                    #print('found ela element')

                else:
                    mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )
    else:
        for i in range (len(mechElemRidges)):
            mechElemRidges[i] = np.hstack( (mechElemRidges[i],  np.array([0])) )

    if (onlyMechNodesConnected):
        print ('MechElems connect only MechNodes. That is ok.')
    else:
        print ('MechElems CONNECT WRONG NODES !!!')

    if (notches!=None ):
        print('Filtering out elements connecting notches...' )
        n = 0
        elementsWithoutNotch = []
        for i in range (len(mechElemRidges)):
            nA = mechElemRidges[i][0]
            nB = mechElemRidges[i][1]

            addElem = True
            n=0
            for notch in notches:
                if ((nA in notch[0] and nB in notch[1]) or ((nB in notch[0] and nA in notch[1]))):
                    #print('%d'%n)
                    addElem = False
                    break
                n+=1

            if addElem == True:
                elementsWithoutNotch.append(mechElemRidges[i])

        mechElemRidges = elementsWithoutNotch

        print('done.')


    if (dim == 2):
        headerLine = 'ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tvrtxAIdx\tvrtxBIdx\tMaterial'
        fl=open(os.path.join(master_folder,mechElemsFile),'w')
        np.savetxt(fl, mechElemRidges, delimiter='\t',fmt='LTCBEAM\t%d\t%d\t%d\t%d\t%d\t%d', header = headerLine )
        fl.close()

    if (dim == 3):
        headerLine = '#ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tverticesIdxs\tMaterial\n'
        fl=open(os.path.join(master_folder,mechElemsFile),'w')
        ro = np.asarray(mechElemRidges[0])
        fl.write(headerLine)
        for i in range (len(mechElemRidges)):
            ro = np.asarray(mechElemRidges[i])
            fmt='LTCBEAM\t%d\t%d\t%d'
            sh = ro.shape[0]
            np.savetxt(fl,  ro.reshape(1, sh) , delimiter='\t', fmt=fmt+'\t%d'*(sh-3)+ '\t0')

    sys.stdout.flush()


def saveTransportElements(master_folder,ridges_out, dim, node_count, vertCount, aux_nodes, maxLim, nodes_out, vertices_out):
    print('Creating TRSPRT elements...', end='')
    sys.stdout.flush()
    transportElements = []
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
                #print('%d ; %d = %d ; %d' %(n,m, ro[n], ro[m]))
                for elem in transportElements:
                    if ( ((elem.vertexA == ro[n]) and (elem.vertexB == ro[m]))
                      or ((elem.vertexA == ro[m]) and (elem.vertexB == ro[n])) ):
                    #if ( elem.vertexA in ro and elem.vertexB in ro ):
                        newPath = False
                        elem.addConnectedNodes(ro)
                        break
                if (newPath == True):
                    connNds = []
                    connNds.clear()
                    connNds.append (ro[0])
                    connNds.append (ro[1])
                    #if (ro[0]>node_count and ro[1]>node_count):
                    #    print ('both aux')
                    transportElements.append (utilitiesMech.transportPath (ro[n], ro[m], connNds.copy(), 1))
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
    updatedElems =0
    newAuxNodesA = []
    if (auxNodesInitLength != 0):
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
                if np.abs(vertexA[d]-vertexB[d])<1e-8:
                    #print(np.abs(vertexA[d]-vertexB[d]))
                    equalCoords +=1
                else:
                    diffIdx = d


            if (elem.connectedNodes[0]>=len(nodes_out) and elem.connectedNodes[len(elem.connectedNodes)-1]>=len(nodes_out)
            and (elem.connectedNodes[0]/node_count)
            != (elem.connectedNodes[len(elem.connectedNodes)-1]/node_count) and equalCoords ==2):

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
                #adding new aux node to connected nodes
                #print('old elem: %s' %elem.connectedNodes)
                elem.addSingleConnectedNode( elem.connectedNodes[len(elem.connectedNodes)-1])
                elem.addSingleConnectedNode( node_count + len(aux_nodes) )
                elem.addSingleConnectedNode( node_count + len(aux_nodes) )
                elem.addSingleConnectedNode( elem.connectedNodes[0] )

                aux_nodes.append(nanode)
                newAuxNodesA.append(nanode)

                """
                #print(ridgeCoords)
                ridgeCoords = []
                lasti = -1
                for i in ((elem.connectedNodes)):
                    if (i!=lasti):
                        print(i)
                        if (i<len(nodes_out)):
                            print(nodes_out[int(i),0:3])
                            ridgeCoords.append(nodes_out[int(i),0:3])
                        if (i>=len(nodes_out)):
                            print(aux_nodes[int(i-len(nodes_out))][0:3])
                            ridgeCoords.append(aux_nodes[int(i-len(nodes_out))][0:3])
                    lasti=i
                ridgeCoords = np.asarray(ridgeCoords)

                fig = plt.figure()
                ax = Axes3D(fig)
                ax.set_aspect('equal')
                ax.plot3D([vertexA[0], vertexB[0]], [vertexA[1], vertexB[1]], [vertexA[2], vertexB[2]], marker='o')
                #ax.plot3D([anodeA[0], anodeB[0]], [anodeA[1], anodeB[1]], [anodeA[2], anodeB[2]], marker='o')
                #ax.scatter3D(elemMidpoint[0], elemMidpoint[1], elemMidpoint[2])
                #ax.scatter3D(anodeA[0], anodeA[1], anodeA[2])
                #ax.scatter3D(anodeB[0], anodeB[1], anodeB[2])
                #ax.scatter3D(ridgeMidpoint[0], ridgeMidpoint[1], ridgeMidpoint[2])
                for r in range (len(ridgeCoords)-1):
                    ax.plot3D([ridgeCoords[r,0], ridgeCoords[r+1,0]], [ridgeCoords[r,1], ridgeCoords[r+1,1]], [ridgeCoords[r,2], ridgeCoords[r+1,2]], marker='x')
                ax.plot3D([ridgeCoords[0,0], ridgeCoords[len(ridgeCoords)-1,0]], [ridgeCoords[0,1], ridgeCoords[len(ridgeCoords)-1,1]], [ridgeCoords[0,2], ridgeCoords[len(ridgeCoords)-1,2]], marker='x')

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
                plt.show()
                """
                #print('new elem: %s' %elem.getString())
                #print('new elem: %s' %elem.getStringyString(len(nodes_out), auxNodesInitLength))
                #print(elem.connectedNodes)
                #print()
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
    print('Saving TRSPRT elements...', end='')
    sys.stdout.flush()
    print('Trsprt elements: %d' %len(transportElements))
    with open(os.path.join(master_folder,trsprtElemsFile), 'w') as f:
        headerLine = '#ElemType\tvrtxAIdx\tvrtxBIdx\tnrOfNodes\tnodesIdx\tMaterial'
        f.write("%s\n" % headerLine )
        for element in transportElements:
            #print ("%s\n" % element.getString() )
            if (dim==2):f.write("%s\n" % element.getString() )
            if (dim==3): f.write("%s\n" % element.getReducedString() )


    return newAuxNodes

def saveRigidPlates(master_folder, dim, rigidPlates, totalNodeCount):
    print('Saving rigid plates...', end='')
    with open(os.path.join(master_folder,constraintFile), 'w') as f:
        headerLine = '#ConstraintType\tGovNodeIdx\tXmin\tXmax\tYmin\tYmax'
        if (dim==3):
            headerLine += '\tZmin\tZmax'
        f.write("%s\n" % headerLine )
        for i in range(len(rigidPlates)):
            rplt = rigidPlates[i]
            rplt.govNodeIdx = totalNodeCount + i
            f.write("%s\n" % rplt.getString() )

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

def saveConstraint(master_folder, dim, govNodes, govNodesMechBC, rigidPlates, totalNodeCount, nodes):
    #saving gov nodes
    saveNodes (master_folder,govNodes, "MasterNode", dim, govNodesFile)
    #saving gov nodes mech BC
    for i in range (len(govNodesMechBC)):
        m = govNodesMechBC[i]
        m.nodeIdx = totalNodeCount + i
    saveMechBC(master_folder,dim, govNodesMechBC, govNodesBC = True)

    #saving force gauges for rigid plates
    for i in range (len(govNodesMechBC)):
        saveForceGauges(master_folder, dim, govNodesMechBC[i].nodeIdx)

    #saving rigid plates
    saveRigidPlates(master_folder, dim, rigidPlates, totalNodeCount)

    for rp in range (len(rigidPlates)):
        print('Nodes affected by Rigid plate #%d:' %rp)
        print(rigidPlates[rp].getNodesAffected(nodes))



def saveForceGauges(master_folder, dimension, nodeIdx, moments=True):
    saveForceGauge(master_folder, 'fx#%d'%nodeIdx , 'fx', nodeIdx )
    saveForceGauge(master_folder, 'fy#%d'%nodeIdx, 'fy', nodeIdx )
    if (dimension==3): saveForceGauge(master_folder, 'fz#%d'%nodeIdx, 'fz', nodeIdx )
    if moments == True:
        if (dimension==3): saveForceGauge(master_folder, 'mx#%d'%nodeIdx , 'mx', nodeIdx )
        if (dimension==3): saveForceGauge(master_folder, 'my#%d'%nodeIdx, 'my', nodeIdx )
        saveForceGauge(master_folder, 'mz#%d'%nodeIdx, 'mz', nodeIdx )

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
    for mg in measuringGauges:
        saveDisplacementGauges(master_folder, dimension, mg.name, mg.coordsA, mg.coordsB, rotations = mg.rotation)


def checkSavedModel(master_folder, dim, activeMechanics, activeTransport):
    print('\nChecking  generated files...')
    allOK = True
    cols = np.arange(1,dim+1)

    print('Loading back node coords...', end='')
    test_nodeCoords = np.genfromtxt(os.path.join(master_folder,nodesFile),  dtype= None, encoding='ascii', usecols=cols)
    print('\t\t %d nodes loaded.' %len(test_nodeCoords))

    print('Loading back aux node coords...', end='')
    test_auxNodeCoords = np.genfromtxt(os.path.join(master_folder,auxNodesFile),  dtype= None, encoding='ascii', usecols=cols)
    print('\t\t %d nodes loaded.' %len(test_auxNodeCoords))

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
                #print(mechElem[3:3+verticesNr])
                #print(vertices)
        if (wrongElems==0):
            print('All faces coplanar. Mech Elems OK. MaxErr: %e' %maxErr)
        else:
            print ('Wrong faces: %d !!!!' %wrongElems)
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
            allCoplanar, val = checkCoplanarity(vertices, 1e-15)
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











#
