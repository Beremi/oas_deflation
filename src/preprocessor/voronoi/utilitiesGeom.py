import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os

import utilitiesMech

master_folder = "coupled_problem"
try:
    if not os.path.exists(master_folder):
        os.makedirs(master_folder)
except:
    print('Please create directory %s! Code Exited.'%master_folder)
    sys.exit()

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


#
#coplanarity test
def equation_plane(pA, pB, pC, pD):
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


def extractGeometry (dim, node_count, maxLim, vor, node_coords, areas):
    if (dim == 2):
        vert_count, verticesIdxDict, vertIdxStart = output2D(node_count,  maxLim, vor, node_coords, areas)
    if (dim == 3):
        vert_count, verticesIdxDict, vertIdxStart = output3D(node_count,  maxLim, vor, node_coords, areas)

    return vert_count, verticesIdxDict, vertIdxStart


#Extract geometry 2d
def output2D(node_count,  maxLim, vor, node_coords, areas):
    dim = 2
    print('Extracting the geometry...', end='')
    sys.stdout.flush()
    nodes_out = np.zeros( (node_count, (2 + 1 + 1 )))
    #nodes_out[:,  0:2] = node_coords[:,  0:2]
    nodes_out[:,  0:2] = vor.points[0:node_count , 0:2]
    nodes_out[:,dim] = 0
    nodes_out[:,dim + 1] = areas[:]

    #relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
    #print ('Area error: %.5e ' %(relAreaError)  )
    ########################################################################################
    validRidgeIdxs = []
    for i in range (vor.ridge_points.shape[0]):
        pr = False
        for p in range (2):
            if (vor.ridge_points[i][p] < node_count):
                pr=True

        if (pr):
           validRidgeIdxs.append(i)

    validRidgeIdxs = np.asarray(validRidgeIdxs)

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

    saveNodes(nodes_out, aux_nodes, dim)
    saveVertices(vertices_out, dim)
    saveMechanicalElements(ridges_out, node_count, dim)
    saveTransportElements(ridges_out,dim, node_count, aux_nodes, maxLim)

    return v_count, verticesIdxDict, vertIdxStart#, nodes_out, aux_nodes, vertices_out, ridges_out




def output3D(node_count, maxLim, vor, node_coords, areas):
    dim = 3
    print('Extracting the geometry...',  end ='')
    sys.stdout.flush()

    printout = False
    # nody: [x,y,z] [powerR] [area]
    nodes_out = np.zeros( (node_count, (dim + 1 +1)))
    nodes_out[:,  0:dim] = node_coords[:,  0:dim]
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

    """
    #it is a question if consider ridges with both nodes out of sample
    #but both are part of at least one ridge with one node in sample
    """

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
        """
        for v in range ( nrVertices-3 ):
            pA = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v]][:]
            pB = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+1]][:]
            pC = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+2]][:]
            pD = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+3]][:]

            tol = 1e-10
            val = equation_plane(pA, pB, pC, pD)
            if ( val > tol):
                allCoplanar = False
                print('Not coplanar!!! Ridge nr. %d, err: %e' %(i, val ))
            #else: print('Coplanar  %d' %i)
        """
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
        if (printout):print ('ALL ridges coplanar OK')
    else:
        print ('!!! NOT ALL RIDGES COPLANAR !!!')
    print('done.')
    vertIdxStart = node_count + len(aux_nodes)
    v_count = len (vertices_out)

    for i in range (len(ridges_out)):
        ln = len(np.asarray(ridges_out[i]) )
        for l in range (3, ln):
            ridges_out[i][l] += vertIdxStart



    newAuxNodes = saveTransportElements(ridges_out,dim, node_count, aux_nodes, maxLim)

    for i in range (len(ridges_out)):
        ln = len(np.asarray(ridges_out[i]) )
        for l in range (3, ln):
            ridges_out[i][l] += newAuxNodes

    saveNodes(nodes_out, aux_nodes,dim)
    saveVertices(vertices_out, dim)
    saveMechanicalElements(ridges_out, node_count, dim)


    return v_count, verticesIdxDict, vertIdxStart



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

    return selectedPointIdxs



def saveTransportIC(transportIC_merged):
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



def saveMechIC(dim, nodes_mechICmerged):
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




def saveMechBC(dim, nodes_mechBCmerged):
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

    print('done.')

def saveMasterInput(dim, solver, solStep, simTime):
     print('Saving master file...', end='')
     sys.stdout.flush()
     fl=open(os.path.join(master_folder,masterFile),'w')

     fl.write("Dimension\t%d\n"%dim)
     if (solver == 0):
            fl.write('Solver\tSteadyStateNonLinearSolver\ttime_step\t%e\ttotal_time\t%f\n' %(solStep, simTime))
     fl.write("NodeFiles\t3\t%s\t%s\t%s\n"%(nodesFile,auxNodesFile,verticesFile))
     fl.write("MatFiles\t1\t%s\n"%materialsFile)
     fl.write("ElemFiles\t2\t%s\t%s\n"%(mechElemsFile,trsprtElemsFile))
     fl.write("BCFiles\t2\t%s\t%s\n"%(mechBCFile,trsprtBCFile))
     fl.write("ICFiles\t2\t%s\t%s\n"%(mechICFile,trsprtICFile))
     fl.write("FunctionFiles\t1\t%s\n"%functionsFile)
     fl.write("ExporterFiles\t1\t%s"%exportersFile)

     #fl.write('INITMECH:\t%s\n' % initConditionsMechFile   )
     #fl.write('INITTRSPRT:\t%s\n' % initConditionsTrsprtFile   )

     fl.close()
     print('done.')

def saveMaterials (materials):
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

def saveFunctions (functions):
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



def saveTransportBC(transportBCmerged, verticesDict, vertIdxStart):
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
    np.savetxt(fl, trsptBC_out, delimiter='\t', fmt='%d\t%d\t%d', header = headerLine)
    fl.close()

    print('done.')


def saveExporters():
    print('Saving exporters...', end='')
    sys.stdout.flush()
    fl=open(os.path.join(master_folder,exportersFile),'w')
    fl.write("TXTNodalExporter translations 2 ux uy\n")
    fl.write("TXTNodalExporter pressure 1 pressure\n")
    fl.write('VTKElementExporter out  time_every 1e-20 0')
    fl.close()

    print('done.')



def saveNodes (nodes_out, aux_nodes, dim):
    print('Saving nodes...', end='')
    sys.stdout.flush()
    #writing nodes
    if (dim == 2):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY\tpowRadius"
        fmt='Particle\t%.12f\t%.12f\t%.12f'
    if (dim == 3):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY\tnodeCrdZ\tpowRadius"
        fmt='Particle\t%.12f\t%.12f\t%.12f\t%.12f'

    fl=open(os.path.join(master_folder,nodesFile),'w')
    np.savetxt(fl,  nodes_out[:,  0:dim+1], delimiter='\t',   fmt=fmt,  header = headerLine)
    fl.close()

    #writing aux nodes
    if (dim == 2):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY"
        fmt='AuxNode\t%.12f\t%.12f'
    if (dim == 3):
        headerLine  = "Type\tnodeCrdX\tnodeCrdY\tnodeCrdZ"
        fmt='AuxNode\t%.12f\t%.12f\t%.12f'

    fl=open(os.path.join(master_folder,auxNodesFile),'w')
    np.savetxt(fl,  aux_nodes, delimiter='\t',   fmt=fmt,  header = headerLine)
    fl.close()
    print('done.')
    sys.stdout.flush()



def saveVertices (vertices_out, dim):
    print('Saving vertices...', end='')
    sys.stdout.flush()
    if (dim == 2):
        headerLine = 'Type\tvrtxCrdX\tvrtxCrdY'
        fmt ='TrsprtNode\t%.12f\t%.12f'
    if (dim == 3):
        headerLine = 'Type\tvrtxCrdX\tvrtxCrdY\tvrtxCrdZ'
        fmt = 'TrsprtNode\t%.12f\t%.12f\t%.12f'

    vertices_print = np.asarray(vertices_out)
    fl=open(os.path.join(master_folder,verticesFile),'w')
    np.savetxt(fl, vertices_print[:, 0:dim], delimiter='\t', fmt = fmt, header = headerLine)
    fl.close()
    print('done.')
    sys.stdout.flush()


def saveMechanicalElements (ridges_out, node_count, dim):
    print('Saving MECH elements...', end ='')
    sys.stdout.flush()
    #filtering ridges to ridges with both nodes in sample -> mech elements
    mechElemRidges = []
    for m in range (len(ridges_out)):
        if (ridges_out[m][0] < node_count and ridges_out[m][1] < node_count and ridges_out[m][0] >=0  and ridges_out[m][1] >= 0):
            mechElemRidges.append( ridges_out[m] )

    if (dim ==2):
        headerLine = 'ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tvrtxAIdx\tvrtxBIdx\tMaterial'
        fl=open(os.path.join(master_folder,mechElemsFile),'w')
        np.savetxt(fl, mechElemRidges, delimiter='\t',fmt='LTCBEAM\t%d\t%d\t%d\t%d\t%d\t0', header = headerLine )
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

    print('done.')
    sys.stdout.flush()


def saveTransportElements(ridges_out, dim, node_count, aux_nodes, maxLim):
    print('Creating TRSPRT elements...', end='')
    sys.stdout.flush()
    transportElements = []
    ridges_out = np.asarray(ridges_out)
    if (dim == 2):
        for i in range (len(ridges_out)):
            connNds = []
            connNds.append (ridges_out[i,0])
            connNds.append (ridges_out[i,1])
            vrtA = ridges_out[i,3]
            vrtB = ridges_out[i,4]
            trp = utilitiesMech.transportPath (vrtA, vrtB, connNds, 1)
            transportElements.append (trp)
    if (dim==3):
        for i in range (len(ridges_out)):
            ro = np.asarray(ridges_out[i])
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
    print('done.')
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


    print('Generating additional aux_nodes...', end='')
    auxNodesInitLength = len (aux_nodes)
    beamMidpoint = maxLim /2
    for elem in transportElements:
        if (elem.connectedNodes[0]>=node_count and elem.connectedNodes[len(elem.connectedNodes)-1]>=node_count):
            #print('corner line: %s' %elem.getReducedString())
            anodeA = np.asarray (aux_nodes[ int(elem.connectedNodes[0]-node_count) ][:])
            anodeB = np.asarray (aux_nodes[ int(elem.connectedNodes[len(elem.connectedNodes)-1]-node_count) ][:])
            #print( anodeA )
            #print( anodeB )
            nanode = np.zeros(3)
            for i in range (dim):
                if ( scipy.spatial.distance.euclidean(anodeB[i], beamMidpoint) >
                     scipy.spatial.distance.euclidean(anodeA[i], beamMidpoint) ):
                    nanode[i] = anodeB[i]
                else:
                    nanode[i] = anodeA[i]
            #print('New aux node: %s' %nanode)
            aux_nodes.append(nanode)
            #
            #adding new aux node to connected nodes
            #print('old elem: %s' %elem.connectedNodes)
            elem.connectedNodes.append( node_count + len(aux_nodes) )
            elem.connectedNodes.append( elem.connectedNodes[0] )
            #print('new elem: %s' %elem.connectedNodes)
            #print()
    print('done.')
    sys.stdout.flush()

    print('Another renumbering of vertices...', end='')
    auxNodesFinalLength = len (aux_nodes)
    newAuxNodes = auxNodesFinalLength - auxNodesInitLength
    for elem in transportElements:
        elem.vertexA += newAuxNodes
        elem.vertexB += newAuxNodes
    print('done.')
    sys.stdout.flush()

    print('Saving TRSPRT elements...', end='')
    sys.stdout.flush()
    with open(os.path.join(master_folder,trsprtElemsFile), 'w') as f:
        headerLine = '#ElemType\tvrtxAIdx\tvrtxBIdx\tnrOfNodes\tnodesIdx\tMaterial'
        f.write("%s\n" % headerLine )
        for element in transportElements:
            #f.write("%s\n" % element.getString() )
            f.write("%s\n" % element.getReducedString() )
    print('done.')

    return newAuxNodes
