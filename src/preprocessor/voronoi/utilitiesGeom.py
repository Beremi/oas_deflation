import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os
from IPython.display import clear_output

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

#check if any number in matrix is lower than
def checkLowerThan (matrix, minDist):
    #return not (all(x < minDist for x in matrix))
    return   all(i >= minDist for i in matrix)


# generates random points no closer to each other than minDist
# into 2d or 3d block
# maxLim: n-d array of dimensions
def generateNodesRect(maxLim, minDist, dim, trials, node_coords):
    if (dim==2):
        print('Generating 2d block segment of size: %f / %f. This may take few minutes. Do not panic. ' %(maxLim[0], maxLim[1]) )

    if (dim==3):
        print('Generating 3d block segment of size: %f / %f / %f. This may take long. Keep calm.' %(maxLim[0], maxLim[1], maxLim[2]) )

    generatedPoints = 0
    tr = 0
    while (tr<trials):
        tr = 0
        #
        distIsGood = False
        while (distIsGood == False):
            coords = np.random.random(dim)
            coords *= maxLim
            #
            distIsGood = True
            #
            #############################################x
            ######old sequential computation of distances
            """
            for p in range (len(node_coords)):
                #if (i!=p):
                    distInt = scipy.spatial.distance.euclidean(node_coords[p], coords)
                    #
                    if (distInt < minDist):
                        distIsGood = False
                        tr += 1
            """
            #############################################x
            ######new using c dist. The best so far.

            ncrds = np.asarray(node_coords)
            #print (ncrds)
            crds = np.asarray(coords)
            crds = np.reshape(crds, (-1, dim))
            #print (crds)
            dists = scipy.spatial.distance.cdist(crds, ncrds , 'euclidean')
            dists = dists.flatten()
            #print(dists)
            distIsGood = checkLowerThan(dists, minDist)
            if (distIsGood == False):
                tr += 1

            #############################################x
            #trying scipy cKDTree for searching for nearest neighbors. About the same performance as cDist so far.
            #does not work in parallel asi it should with n_jobs = -1
            ##############################################
            """
            #node_coords.append(coords)
            crds = np.asarray(coords)
            pts = np.asarray (node_coords)
            tree = scipy.spatial.cKDTree ( pts,  leafsize=5 )
            #
            violatingPoints = tree.query_ball_point (x = crds,  r = minDist, n_jobs = -1 )
            #print (violatingPoints)
            if ( len(violatingPoints) != 0):
                distIsGood = False
                tr += 1
            #else:
            #    print('GOOD POINT')
            #print(trials)
            """

            if (tr > trials): break
        if (tr > trials): break
        #
        #Adding node coords:
        if (tr < trials):
            node_coords.append(coords)
            generatedPoints  += 1
        #print(generatedPoints)
        #


#generates random points onto a set 3d line. No closer than minDst
#catch corners: samples the boundary points first
def generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners):
    print('Generating 3d line segment from [%f; %f; %f] to [%f; %f; %f] '
     %(nodeA[0], nodeA[1],nodeA[2],nodeB[0], nodeB[1],nodeB[2]) )
    generatedPoints = 0

    if(catchCorners):
        node_coords.append(np.copy(nodeA))
        node_coords.append(np.copy(nodeB))
        generatedPoints  += 2

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = np.zeros(dim)
            #
            r = np.random.uniform()
            coords[0] = (nodeB[0] - nodeA[0])*r  + nodeA[0]
            coords[1] = (nodeB[1] - nodeA[1])*r  + nodeA[1]
            coords[2] = (nodeB[2] - nodeA[2])*r  + nodeA[2]
            #
            distIsGood = True
            #
            for p in range (len(node_coords)):
                distInt = scipy.spatial.distance.euclidean(node_coords[p], coords)
                #
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    #clear_output(True)
                    #print('Line Pt: %d Tr: %d' %(generatedPoints, tr))
                    break
                if (tr > trials): break
            if (tr > trials): break
        #
        #Adding node coords
        if (tr < trials):
            node_coords.append(coords)
       # node_coords [i,:] = coords
        generatedPoints  += 1
        #



def generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials):
    print('Generating 3d surface segment from [%f; %f; %f] to [%f; %f; %f]'
     %(nodeA[0], nodeA[1],nodeA[2],nodeB[0], nodeB[1],nodeB[2]) )
    generatedPoints = 0

    tr=0
    while (tr<trials):
        tr = 0;
        #
        distIsGood = False
        while (distIsGood == False):
            coords = np.zeros(dim)
            #
            for c in range (dim):
                if (nodeA[c] == nodeB[c]):
                    coords[c] = nodeA[c]
                else:
                    coords[c] = (nodeB[c] - nodeA[c])*np.random.uniform()  + nodeA[c]

            distIsGood = True
            #
            for p in range (len(node_coords)):
                distInt = scipy.spatial.distance.euclidean(node_coords[p], coords)
                #
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    #clear_output(True)
                    #print('Line Pt: %d Tr: %d' %(generatedPoints, tr))
                    break
                if (tr > trials): break
            if (tr > trials): break
        #
        #Adding node coords
        #
        if (tr < trials):
            node_coords.append(coords)
       # node_coords [i,:] = coords
        generatedPoints  += 1


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
    generatedPoints = 0

    if(catchCorners):
        node_coords.append(np.copy(nodeA))
        node_coords.append(np.copy(nodeB))
        generatedPoints  += 2

    if (not equidist):
        tr=0
        while (tr<trials):
            tr = 0;
            #
            distIsGood = False
            while (distIsGood == False):

                coords = np.zeros(dim)
                r = np.random.uniform()
                coords[0] = (nodeB[0] - nodeA[0])*r  + nodeA[0]
                coords[1] = (nodeB[1] - nodeA[1])*r  + nodeA[1]

                distIsGood = True
                #
                for p in range (len(node_coords)):
                    #if (i!=p):
                    distInt = scipy.spatial.distance.euclidean(node_coords[p], coords)
                    #
                    if (distInt < minDist):
                        distIsGood = False
                        tr += 1
                        break

                    if (tr > trials): break
                if (tr > trials): break
            #
            #Adding node coords
            #
            if (tr < trials):
                node_coords.append(coords)
                generatedPoints  += 1
    else:
        #print('Equid')
        mD = minDist * 1.6
        length = np.linalg.norm(nodeA - nodeB)
        nodeNr = int (length / mD)
        indnt = (length-(nodeNr-1)*mD) / 2
        for i in range (nodeNr):
            #print(i)
            coords = np.zeros(2)
            coords[0] = (nodeB[0] - nodeA[0])*indnt/length +(nodeB[0] - nodeA[0])*mD/length*i  + nodeA[0]
            coords[1] = (nodeB[1] - nodeA[1])*indnt/length +(nodeB[1] - nodeA[1])*mD/length*i  + nodeA[1]
            node_coords.append(coords)
            #print (len(node_coords))
            generatedPoints  += 1


#OUTPUT METHODS
def output2D(node_count, dim, maxLim, vor, node_coords, areas, reOrderedIdxs,  diagonalize):
    print('Extracting the geometry...', end='')
    sys.stdout.flush()
    nodes_out = np.zeros( (node_count, (2 + 1 + 1 + 1 +1)))
    nodes_out[:,  0:2] = node_coords[:,  0:2]
    nodes_out[:,dim] = 0
    nodes_out[:,dim + 1] = areas[:]

    relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
    #print (np.sum(areas))
    #print ('Area Error: %.5E ' %(relAreaError) )

    for n in range (len(node_coords) ):
        nodes_out[n, dim+2] = reOrderedIdxs[n]

        mechBCidx = -1
        nodes_out[n, dim+3]  =  mechBCidx

    ########################################################################################
    validRidgeIdxs = []

    #print('ridge points')
    for i in range (vor.ridge_points.shape[0]):
        pr = False
        for p in range (2):
            if (vor.ridge_points[i][p] < node_count):
                pr=True

        if (pr):
           # print(vor.ridge_points[i,:])
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
        #Array for ridge nAidx, nBidx, trBc, vertAidx, vertBidx
        rdg = np.zeros ( (2 + 1 +  2) )

        #indices of two nodes that are divided by the ridge
        pointA = vor.ridge_points[validRidgeIdxs[i],0]
        pointB = vor.ridge_points[validRidgeIdxs[i],1]

        #creating auxiliary nodes if one of nodes is outside
        if(pointA >= node_count and pointB<node_count):
            ptA = np.zeros((2))
            ptA[0] = (vor.points[pointB, 0] + vor.points[pointA, 0]  ) /2
            ptA[1] = (vor.points[pointB, 1] + vor.points[pointA, 1]  ) /2

            pointA = node_count + len(aux_nodes)
            aux_nodes.append(ptA)

        if(pointB >= node_count  and pointA<node_count):
            ptB = np.zeros((2))
            ptB[0] = (vor.points[pointB, 0] + vor.points[pointA, 0]  ) /2
            ptB[1] = (vor.points[pointB, 1] + vor.points[pointA, 1]  ) /2

            pointB = node_count + len(aux_nodes)
            aux_nodes.append(ptB)
        #
        rdg[0] = pointA
        rdg[1] = pointB

        #number of vertices
        rdg[2] = 2

        #indices of vertices
        rdg[3] = verticesIdxDict[vertA] #vrtxA [dim] #verticesIdxDict[vertA]
        rdg[4] = verticesIdxDict[vertB] #vrtxB [dim] #verticesIdxDict[vertA]

        #adding the ridge into the list of ridges
        ridges_out.append(rdg)

        #print ('Ridge nr %d: vertices %d and %d' %(validRidgeIdxs[i], vertA, vertB ) )
        #print ('btw pts: %d and %d' %( pointA, pointB ) )

    print('done.')
    sys.stdout.flush()

    print('Saving nodes...', end='')
    sys.stdout.flush()
    #writing nodes
    ##############################################
    headerLine  = "Type\tnodeCrdX\tnodeCrdY\tpowRadius"
    fmt='Particle\t%.12f\t%.12f\t%.12f'

    if (diagonalize):
        nodes_backup = np.copy(nodes_out)
        #
        for x in range (len(nodes_out)):
            nodes_out[x, :] = nodes_backup[order[x],:]

    fl=open(os.path.join(master_folder,nodesFile),'w')
    np.savetxt(fl,  nodes_out[:,  0:3], delimiter='\t',   fmt=fmt,  header = headerLine)
    fl.close()

    #writing aux nodes
    headerLine  = "Type\tnodeCrdX\tnodeCrdY"
    fmt='AuxNode\t%.12f\t%.12f'
    fl=open(os.path.join(master_folder,auxNodesFile),'w')
    np.savetxt(fl,  aux_nodes, delimiter='\t',   fmt=fmt,  header = headerLine)
    fl.close()
    print('done.')
    sys.stdout.flush()

    print('Saving vertices...', end='')
    sys.stdout.flush()
    #writing vertices
    ###############################################
    headerLine = 'Type\tvrtxCrdX\tvrtxCrdY'

    vertices_print = np.asarray(vertices_out)

    v_count = len (vertices_out)

    fl=open(os.path.join(master_folder,verticesFile),'w')
    np.savetxt(fl, vertices_print[:, 0:2], delimiter='\t', fmt='TrsprtNode\t%.12f\t%.12f', header = headerLine)
    fl.close()

    vertIdxStart = node_count + len(aux_nodes)

    for i in range (len(ridges_out)):
        ridges_out[i][3] += vertIdxStart
        ridges_out[i][4] += vertIdxStart
        #
        nA = ridges_out[i][0]
        nB = ridges_out[i][1]
        vA = ridges_out[i][3]
        vB = ridges_out[i][4]
        #
        ridges_out[i][0] = vA
        ridges_out[i][1] = vB
        ridges_out[i][2] = 2
        ridges_out[i][3] = nA
        ridges_out[i][4] = nB
        #
    print('done.')
    sys.stdout.flush()

    print('Saving TRSPRT elements...', end='')
    sys.stdout.flush()
    #writing ridges - transport elements
    ############################################### ridges: idx noduA, idx noduB, transport okr. podminka, idx vertexu
    headerLine = 'ElemType\tvrtxAIdx\tvrtxBIdx\tnrOfNodes\tnodeAidx\tnodeBidx\tMaterial'
    fl=open(os.path.join(master_folder,trsprtElemsFile),'w')
    np.savetxt(fl, ridges_out, delimiter='\t',fmt='LTCTRSP\t%d\t%d\t%d\t%d\t%d\t1',
          header = headerLine )
    fl.close()

    for i in range (len(ridges_out)):
        #
        nA = ridges_out[i][0]
        nB = ridges_out[i][1]
        vA = ridges_out[i][3]
        vB = ridges_out[i][4]
        #
        ridges_out[i][0] = vA
        ridges_out[i][1] = vB
        ridges_out[i][2] = 2
        ridges_out[i][3] = nA
        ridges_out[i][4] = nB
    print('done.')
    sys.stdout.flush()

    print('Saving MECH elements...', end ='')
    sys.stdout.flush()
    #filtering ridges to ridges with both nodes in sample
    mechElemRidges = []
    for m in range (len(ridges_out)):
        if (ridges_out[m][0] < node_count and ridges_out[m][1] < node_count and ridges_out[m][0] >=0  and ridges_out[m][1] >= 0):
            mechElemRidges.append( ridges_out[m] )

  #  if (diagonalize):
     #   for x in range (len(mechElemRidges)):
     #       tady se to prehazi
    headerLine = 'ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tvrtxAIdx\tvrtxBIdx\tMaterial'
    fl=open(os.path.join(master_folder,mechElemsFile),'w')
    np.savetxt(fl, mechElemRidges, delimiter='\t',fmt='LTCBEAM\t%d\t%d\t%d\t%d\t%d\t0', header = headerLine )
    fl.close()
    print('done.')
    sys.stdout.flush()

    return v_count, verticesIdxDict, vertIdxStart





def output3D(node_count, dim, maxLim, vor, node_coords, areas, reOrderedIdxs,  mechBC_merged,  materials):
    ############################################################################################
    ############################################################################################
    ###################################### SAVING LATTICE MODEL GEOMETRY #######################
    ############################################################################################
    print('Extracting the geometry...',  end ='')
    sys.stdout.flush()

    printout = False
     # nody: [x,y,z] [powerR] [area]
    nodes_out = np.zeros( (node_count, (dim + 1 + 1 +1+1)))

    for d in range (dim):
        nodes_out[:,d] = node_coords[:,d]

    nodes_out[:,  0:dim] = node_coords[:,  0:dim]

    nodes_out[:,dim] = -1
    nodes_out[:,dim + 1] = areas[:]

    relAreaError = (np.sum(areas) - np.product(maxLim)) / np.product(maxLim)
   # print (np.sum(areas))
    if (printout): print ('Area Error: %.5E ' %(relAreaError) )

    ########################################################################################################
    # ridges, ktere maji nody ve vzorku
    validRidgeIdxs = []

    #print('ridge points')
    for i in range (vor.ridge_points.shape[0]):
        pr = False
        for p in range (2):
            if (vor.ridge_points[i][p] < node_count):
                pr=True

        if (pr):
            #print(vor.ridge_points[i,:])
            validRidgeIdxs.append(i)

    validRidgeIdxs = np.asarray(validRidgeIdxs)


    ########################################################################################################
    # vertices: [xA,yA,zA] [origIdx]
    vertices_out = []

    # slovnik originalnich a novych indexu vertexu
    verticesIdxDict = {}

    # ridges: idx noduA, idx noduB, transport okr. podminka, idx vertexu
    ridges_out = []

    #auxiliary nodes
    aux_nodes = []

    #vertices
    ########################################################################################################
    #list of valid vertices
    allCoplanar = True
    for i in range (validRidgeIdxs.size):
        #pole pro dva vertexy A a B
        #vrtxA = np.zeros ( (dim + 1 +1 ) )
        #vrtxB = np.zeros ( (dim + 1 +1 ) )
       # vrtxs = []

        rdge = vor.ridge_vertices[validRidgeIdxs[i]]
        #indexy vsech vertexu tvoricich plosny ridge
        for j in range (len(rdge)):
            vrtx = np.zeros ( (dim + 1 +1 +1) ) # vor.ridge_vertices[validRidgeIdxs[i]]  [j]
            #
            for d in range (dim):
                vrtx [d] = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][j]][d]
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
                vrtx [dim +2] = 0 #utilitiesMech.getVertexBC ( vrtx[0:dim])
                vertices_out.append(vrtx)

        #ridges
        ########################################################
        #pole pro ridge nAidx, nBidx, trBc, pocetVertexu, newIdxVertexu
        rdg = np.zeros ( (2 + 0 + 1 + len(vor.ridge_vertices[validRidgeIdxs[i]])  ) )

        #indexy dvou nodu, ktere ridge rozdeluje
        pointA = vor.ridge_points[validRidgeIdxs[i]][0]
        pointB = vor.ridge_points[validRidgeIdxs[i]][1]
        #
        #if(pointA >= node_count):
         #   pointA = -1
        #if(pointB >= node_count):
         #   pointB = -1
        if(pointA >= node_count and pointB<node_count):
            ptA = np.zeros((3))
            ptA[0] = (vor.points[pointB, 0] + vor.points[pointA, 0]  ) /2
            ptA[1] = (vor.points[pointB, 1] + vor.points[pointA, 1]  ) /2
            ptA[2] = (vor.points[pointB, 2] + vor.points[pointA, 2]  ) /2

            pointA = node_count + len(aux_nodes)
            aux_nodes.append(ptA)


        if(pointB >= node_count  and pointA<node_count):
            ptB = np.zeros((3))
            ptB[0] = (vor.points[pointB, 0] + vor.points[pointA, 0]  ) /2
            ptB[1] = (vor.points[pointB, 1] + vor.points[pointA, 1]  ) /2
            ptB[2] = (vor.points[pointB, 2] + vor.points[pointA, 2]  ) /2

            pointB = node_count + len(aux_nodes)
            aux_nodes.append(ptB)

        #
        rdg[0] = pointA
        rdg[1] = pointB

        #print ('%d \t %d'  %(pointA,pointB) )


        #pocet vertexu
        nrVertices = len(vor.ridge_vertices[validRidgeIdxs[i]])
        rdg[2] =  nrVertices
        #print (rdg[3])
        #
        #pridani indexu vertexu
        for v in range ( len(vor.ridge_vertices[validRidgeIdxs[i]]) ):
            rdg[2+1+v] =  verticesIdxDict[ vor.ridge_vertices[validRidgeIdxs[i]][v] ]

            #print (rdg[3+v] )


        #kontrola, ze vsechny body v ridge jsou koplanarni
        for v in range ( nrVertices-3 ):
            pA = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v]][:]
            pB = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+1]][:]
            pC = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+2]][:]
            pD = vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][v+3]][:]
            #pD[2] = 10
            tol = 1e-12
            val = equation_plane(pA, pB, pC, pD)
            if ( val > tol):
                allCoplanar = False
                print('Not coplanar!!! Ridge nr. %d, err: %e' %(i, val ))
            #else: print('Coplanar  %d' %i)

        #normala plochy z prvnich trech vertexu, normalizovana
        planeNormal = getPlaneNormalVector(vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][0]][:],
                                     vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][1]][:],
                                     vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][2]][:])
        #if (printout): print ('planeNormal:')
        #if (printout): print (planeNormal)

        # spojnice mezi generujicimi body
        # mela by byt totozna s normalou plochy. Kdyz ne, tak se pak musi poradi vertexu otocit
        pointNormal = vor.points[pointB] - vor.points[pointA]
        pointNormal /= np.linalg.norm(pointNormal)
        #if (printout): print ('pointNormal')
        #if (printout): print (pointNormal)

        #print ('diff:')
        diff = np.linalg.norm(planeNormal - pointNormal)
        if (diff < 1e-10):
            if (printout): print ('Direction of plane normal OK')
        else:
            if (printout): print ('Direction of plane normal REVERSE')
            #np.flip( (rdg[4:]), axis = 0)
            rdg[nrVertices:] = rdg[:nrVertices-1:-1]

        ##############prerazeni bodu podle uhlu atan2((Vb x Va) . Vn, Va . Vb)##############
        #prumerny bod ze souradnic vertexu
        avgPoint = np.zeros(3)
        for d in range (3):
            for l in range ( nrVertices ):
                avgPoint [d] += vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][l]][d]
            avgPoint[d] /= len(vor.ridge_vertices[validRidgeIdxs[i]])
        #
        #print('avgPoint')
        #print(avgPoint)

        #pole vzajemnych uhlu mezi body a 'stredem'
        angles = np.zeros(nrVertices)

        #vektor, od ktereho se meri uhly (centrum a prvi vertex)
        referenceVector =  vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][0]][:] - avgPoint
        # print ('refvec')
        # print (referenceVector)

        # vypocet uhlu mezi vsemi vektory tvorenymi stredem a vertexy
        # atan2((Vb x Va) . Vn, Va . Vb)
        for l in range ( nrVertices ):
            currVector =  vor.vertices[vor.ridge_vertices[validRidgeIdxs[i]][l]][:] - avgPoint
            angles [l] = np.degrees( np.arctan2(  np.dot( np.cross( referenceVector, currVector ), planeNormal),
                                                np.dot(referenceVector,currVector)   ) )
            if (angles [l] < 0):
                angles [l] = 360 - (-angles [l])

        #if (printout): print (angles)
        #if (printout): print('\n')

        ridges_out.append(rdg)

    if (allCoplanar):
        if (printout):print ('ALL ridges coplanar OK')
        #print ('ALL ridges coplanar OK. Model seems ok.')
    else:
        #if (printout):print ('NOT ALL RIDGES COPLANAR !!!')
        print ('!!! NOT ALL RIDGES COPLANAR !!!')
    print('done.')

    print('Saving nodes...', end ='')
    #writing nodes
    ########################################################################################################
   # headerLine =  "nodeCrdX \t nodeCrdY \t nodeCrdZ \t powRadius \t vorArea \t bcTransX \t bcTransY \t bcTransZ \t bcRotX \t bcRotY \t bcRotZ"
    headerLine =  "Type\tnodeCrdX\tnodeCrdY\tnodeCrdZ\tpowRadius"


    vertIdxStart = node_count + len(aux_nodes)

    fl=open(os.path.join(master_folder,nodesFile),'w')
    np.savetxt(fl,  nodes_out[:,  0:4], delimiter='\t',  fmt = 'Particle\t%.12f\t%.12f\t%.12f\t%.12f',  header = headerLine)
    fl.close()

    #writing auxNodes
    headerLine  = "Type\tnodeCrdX\tnodeCrdY\tnodeCrdZ"
    fmt='AuxNode\t%.12f\t%.12f\t%.12f'
    fl=open(os.path.join(master_folder,auxNodesFile),'w')
    np.savetxt(fl,  aux_nodes, delimiter='\t',   fmt=fmt,  header = headerLine)
    fl.close()
    print('done.')

    print('Saving vertices...', end ='')
    #writing vertices
    ########################################################################################################
    headerLine = 'Type\tvrtxCrdX\tvrtxCrdY\tvrtxCrdZ'

    vertices_print = np.asarray(vertices_out)
    v_count = len(vertices_print)

    fl=open(os.path.join(master_folder,verticesFile),'w')
    np.savetxt(fl, vertices_print[:, 0:3], delimiter='\t', fmt='TrsprtNode%.12f\t%.12f\t%.12f', header = headerLine)
    fl.close()


    ridges_out_trsprt = []
    for i in range (len(ridges_out)):
        ro = np.copy(ridges_out[i])
        ro = np.asarray(ro)
        #
        sh = ro.shape[0]
        #
        #
        for j in range (sh-3):
            #ro[j+3] += 500
            ro[j+3] += node_count + len(aux_nodes)

        nrVert = int( ro[2] )
        #print (nrVert)
        nwarr = np.copy (ro)
        nwarr [0:nrVert] = ro [2: (nrVert+2)]
        nwarr [nrVert+1: (nrVert+3)] = ro [0:2]

        ridges_out_trsprt.append(nwarr)
    print('done.')

    print('Saving TRSPRT elements...', end ='')
    #writing ridges
    ############################################### ridges: idx noduA, idx noduB, transport okr. podminka, idx vertexu
    headerLine = '#ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tverticesIdxs\tMaterial\n'

    fl=open(os.path.join(master_folder,trsprtElemsFile),'w')
    ro = np.asarray(ridges_out_trsprt[0])
    #print (ro)
    fl.write(headerLine)
    #np.savetxt(fl, ro.reshape(1, ro.shape[0]), delimiter='\t', header = headerLine )
    for i in range (len(ridges_out_trsprt)):
        ro = np.asarray(ridges_out_trsprt[i])
        fmt='LTCTRSP\t%d\t%d\t%d'
        sh = ro.shape[0]
        #for j in range (sh-3):
        #    ro[j+3] += node_count

        np.savetxt(fl,  ro.reshape(1, sh) , delimiter='\t', fmt=fmt+'\t%d'*(sh-3)+ '\t1')
        #print (ro)
    fl.close()
    print ('done.')
   # print (ridges_out)


    print('Saving MECH elements...', end ='')
    headerLine = '#ElemType\tnodeAidx\tnodeBidx\tnrOfVertices\tverticesIdxs\tMaterial\n'
    #filtering ridges to ridges with both nodes in sample
    mechElemRidges = []
    for m in range (len(ridges_out)):
        if (ridges_out[m][0] < node_count and ridges_out[m][1] < node_count and ridges_out[m][0] >0  and ridges_out[m][1] > 0):
            mechElemRidges.append( ridges_out[m] )
           # print (mechElemRidges[0])



    fl=open(os.path.join(master_folder,mechElemsFile),'w')
    ro = np.asarray(mechElemRidges[0])
    #print (ro)
    fl.write(headerLine)
    #np.savetxt(fl, ro.reshape(1, ro.shape[0]), delimiter='\t', header = headerLine )
    for i in range (len(mechElemRidges)):
        ro = np.asarray(mechElemRidges[i])
        fmt='LTCBEAM\t%d\t%d\t%d'
        sh = ro.shape[0]
        np.savetxt(fl,  ro.reshape(1, sh) , delimiter='\t', fmt=fmt+'\t%d'*(sh-3)+ '\t0')
        #print (ro)


    print ('done.')
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



class transportIC:
    def __init__(self, vrtxIdx, pressure):
        self.vrtxIdx = vrtxIdx
        self.pressure = pressure

    def getVrtxIdx(self):
        return self.vrtxIdx

    def getPressure(self):
        return self.pressure

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

class mechanicalIC:
    def __init__(self, dim, nodeIdx, mechICArray):
        self.mechICArray = mechICArray
        self.dim = dim
        self.nodeIdx = nodeIdx

    def getNodeIdx(self):
        return self.nodeIdx
    def getMechIC(self):
        return self.mechICArray



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

class mechanicalBC:
    def __init__(self, dim, nodeIdx, mechBCarray):
        self.mechBCarray = mechBCarray
        self.dim = dim
        self.nodeIdx = nodeIdx

    def getDim(self):
        return self.dim
    def getMechBC(self):
        return self.mechBCarray
    def getNodeIdx(self):
        return self.nodeIdx

    def printProps(self):
        print ('Node %d' %(self.nodeIdx))

        if (self.dim ==2):
            print('TransX: %d' %self.mechBCarray[0])
            print('TransY: %d' %self.mechBCarray[1])
            print('RotZ: %d' %self.mechBCarray[2])

        if (self.dim ==3):
            print('TransX: %d' %self.mechBCarray[0])
            print('TransY: %d' %self.mechBCarray[1])
            print('TransZ: %d' %self.mechBCarray[2])
            print('RotX: %d' %self.mechBCarray[3])
            print('RotY: %d' %self.mechBCarray[4])
            print('RotZ: %d' %self.mechBCarray[5])





def saveMechBC(dim, nodes_mechBCmerged):
    print('Saving MECH boundary conditions...', end='')
    sys.stdout.flush()

    #print (len(nodes_mechBCmerged))
    mechBC_out = []

    for i in range (len(nodes_mechBCmerged)):
        if (dim == 2):
            bc = np.zeros (1 + 6)
            #
            bc[0] = nodes_mechBCmerged[i].getNodeIdx()
            #
            bc[1:] = nodes_mechBCmerged[i].getMechBC()
            #
            #print(bc)
        elif (dim == 3):
            bc = np.zeros (1 + 12)
            #
            bc[0] = nodes_mechBCmerged[i].getNodeIdx()
            #
            bc[1:] = nodes_mechBCmerged[i].getMechBC()
            #
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

def saveMasterInput(dim, solver, solStep):
     print('Saving master file...', end='')
     sys.stdout.flush()
     fl=open(os.path.join(master_folder,masterFile),'w')

     fl.write("Dimension\t%d\n"%dim)
     if (solver == 0):
            fl.write('Solver\tSteadyStateLinearSolver\ttime_step\t%e\ttotal_time\t50\n'%solStep)
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

class transportBC:
    def __init__(self,  nodeIdx, transportBCarray):
        self.transportBCarray = transportBCarray
        self.nodeIdx = nodeIdx


    def getTrsprtBC(self):
        return self.transportBCarray
    def getNodeIdx(self):
        return self.nodeIdx


def saveTransportBC(vertices_transportBCmerged, verticesDict, vertIdxStart):
    print('Saving TRSPRT boundary conditions...', end = '')
    sys.stdout.flush()
    trsptBC_out = []

    for i in range (len(vertices_transportBCmerged)):
        bc = np.zeros ((1 + 1+1))
        #
        bc[0] = verticesDict[vertices_transportBCmerged[i].getNodeIdx()] + vertIdxStart
        bc[1:] = vertices_transportBCmerged[i].getTrsprtBC()
        #
        #print (len( bc))


        trsptBC_out.append(bc)

    #
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
    fl.write("TXTNodalExporter pressure 1 pressure")
    fl.close()

    print('done.')
