import numpy as np
import random
import utilitiesGeom
import utilitiesMech

def create2dCantileverUniTens(maxLim, minDist, trials ):
    ### sampling of node point_size
    node_coords,node_mechBC, mechBC_merged  = assemble2DRectangle(maxLim, minDist, trials );


    return node_coords,node_mechBC, mechBC_merged


#
######## FUNCTION FOR CREATING OF A 2D SUPPORTED RECTANGLE MODEL
def assemble2DRectangle (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    node_mechBC = []
    mechBC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    ###############generating of nodes, supported line left vertical ###############
    #mech bc
    lineBC = np.array([0,0,0,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of nodes, a single point top right (a line of zero length) ###############
    lineBC = np.array([-1,-1,-1,0,1,0])

    #defining points of the line
    nodeA = np.array([maxLim[0] - indent, maxLim[1] - indent])
    nodeB = np.array([maxLim[0] - indent, maxLim[1] - indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])

    #rect
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, rectBC)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords, node_mechBC, mechBC_merged




######## FUNCTION FOR CREATING OF A 3D SUPPORTED RECTANGE
def assemble3Dblock():
    indent = 1e-5

    transtBC = np.array([0,1])
    tsptBC = utilitiesGeom.transportBC(20, transtBC)
    transportBC_merged.append(tsptBC)

    ###############generating of points supported line bottom left ###############
    mechBC = np.array([0,0,0,0,0,0,0,0,0,0,0,0])

    nodeA = np.array([indent , indent, indent])
    nodeB = np.array([indent , maxLim[1] - indent, indent])

    oldLen = len(node_coords)
   # utilitiesGeom.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials, True)
    #
    nrOfPoints =  (len(node_coords)) - oldLen

    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of points supported line top right ###############
    mechBC = np.array([-1,-1,-1, -1,-1,-1,    0,1,0, 0,0,0])

    nodeA = np.array([maxLim[0] - indent , indent, maxLim[2] -indent])
    nodeB = np.array([maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials, True)
    #
    nrOfPoints =  (len(node_coords)) - oldLen
     #print (nrOfPoints)

    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')


    ###############generating of points supported surface  left face ###############
    mechBC = np.array([0,0,0,0,0,0,-1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')


    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    #kvadr
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, mechBC)
    newLen = len(node_coords)-1
    nrOfPoints =  (len(node_coords)) - oldLen
   # for n in range ( nrOfPoints ):
       # mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
       # mechBC_merged.append(mBC)
       # print('adding')
    ####################################################################################################
